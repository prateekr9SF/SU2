#!/usr/bin/env python3
"""
Single tapered "box" with a geometry-following lattice.

- Samples many spanwise stations from the surface mesh
- Estimates xLE/xTE/zlo/zhi at each station via robust percentiles
- Applies padding => xmin(y), xmax(y), zmin(y), zmax(y)
- Fits smooth splines for xmin/xmax/zmin/zmax as functions of y
- Builds a lattice with:
    x(i,j,k) = xmin(y_j) + xi_i * (xmax(y_j) - xmin(y_j))
    z(i,j,k) = zmin(y_j) + zeta_k * (zmax(y_j) - zmin(y_j))
  so x follows the underlying geometry with y.

NOTE (SU2): If you only provide 8 corners to SU2, SU2 will create trilinear control points
internally (linear variation in y between end faces). The curved lattice here is ideal for:
  - visualizing how well the lattice fits
  - exporting control points for workflows that can set them explicitly
  - deciding how many spanwise boxes you need if you want SU2 to match it piecewise.
"""

import numpy as np
import pyvista as pv

# Optional: spline (recommended). Falls back to linear interpolation if SciPy not available.
try:
    from scipy.interpolate import UnivariateSpline
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

# -----------------------------
# Inputs
# -----------------------------
VTU_FILE = "surface_deformed.vtu"

BOX_PREFIX = "WING"
BOX_NAME = f"{BOX_PREFIX}_SINGLE"

# Lattice density (uniform in parametric i/j)
DEG_X = 8    # nx = DEG_X + 1
DEG_Y = 15   # ny = DEG_Y + 1
DEG_Z = 1    # nz = DEG_Z + 1 (usually 1 => 2 planes)

# Spanwise sampling resolution for fitting bounds
N_PROFILE_STATIONS = 101  # increase for more fidelity

# Padding relative to local chord and local z-thickness
PAD_LE = 0.03
PAD_TE = 0.05
PAD_Z  = 0.15

# Neighborhood for estimating profiles at each station
WINDOW_FRAC = 1.5

# Robust percentiles
LE_PCTL  = 1.0
TE_PCTL  = 99.0
ZLO_PCTL = 1.0
ZHI_PCTL = 99.0

# Spline smoothing (larger => smoother). Only used if SciPy available.
SPLINE_SMOOTH = 0.0  # 0 => interpolate; try 1e-4..1e-2 if noisy

# Tecplot outputs
WRITE_TEC_CTRL_DAT = True
TEC_CTRL_DAT_FILE = f"{BOX_PREFIX}_single_curved_ctrl_points.dat"

WRITE_TEC_WIREFRAME_DAT = True
TEC_BOX_DAT_FILE = f"{BOX_PREFIX}_single_curved_box_wireframe.dat"

# Visualization
SHOW_MESH_EDGES = False
MESH_OPACITY = 0.60
DRAW_BOUND_WIREFRAME = True
DRAW_CTRL_POINTS = True
CTRL_POINT_SIZE = 8

WRITE_SU2_CTRL_CSV = True
SU2_CTRL_CSV_FILE = "lattice.csv"



# -----------------------------
# Helpers
# -----------------------------
def station_profile(points, yj, half_window):
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    mask = np.abs(y - yj) <= half_window
    if mask.sum() < 50:
        mask = np.abs(y - yj) <= 2.5 * half_window
    if mask.sum() < 10:
        mask = np.ones_like(y, dtype=bool)

    xs = x[mask]
    zs = z[mask]

    xLE = np.percentile(xs, LE_PCTL)
    xTE = np.percentile(xs, TE_PCTL)
    zlo = np.percentile(zs, ZLO_PCTL)
    zhi = np.percentile(zs, ZHI_PCTL)

    xTE = max(xTE, xLE + 1e-12)
    return xLE, xTE, zlo, zhi


def fit_curve(y, v, smooth=0.0):
    """Return a callable v(yq) based on either spline (preferred) or linear interp."""
    y = np.asarray(y, float)
    v = np.asarray(v, float)

    # Ensure strictly increasing y for interpolators
    order = np.argsort(y)
    y = y[order]
    v = v[order]

    # Remove duplicates in y (keep average)
    yu, idx = np.unique(y, return_inverse=True)
    if len(yu) != len(y):
        vv = np.zeros_like(yu)
        cnt = np.zeros_like(yu)
        for i, k in enumerate(idx):
            vv[k] += v[i]
            cnt[k] += 1
        v = vv / np.maximum(cnt, 1.0)
        y = yu

    if HAVE_SCIPY:
        spl = UnivariateSpline(y, v, s=float(smooth))
        return lambda yq: spl(yq)
    else:
        return lambda yq: np.interp(yq, y, v)


def write_tecplot_points(filename, title, xyz):
    xyz = np.asarray(xyz, float)
    n = xyz.shape[0]
    with open(filename, "w") as f:
        f.write(f'TITLE = "{title}"\n')
        f.write('VARIABLES = "X", "Y", "Z"\n')
        f.write(f'ZONE T="{title}", N={n}, E={n}, ZONETYPE=FELINESEG, DATAPACKING=POINT\n')
        for p in xyz:
            f.write(f"{p[0]:.15e} {p[1]:.15e} {p[2]:.15e}\n")
        for i in range(1, n + 1):
            f.write(f"{i} {i}\n")
    print(f"# Wrote Tecplot points: {filename} (N={n})")


def write_tecplot_segments(filename, title, segments):
    pts = []
    conn = []
    idx = 1
    for a, b in segments:
        pts.append(a); pts.append(b)
        conn.append((idx, idx + 1))
        idx += 2
    pts = np.asarray(pts, float)
    with open(filename, "w") as f:
        f.write(f'TITLE = "{title}"\n')
        f.write('VARIABLES = "X", "Y", "Z"\n')
        f.write(f'ZONE T="{title}", N={len(pts)}, E={len(conn)}, ZONETYPE=FELINESEG, DATAPACKING=POINT\n')
        for p in pts:
            f.write(f"{p[0]:.15e} {p[1]:.15e} {p[2]:.15e}\n")
        for i1, i2 in conn:
            f.write(f"{i1} {i2}\n")
    print(f"# Wrote Tecplot segments: {filename} (segments={len(conn)})")


def write_su2_ffd_control_points_csv(csv_path, ctrl):
    """
    Write SU2 control-point lattice to CSV in the explicit indexed format:

        i,j,k,x,y,z

    where ctrl has shape (nx, ny, nz, 3) and
      i = 0..nx-1 corresponds to x-direction (lattice iDegree)
      j = 0..ny-1 corresponds to y-direction (lattice jDegree)
      k = 0..nz-1 corresponds to z-direction (lattice kDegree)

    This matches SU2's Coord_Control_Points[i][j][k][dim] layout.
    """
    ctrl = np.asarray(ctrl, float)
    if ctrl.ndim != 4 or ctrl.shape[-1] != 3:
        raise ValueError(f"ctrl must have shape (nx, ny, nz, 3); got {ctrl.shape}")

    nx, ny, nz, _ = ctrl.shape

    with open(csv_path, "w") as f:
        f.write("i,j,k,x,y,z\n")
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x, y, z = ctrl[i, j, k, :]
                    f.write(f"{i},{j},{k},{x:.15e},{y:.15e},{z:.15e}\n")

    print(f"# Wrote SU2 FFD control points CSV: {csv_path} (nx={nx}, ny={ny}, nz={nz}, rows={nx*ny*nz})")


def main():
    mesh = pv.read(VTU_FILE)
    pts = mesh.points

    # Span range
    y_min, y_max = float(pts[:, 1].min()), float(pts[:, 1].max())
    span = y_max - y_min
    if span <= 0:
        raise RuntimeError("Degenerate span: y_max == y_min")

    # Dense profile stations for fitting bounds
    y_prof = np.linspace(y_min, y_max, N_PROFILE_STATIONS)

    # Neighborhood size (based on avg station spacing in the profile sampling)
    dy_prof = span / (N_PROFILE_STATIONS - 1)
    half_window = 0.5 * WINDOW_FRAC * dy_prof

    # Compute profiles
    xLE = np.zeros_like(y_prof)
    xTE = np.zeros_like(y_prof)
    zlo = np.zeros_like(y_prof)
    zhi = np.zeros_like(y_prof)

    for j, yj in enumerate(y_prof):
        xLE[j], xTE[j], zlo[j], zhi[j] = station_profile(pts, yj, half_window)

    chord = np.maximum(xTE - xLE, 1e-12)
    zth   = np.maximum(zhi - zlo, 1e-12)

    xmin = xLE - PAD_LE * chord
    xmax = xLE + (1.0 + PAD_TE) * chord
    zmin = zlo - PAD_Z * zth
    zmax = zhi + PAD_Z * zth

    # Fit smooth bounds
    fxmin = fit_curve(y_prof, xmin, smooth=SPLINE_SMOOTH)
    fxmax = fit_curve(y_prof, xmax, smooth=SPLINE_SMOOTH)
    fzmin = fit_curve(y_prof, zmin, smooth=SPLINE_SMOOTH)
    fzmax = fit_curve(y_prof, zmax, smooth=SPLINE_SMOOTH)

    # Lattice param coordinates (uniform in i and j)
    nx, ny, nz = DEG_X + 1, DEG_Y + 1, DEG_Z + 1
    xi   = np.linspace(0.0, 1.0, nx)
    eta  = np.linspace(0.0, 1.0, ny)
    zeta = np.linspace(0.0, 1.0, nz)

    # Physical y locations for lattice lines
    y_lat = y_min + eta * (y_max - y_min)

    # Build curved lattice: ctrl[i,j,k,:]
    ctrl = np.zeros((nx, ny, nz, 3), dtype=float)
    for j, yj in enumerate(y_lat):
        x0 = float(fxmin(yj))
        x1 = float(fxmax(yj))
        z0 = float(fzmin(yj))
        z1 = float(fzmax(yj))

        # guard
        x1 = max(x1, x0 + 1e-12)
        z1 = max(z1, z0 + 1e-12)

        for i, s in enumerate(xi):
            xij = x0 + s * (x1 - x0)
            for k, t in enumerate(zeta):
                zij = z0 + t * (z1 - z0)
                ctrl[i, j, k, :] = (xij, yj, zij)

    ctrl_pts = ctrl.reshape(-1, 3)

    if WRITE_SU2_CTRL_CSV:
        write_su2_ffd_control_points_csv(SU2_CTRL_CSV_FILE, ctrl)


    # Build a “wireframe” of the fitted bounds (for visualization):
    #  - 4 curves: (xmin,zmin), (xmax,zmin), (xmax,zmax), (xmin,zmax) sampled along y_prof
    wire = []
    for a, b, c, d in [
        (xmin, zmin, "xmin", "zmin"),
        (xmax, zmin, "xmax", "zmin"),
        (xmax, zmax, "xmax", "zmax"),
        (xmin, zmax, "xmin", "zmax"),
    ]:
        # connect consecutive points along y_prof
        for j in range(len(y_prof) - 1):
            pA = np.array([a[j],     y_prof[j],     b[j]])
            pB = np.array([a[j + 1], y_prof[j + 1], b[j + 1]])
            wire.append((pA, pB))

    # Tecplot
    if WRITE_TEC_CTRL_DAT:
        write_tecplot_points(TEC_CTRL_DAT_FILE, f"{BOX_NAME} Curved Control Points", ctrl_pts)
    if WRITE_TEC_WIREFRAME_DAT:
        write_tecplot_segments(TEC_BOX_DAT_FILE, f"{BOX_NAME} Fitted Bounds Wireframe", wire)

    # Print the *single* SU2 8-corner box (end faces only)
    # (this is the best a single-corner-defined SU2 box can do)
    y0, y1 = y_min, y_max
    p1 = np.array([float(fxmin(y0)), y0, float(fzmin(y0))])
    p2 = np.array([float(fxmax(y0)), y0, float(fzmin(y0))])
    p3 = np.array([float(fxmax(y1)), y1, float(fzmin(y1))])
    p4 = np.array([float(fxmin(y1)), y1, float(fzmin(y1))])
    p5 = np.array([float(fxmin(y0)), y0, float(fzmax(y0))])
    p6 = np.array([float(fxmax(y0)), y0, float(fzmax(y0))])
    p7 = np.array([float(fxmax(y1)), y1, float(fzmax(y1))])
    p8 = np.array([float(fxmin(y1)), y1, float(fzmax(y1))])

    def fmt(p): return f"{p[0]:.10e}, {p[1]:.10e}, {p[2]:.10e}"

    print("\n# --------------------------------------------")
    print("# SU2 SINGLE box (8 corners). Interior CPs in SU2 will be trilinear.")
    print("# Curved lattice written to Tecplot for inspection / external CP workflows.")
    print("# --------------------------------------------\n")
    print(f"FFD_DEGREE = ({DEG_X}, {DEG_Y}, {DEG_Z})")
    print("FFD_DEFINITION = \\")
    print(
        "  "
        f"({BOX_NAME}, "
        f"{fmt(p1)}, {fmt(p2)}, {fmt(p3)}, {fmt(p4)}, "
        f"{fmt(p5)}, {fmt(p6)}, {fmt(p7)}, {fmt(p8)})\n"
    )

    # Visualize
    pl = pv.Plotter()
    pl.add_mesh(mesh, opacity=MESH_OPACITY, show_edges=SHOW_MESH_EDGES)

    if DRAW_BOUND_WIREFRAME:
        for a, b in wire:
            pl.add_lines(np.vstack([a, b]), color="black", width=2)

    if DRAW_CTRL_POINTS:
        pl.add_points(ctrl_pts, render_points_as_spheres=True, point_size=CTRL_POINT_SIZE, color="red")

    pl.add_axes()
    pl.show_grid()
    pl.show()


if __name__ == "__main__":
    main()
