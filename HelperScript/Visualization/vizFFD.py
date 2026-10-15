#!/usr/bin/env python3
"""
Plot a translucent VTU surface and overlay X,Y,Z points as sphere markers.
Also overlays tangent and normal vectors.

Files expected:
  - surface_deformed.vtu
  - tangent_normals_WING_TST.csv
"""

import numpy as np
import pandas as pd
import pyvista as pv

# -----------------------------
# File paths
# -----------------------------
VTU_FILE = "surface_deformed.vtu"
CSV_FILE = "tangent_normals_WING_TST.csv"

# -----------------------------
# Load VTU surface mesh
# -----------------------------
mesh = pv.read(VTU_FILE)
print(mesh)

# -----------------------------
# Load point data (X, Y, Z)
# -----------------------------
df = pd.read_csv(CSV_FILE)

required_cols = [
    "X", "Y", "Z",
    "Tangent_X", "Tangent_Y", "Tangent_Z",
    "Normal_Y", "Normal_Z"
]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"CSV is missing required columns: {missing}")

points = np.column_stack(
    (df["X"].values,
     df["Y"].values,
     df["Z"].values)
)

# Create PolyData for the points
points_pd = pv.PolyData(points)

# --- NEW: attach tangent and normal vectors ---
points_pd["tangent"] = np.column_stack(
    (df["Tangent_X"].values,
     df["Tangent_Y"].values,
     df["Tangent_Z"].values)
)

points_pd["normal"] = np.column_stack(
    (np.zeros(len(df)),              # Normal_X assumed = 0
     df["Normal_Y"].values,
     df["Normal_Z"].values)
)

# --- NEW: build glyphs ---
TANGENT_SCALE = 1.5
NORMAL_SCALE  = 1.5

tangent_glyphs = points_pd.glyph(
    orient="tangent",
    scale=False,
    factor=TANGENT_SCALE
)

normal_glyphs = points_pd.glyph(
    orient="normal",
    scale=False,
    factor=NORMAL_SCALE
)

# -----------------------------
# Set up plotter
# -----------------------------
plotter = pv.Plotter()

# Add translucent VTU surface
plotter.add_mesh(
    mesh,
    color="lightgray",
    opacity=0.3,
    smooth_shading=True,
    show_edges=False,
    name="surface"
)

# Add points as sphere-like markers
plotter.add_mesh(
    points_pd,
    render_points_as_spheres=True,
    point_size=12,
    color="red",
    name="markers"
)

# --- NEW: add tangent & normal vectors ---
plotter.add_mesh(
    tangent_glyphs,
    color="orange",
    name="tangents"
)

plotter.add_mesh(
    normal_glyphs,
    color="blue",
    name="normals"
)

plotter.add_axes()
plotter.add_legend(
    labels=[
        ("Surface", "lightgray"),
        ("Points", "red"),
        ("Tangents", "orange"),
        ("Normals", "blue"),
    ],
    bcolor="white"
)

# -----------------------------
# Show interactive window
# -----------------------------
plotter.show()