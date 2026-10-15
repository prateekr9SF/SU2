#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# -----------------------------
# Surface CSV files (0 → 6)
# -----------------------------
surface_csvs = [f"surface_coords_{i}.csv" for i in range(7)]
surface_labels = [f"MPI rank {i}" for i in range(7)]

# Marker
marker_list = ["o"]
surface_markers = [marker_list[i % len(marker_list)] for i in range(len(surface_csvs))]

# -----------------------------
# Create figure
# -----------------------------
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# We DO NOT rely on ax.dist (ignored in many Matplotlib versions)
# ax.dist = 7

# -----------------------------
# Colormap for MPI ranks
# -----------------------------
cmap = mpl.colormaps.get_cmap("Dark2").resampled(len(surface_csvs))

# -----------------------------
# Plot scatter for all ranks
# -----------------------------
for i, (csv, marker) in enumerate(zip(surface_csvs, surface_markers)):
    df = pd.read_csv(csv)
    ax.scatter(
        df["X"], df["Y"], df["Z"],
        s=5,
        marker=marker,
        alpha=0.5,
        color=cmap(i),
    )

# -----------------------------
# Equal aspect ratio + zoom
# -----------------------------
xlim = ax.get_xlim()
ylim = ax.get_ylim()
zlim = ax.get_zlim()

max_range = max(
    xlim[1] - xlim[0],
    ylim[1] - ylim[0],
    zlim[1] - zlim[0],
)

x_mid = 0.5 * (xlim[0] + xlim[1])
y_mid = 0.5 * (ylim[0] + ylim[1])
z_mid = 0.5 * (zlim[0] + zlim[1])

# ---- Set your zoom level here ----
# <1 zooms in, >1 zooms out
zoom_factor = 0.1   # try 0.3 for very strong zoom-in
max_range *= zoom_factor
# ----------------------------------

#ax.set_xlim(x_mid - max_range / 2, x_mid + max_range / 2)
#ax.set_ylim(y_mid - max_range / 2, y_mid + max_range / 2)
#ax.set_zlim(z_mid - max_range / 2, z_mid + max_range / 2)

# -----------------------------
# Remove axes completely
# -----------------------------
ax.set_axis_off()

# -----------------------------
# Colorbar for MPI ranks (BOTTOM)
# -----------------------------
norm = Normalize(vmin=0, vmax=len(surface_csvs) - 1)
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

cbar = plt.colorbar(
    sm,
    ax=ax,
    orientation='horizontal',
    fraction=0.06,
    pad=0.10
)

cbar.set_label("MPI Rank", fontsize=18, family="Times New Roman")

cbar.set_ticks(np.arange(len(surface_csvs)))
cbar.set_ticklabels([str(i) for i in range(len(surface_csvs))])

for t in cbar.ax.get_xticklabels():
    t.set_fontsize(14)
    t.set_family("Times New Roman")

# -----------------------------
# Figure adjustments
# -----------------------------
F = plt.gcf()
Size = F.get_size_inches()
F.set_size_inches(Size[0] * 1.5, Size[1] * 1.5, forward=True)

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

plt.tight_layout()
plt.show()