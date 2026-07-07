#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute aircraft CG locations corresponding to selected
percentages of the Mean Aerodynamic Chord (MAC).

Assumptions:
- MAC_X, MAC_Y, MAC_Z correspond to the leading edge
  of the Mean Aerodynamic Chord (LEMAC).
- Aircraft x-axis points downstream.
"""

import pandas as pd

# ==========================================================
# User Inputs
# ==========================================================

MAC   = 3.19     # Mean Aerodynamic Chord [m]
MAC_X = 19.64727    # LEMAC x-coordinate [m]
MAC_Y = 7.72524     # LEMAC y-coordinate [m]
MAC_Z = 0.00366    # LEMAC z-coordinate [m]

# CG locations to evaluate (% MAC)
cg_percentages = [10, 15, 25, 35, 45]

# ==========================================================
# Compute CG Coordinates
# ==========================================================

results = []

for pct in cg_percentages:

    x_cg = MAC_X + (pct / 100.0) * MAC
    y_cg = MAC_Y
    z_cg = MAC_Z

    results.append(
        {
            "% MAC": pct,
            "X_CG (m)": x_cg,
            "Y_CG (m)": y_cg,
            "Z_CG (m)": z_cg,
        }
    )

df = pd.DataFrame(results)

# ==========================================================
# Output
# ==========================================================

print("\nAircraft CG Locations")
print("=" * 70)
print(f"MAC    = {MAC:.4f} m")
print(f"LEMAC  = ({MAC_X:.4f}, {MAC_Y:.4f}, {MAC_Z:.4f}) m")
print("=" * 70)

print(
    df.to_string(
        index=False,
        formatters={
            "X_CG (m)": "{:.4f}".format,
            "Y_CG (m)": "{:.4f}".format,
            "Z_CG (m)": "{:.4f}".format,
        },
    )
)

print("=" * 70)