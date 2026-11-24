import numpy as np

"""
Generate the coordinate and velocity transformation equations for rotating
a CFD solution from the aircraft body-fixed coordinate system to the
wind (freestream) coordinate system.

The script accepts a user-defined angle of attack (AoA), computes the
corresponding rotation matrix about the body Y-axis, and prints the
transformation equations for both mesh coordinates and velocity components.
The output is formatted for direct use in Tecplot 360 equation definitions.

Coordinate transformation:
    X_fs =  cos(alpha) * X + sin(alpha) * Z
    Y_fs =  Y
    Z_fs = -sin(alpha) * X + cos(alpha) * Z

Velocity transformation:
    U_fs =  cos(alpha) * Velocity_X + sin(alpha) * Velocity_Z
    V_fs =  Velocity_Y
    W_fs = -sin(alpha) * Velocity_X + cos(alpha) * Velocity_Z

This transformation is commonly used to express CFD results in the
freestream reference frame for wake visualization and post-processing,
ensuring that wake convection is aligned with the aircraft flight direction.
"""

# User input: angle of attack in degrees
aoa_deg = 4

alpha = np.deg2rad(aoa_deg)
c = np.cos(alpha)
s = np.sin(alpha)

print(f"# AoA = {aoa_deg:.6f} deg")
print(f"# cos(alpha) = {c:.9f}")
print(f"# sin(alpha) = {s:.9f}")
print()

print("# Coordinate transformation: body axes -> wind/freestream axes")
print(f"{{X_fs}} = {c:.9f}*{{X}} + {s:.9f}*{{Z}}")
print(f"{{Y_fs}} = {{Y}}")
print(f"{{Z_fs}} = -{s:.9f}*{{X}} + {c:.9f}*{{Z}}")
print()

print("# Velocity transformation: body axes -> wind/freestream axes")
print(f"{{U_fs}} = {c:.9f}*{{Velocity_X}} + {s:.9f}*{{Velocity_Z}}")
print(f"{{V_fs}} = {{Velocity_Y}}")
print(f"{{W_fs}} = -{s:.9f}*{{Velocity_X}} + {c:.9f}*{{Velocity_Z}}")