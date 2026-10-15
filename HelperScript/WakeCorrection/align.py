import numpy as np

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