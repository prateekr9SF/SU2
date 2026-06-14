#!/usr/bin/env python3

import re

g = 9.80665
filename = "forces_breakdown.dat"

Sref = 229.445
Mach_target = 0.85
gamma = 1.4
R = 287.058

with open(filename, "r") as f:
    text = f.read()


def extract(pattern, text, name):
    match = re.search(pattern, text, re.IGNORECASE)

    if match is None:
        print(f"Warning: Could not find {name}")
        return None

    value = match.group(1).strip().rstrip(".,;:")
    return float(value)


force_factor = extract(
    r"non-dimensional factor:\s*([0-9Ee+\-\.]+)",
    text,
    "force multiplier"
)

CL = extract(
    r"Total\s+CL:\s*([0-9Ee+\-\.]+)",
    text,
    "CL"
)

Mach = extract(
    r"Mach\s+number:\s*([0-9Ee+\-\.]+)",
    text,
    "Mach number"
)

Re = extract(
    r"Reynolds\s+number:\s*([0-9Ee+\-\.]+)",
    text,
    "Reynolds number"
)

rho = extract(
    r"Density:\s*([0-9Ee+\-\.]+)",
    text,
    "Density"
)

velocity = extract(
    r"Velocity:\s*([0-9Ee+\-\.]+)",
    text,
    "Velocity"
)

qinf = extract(
    r"Dynamic pressure:\s*([0-9Ee+\-\.]+)",
    text,
    "Dynamic pressure"
)

if CL is None:
    raise RuntimeError("Could not find CL")

if force_factor is None:
    raise RuntimeError("Could not find force multiplier")

lift = CL * force_factor
mass = lift / g

print()
print("=" * 70)
print("SU2 FORCE BREAKDOWN ANALYSIS")
print("=" * 70)

if Mach is not None:
    print(f"Mach Number           : {Mach:.4f}")

if Re is not None:
    print(f"Reynolds Number       : {Re:.3e}")

if rho is not None:
    print(f"Density (kg/m^3)      : {rho:.6f}")

if velocity is not None:
    print(f"Velocity (m/s)        : {velocity:.3f}")

if qinf is not None:
    print(f"Dynamic Pressure (Pa) : {qinf:.3f}")

print()
print(f"Lift Coefficient CL   : {CL:.6f}")
print(f"Force Multiplier (N)  : {force_factor:,.3f}")

print()
print("-" * 70)
print("DIMENSIONAL RESULTS")
print("-" * 70)
print(f"Lift Force (N)        : {lift:,.3f}")
print(f"Aircraft Weight (N)   : {lift:,.3f}")
print(f"Aircraft Mass (kg)    : {mass:,.3f}")

flight_levels = {
    "FL350": {"rho": 0.3796, "T": 218.81},
    "FL450": {"rho": 0.2371, "T": 216.65},
    "FL550": {"rho": 0.1416, "T": 216.65},
}

print()
print("-" * 70)
print("CL REQUIRED FOR SAME AIRCRAFT MASS")
print("-" * 70)

for fl, atm in flight_levels.items():
    rho_fl = atm["rho"]
    T_fl = atm["T"]

    a_fl = (gamma * R * T_fl) ** 0.5
    V_fl = Mach_target * a_fl
    q_fl = 0.5 * rho_fl * V_fl**2

    CL_required = lift / (q_fl * Sref)

    print(f"{fl}")
    print(f"  Temperature (K)     : {T_fl:.2f}")
    print(f"  Density (kg/m^3)    : {rho_fl:.6f}")
    print(f"  Velocity (m/s)      : {V_fl:.3f}")
    print(f"  Dynamic Pressure Pa : {q_fl:.3f}")
    print(f"  Required CL         : {CL_required:.6f}")
    print()

print("-" * 70)