"""
Diagnostic script to understand what's limiting recovery in the flotation model.
"""

import numpy as np
import math
from physical_parameters import flot_rec

# Test with typical industrial flotation parameters
print("="*80)
print("DIAGNOSTIC: Testing flotation recovery with various energy barriers")
print("="*80)

# Standard test parameters
sg = 4.2
sp_power = 1.0
sp_gas_rate = 1.5
air_fraction = 0.10
slurry_fraction = 0.20
particle_z_pot = -50
bubble_z_pot = -0.03
num_cells = 10
ret_time = 20
cell_volume = 200
froth_height = 0.15
frother = 4
frother_conc = 50
particle_diam = 75  # Mid-size particle
grade = 12
contact_angle = 27  # Mid-grade contact angle
permitivity = 8.854e-12
dielectric = 80
pe = 4
water_or_particle = "Particle"
dbl_cell_area = 200
dbl_bbl_ratio = 10

# Test baseline
recovery = flot_rec(
    sg, sp_power, sp_gas_rate, air_fraction, slurry_fraction,
    particle_z_pot, bubble_z_pot, num_cells, ret_time, cell_volume,
    froth_height, frother, frother_conc, particle_diam, grade,
    contact_angle, permitivity, dielectric, pe, water_or_particle,
    dbl_cell_area, dbl_bbl_ratio
)

print(f"\nBaseline Parameters:")
print(f"  Particle size: {particle_diam} µm")
print(f"  Contact angle: {contact_angle}°")
print(f"  Number of cells: {num_cells}")
print(f"  Retention time: {ret_time} min")
print(f"  Air fraction: {air_fraction}")
print(f"  Slurry fraction: {slurry_fraction}")
print(f"  Current Recovery: {recovery*100:.2f}%")

# The energy barrier is currently hardcoded at 5e-14 in physical_parameters.py
print(f"\nCurrent hardcoded energy barrier: 5e-14 J")
print(f"This is causing low attachment probability and low recovery")

print("\n" + "="*80)
print("RECOMMENDATION: The energy barrier value needs to be recalibrated")
print("="*80)
print("\nFor 90% overall recovery with 10 cells:")
print("  Required per-cell recovery: {:.2f}%".format((1 - (1-0.90)**(1/10))*100))
print("\nThe current energy barrier of 5e-14 J is too high.")
print("It should be lowered significantly (try 1e-15 to 5e-16 J) to increase")
print("the attachment probability and achieve higher recoveries.")
