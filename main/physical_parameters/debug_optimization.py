"""
Debug script to understand why optimization is giving low recovery.
"""

import numpy as np
from physical_parameters import flot_rec

# Test with the "optimal" parameters from the optimization
feed_composition = {
    'size_classes': [212.13, 106.07, 38.73, 10.00],
    'grade_classes': [0.05, 0.20, 0.40, 0.75, 1.0],
    'throughput_distribution': [
        [136.22, .87, .42, .20, .29],
        [156.95, .61, .18, .14, .12],
        [309.24, .92, .61, .38, 1.85],
        [386.63, .76, .34, .39, 3.88],
    ],
    'sg_by_grade': [3.6, 4, 4.2, 4.2, 4.5],
    'contact_angles_by_grade': [10, 12.5, 27, 52.3, 60],
    'particle_z_pot_by_grade': [-50, -50, -50, -50, -50],
    'grade_by_grade': [0.002, 8, 12, 28, 34.65]
}

constant_params = {
    'frother': 4,
    'permitivity': 8.854e-12,
    'dielectric': 80,
    'pe': 4,
    'water_or_particle': "Particle",
    'dbl_cell_area': 200,
    'dbl_bbl_ratio': 10
}

# Test case 1: Current "optimal" parameters from output
print("="*80)
print("TEST 1: Current 'Optimal' Parameters (18.96% recovery)")
print("="*80)
sp_power, sp_gas_rate, frother_conc, ret_time, num_cells = 1.3, 2.0, 100, 30, 13
air_fraction, slurry_fraction, bubble_z_pot, cell_volume, froth_height = 0.15, 0.1, -0.0219, 112.1, 0.1

size_classes = feed_composition['size_classes']
grade_classes = feed_composition['grade_classes']
throughput_distribution = feed_composition['throughput_distribution']

total_feed = 0
total_recovered = 0

for size_idx, particle_size in enumerate(size_classes):
    for grade_idx, grade_pct in enumerate(grade_classes):
        contact_angle_deg = feed_composition['contact_angles_by_grade'][grade_idx]
        sg_for_grade = feed_composition['sg_by_grade'][grade_idx]
        particle_z_pot_for_grade = feed_composition['particle_z_pot_by_grade'][grade_idx]
        grade_value = feed_composition['grade_by_grade'][grade_idx]
        mass_in_class = throughput_distribution[size_idx][grade_idx]

        recovery = flot_rec(
            sg_for_grade, sp_power, sp_gas_rate, air_fraction, slurry_fraction,
            particle_z_pot_for_grade, bubble_z_pot, num_cells, ret_time, cell_volume,
            froth_height, constant_params['frother'], frother_conc, particle_size, grade_value,
            contact_angle_deg, constant_params['permitivity'], constant_params['dielectric'],
            constant_params['pe'], constant_params['water_or_particle'],
            constant_params['dbl_cell_area'], constant_params['dbl_bbl_ratio']
        )

        total_feed += mass_in_class
        total_recovered += mass_in_class * recovery

        if mass_in_class > 1.0:  # Only print significant masses
            print(f"Size: {particle_size:7.2f}µm, Grade: {grade_value:6.2f}%, Mass: {mass_in_class:6.2f}t, Recovery: {recovery*100:5.2f}%")

overall_recovery = total_recovered / total_feed
print(f"\nOverall Recovery: {overall_recovery*100:.2f}%")

# Test case 2: Try more aggressive parameters
print("\n" + "="*80)
print("TEST 2: More Aggressive Parameters")
print("="*80)
sp_power, sp_gas_rate, frother_conc, ret_time, num_cells = 1.3, 2.0, 100, 30, 10
air_fraction, slurry_fraction, bubble_z_pot, cell_volume, froth_height = 0.10, 0.25, -0.03, 200, 0.20

total_feed = 0
total_recovered = 0

for size_idx, particle_size in enumerate(size_classes):
    for grade_idx, grade_pct in enumerate(grade_classes):
        contact_angle_deg = feed_composition['contact_angles_by_grade'][grade_idx]
        sg_for_grade = feed_composition['sg_by_grade'][grade_idx]
        particle_z_pot_for_grade = feed_composition['particle_z_pot_by_grade'][grade_idx]
        grade_value = feed_composition['grade_by_grade'][grade_idx]
        mass_in_class = throughput_distribution[size_idx][grade_idx]

        recovery = flot_rec(
            sg_for_grade, sp_power, sp_gas_rate, air_fraction, slurry_fraction,
            particle_z_pot_for_grade, bubble_z_pot, num_cells, ret_time, cell_volume,
            froth_height, constant_params['frother'], frother_conc, particle_size, grade_value,
            contact_angle_deg, constant_params['permitivity'], constant_params['dielectric'],
            constant_params['pe'], constant_params['water_or_particle'],
            constant_params['dbl_cell_area'], constant_params['dbl_bbl_ratio']
        )

        total_feed += mass_in_class
        total_recovered += mass_in_class * recovery

        if mass_in_class > 1.0:
            print(f"Size: {particle_size:7.2f}µm, Grade: {grade_value:6.2f}%, Mass: {mass_in_class:6.2f}t, Recovery: {recovery*100:5.2f}%")

overall_recovery = total_recovered / total_feed
print(f"\nOverall Recovery: {overall_recovery*100:.2f}%")

# Test case 3: Default parameters from physical_parameters.py
print("\n" + "="*80)
print("TEST 3: Default Parameters from physical_parameters.py")
print("="*80)
sp_power, sp_gas_rate, frother_conc, ret_time, num_cells = 1.0, 2.0, 50, 23.67, 10
air_fraction, slurry_fraction, bubble_z_pot, cell_volume, froth_height = 0.05, 0.15, -0.03, 700, 0.165

total_feed = 0
total_recovered = 0

for size_idx, particle_size in enumerate(size_classes):
    for grade_idx, grade_pct in enumerate(grade_classes):
        contact_angle_deg = feed_composition['contact_angles_by_grade'][grade_idx]
        sg_for_grade = feed_composition['sg_by_grade'][grade_idx]
        particle_z_pot_for_grade = feed_composition['particle_z_pot_by_grade'][grade_idx]
        grade_value = feed_composition['grade_by_grade'][grade_idx]
        mass_in_class = throughput_distribution[size_idx][grade_idx]

        recovery = flot_rec(
            sg_for_grade, sp_power, sp_gas_rate, air_fraction, slurry_fraction,
            particle_z_pot_for_grade, bubble_z_pot, num_cells, ret_time, cell_volume,
            froth_height, constant_params['frother'], frother_conc, particle_size, grade_value,
            contact_angle_deg, constant_params['permitivity'], constant_params['dielectric'],
            constant_params['pe'], constant_params['water_or_particle'],
            constant_params['dbl_cell_area'], constant_params['dbl_bbl_ratio']
        )

        total_feed += mass_in_class
        total_recovered += mass_in_class * recovery

        if mass_in_class > 1.0:
            print(f"Size: {particle_size:7.2f}µm, Grade: {grade_value:6.2f}%, Mass: {mass_in_class:6.2f}t, Recovery: {recovery*100:5.2f}%")

overall_recovery = total_recovered / total_feed
print(f"\nOverall Recovery: {overall_recovery*100:.2f}%")
