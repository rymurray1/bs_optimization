"""
Equation Flow: Connects grinding → flotation → leaching → solvent extraction
"""
import numpy as np
from grind_break import grind_break
from flotation_recovery import flot_rec, FITTING_PARAMETERS, FLOTATION_CONSTANTS
from leaching import copper_leaching_rate
from solvent_extraction import solvent_extraction


if __name__ == "__main__":
    # Step 1: GRINDING
    print("=" * 80)
    print("STEP 1: GRINDING")
    print("=" * 80)

    # Grinding inputs
    tonnage_input = [100, 150, 50, 10, 3]
    top_size_input = [600, 300, 150, 75, 37.5]
    bottom_size_input = [300, 150, 75, 37.5, 18.75]
    selection_input = [0.8, 0.6, 0.4, 0.2, 0.1]
    break_intensity = 5

    print(f"Feed: {sum(tonnage_input)} tonnes")
    print(f"Break intensity: {break_intensity}")
    print()

    # Run grinding
    ground_tonnage = grind_break(tonnage_input, top_size_input, bottom_size_input, selection_input, break_intensity)

    print("Ground product distribution:")
    for i, (tonnage, top, bottom) in enumerate(zip(ground_tonnage, top_size_input, bottom_size_input)):
        print(f"  Class {i+1} ({bottom:.1f}-{top:.1f} um): {tonnage:.2f} tonnes")
    print(f"Total: {np.sum(ground_tonnage):.2f} tonnes")
    print()

    # Step 2: FLOTATION
    print("=" * 80)
    print("STEP 2: FLOTATION (BY SIZE CLASS)")
    print("=" * 80)

    # Flotation parameters (same for all size classes except particle_size)
    base_optimizable_vars = {
        'contact_angle': 52.3,
        'bubble_z_pot': -0.03,
        'particle_z_pot': -50,
        'sp_power': 1,
        'sp_gas_rate': 2,
        'air_fraction': 0.05,
        'slurry_fraction': 0.15,
        'num_cells': 1,
        'ret_time': 23.67,
        'cell_volume': 700,
        'froth_height': 0.165,
        'frother_conc': 50
    }

    test_inputs = {
        'copper_grade': 28,
        'frother_type': 2,  # MIBC
        'water_or_particle': 'Particle'
    }

    print(f"{'Class':>6} | {'Size (um)':>12} | {'Feed (t)':>10} | {'Recovery':>10} | {'Conc (t)':>10}")
    print("-" * 80)

    total_feed_flot = 0
    total_concentrate = 0

    for i, (tonnage, top, bottom) in enumerate(zip(ground_tonnage, top_size_input, bottom_size_input)):
        # Calculate geometric mean particle size for this class
        particle_size = np.sqrt(top * bottom)

        # Set particle size for this class
        optimizable_vars = base_optimizable_vars.copy()
        optimizable_vars['particle_size'] = particle_size

        # Calculate recovery for this size class
        recovery = flot_rec(
            FITTING_PARAMETERS,
            FLOTATION_CONSTANTS,
            optimizable_vars,
            test_inputs
        )

        # Calculate concentrate tonnage
        concentrate = tonnage * recovery

        total_feed_flot += tonnage
        total_concentrate += concentrate

        print(f"{i+1:6d} | {particle_size:12.2f} | {tonnage:10.2f} | {recovery*100:9.2f}% | {concentrate:10.2f}")

    # Overall mass-weighted recovery
    overall_recovery = total_concentrate / total_feed_flot if total_feed_flot > 0 else 0

    print("-" * 80)
    print(f"{'TOTAL':>6} | {'':>12} | {total_feed_flot:10.2f} | {overall_recovery*100:9.2f}% | {total_concentrate:10.2f}")
    print()

    # Step 3: SUMMARY
    print("=" * 80)
    print("COMBINED GRIND-FLOTATION SUMMARY")
    print("=" * 80)
    print(f"Original feed:        {sum(tonnage_input):8.2f} tonnes @ 28% Cu")
    print(f"After grinding:       {np.sum(ground_tonnage):8.2f} tonnes")
    print(f"After flotation:      {total_concentrate:8.2f} tonnes")
    print(f"Overall recovery:     {overall_recovery*100:8.2f}%")
    print()
    print(f"Tailings:             {total_feed_flot - total_concentrate:8.2f} tonnes")
    print(f"Tailings loss rate:   {(1-overall_recovery)*100:8.2f}%")
    print()

    # Step 3: LEACHING
    print("=" * 80)
    print("STEP 3: COPPER LEACHING")
    print("=" * 80)

    # Leaching parameters
    Ea = 80  # kJ/mol
    T = 190  # °C
    P_O2 = 12  # bar
    n = 0.75
    Phi = 10  # kg/L
    H_plus = 1.0  # -log(pH)
    MFeS2 = 0.1
    k0 = 1.0

    print(f"Concentrate feed to leaching: {total_concentrate:.2f} tonnes @ 28% Cu")
    print(f"Leaching conditions:")
    print(f"  Temperature: {T}°C")
    print(f"  Oxygen pressure: {P_O2} bar")
    print(f"  Activation energy: {Ea} kJ/mol")
    print()

    # Calculate leaching rate at different conversion levels
    print(f"{'Time (hr)':>10} | {'Cu Conv.':>10} | {'Rate (1/hr)':>15} | {'Cu Leached (t)':>18}")
    print("-" * 80)

    # Calculate copper content in concentrate
    cu_content_tonnes = total_concentrate * 0.28  # 28% Cu grade

    time_hours = [0, 4, 8, 16, 24, 48]
    for t in time_hours:
        # Simple approximation of conversion over time
        X_Cu = min(0.95, t / 50)
        rate = copper_leaching_rate(Ea, T, P_O2, n, Phi, H_plus, MFeS2, X_Cu, k0)
        cu_leached = cu_content_tonnes * X_Cu

        print(f"{t:10.1f} | {X_Cu*100:9.2f}% | {rate*3600:15.6e} | {cu_leached:18.2f}")

    print()
    print("FINAL LEACHING SUMMARY:")
    final_conversion = 0.95
    final_cu_leached = cu_content_tonnes * final_conversion
    print(f"  Total Cu in concentrate: {cu_content_tonnes:.2f} tonnes")
    print(f"  Final conversion: {final_conversion*100:.1f}%")
    print(f"  Cu leached: {final_cu_leached:.2f} tonnes")
    print(f"  Cu in residue: {cu_content_tonnes - final_cu_leached:.2f} tonnes")
    print()

    # Step 4: SOLVENT EXTRACTION
    print("=" * 80)
    print("STEP 4: SOLVENT EXTRACTION (SX)")
    print("=" * 80)

    # Convert Cu tonnage to moles for SX calculation
    # Atomic mass of Cu = 63.546 g/mol
    cu_molar_mass = 63.546  # g/mol
    cu_leached_moles = (final_cu_leached * 1e6) / cu_molar_mass  # Convert tonnes to grams then to moles

    print(f"Cu in PLS (pregnant leach solution): {final_cu_leached:.2f} tonnes")
    print(f"Cu in PLS: {cu_leached_moles:.2e} moles")
    print()

    # SX parameters
    K_ex = 100  # Extraction equilibrium constant
    RH = 0.5  # Extractant concentration (M)
    H_plus = 0.01  # H+ concentration (M), pH ~2
    O_A = 1.0  # Organic/aqueous ratio

    print("SX operating conditions:")
    print(f"  Extraction constant (K_ex): {K_ex}")
    print(f"  Extractant concentration [RH]: {RH} M")
    print(f"  H+ concentration: {H_plus} M (pH = {-np.log10(H_plus):.1f})")
    print(f"  Organic/Aqueous ratio (O/A): {O_A}")
    print()

    # Calculate Cu remaining in raffinate (aqueous phase after extraction)
    cu_raffinate_moles = solvent_extraction(cu_leached_moles, K_ex, RH, H_plus, O_A)
    cu_extracted_moles = cu_leached_moles - cu_raffinate_moles

    # Convert back to tonnes
    cu_raffinate_tonnes = (cu_raffinate_moles * cu_molar_mass) / 1e6
    cu_extracted_tonnes = (cu_extracted_moles * cu_molar_mass) / 1e6

    extraction_efficiency = (cu_extracted_moles / cu_leached_moles) * 100

    print("SX RESULTS:")
    print(f"  Cu extracted to organic: {cu_extracted_tonnes:.2f} tonnes ({extraction_efficiency:.2f}%)")
    print(f"  Cu in raffinate (lost):  {cu_raffinate_tonnes:.2f} tonnes ({100-extraction_efficiency:.2f}%)")
    print()

    print("OVERALL PROCESS SUMMARY")
    print("=" * 80)
    original_cu = sum(tonnage_input) * 0.28
    print(f"Original Cu in ore:           {original_cu:.2f} tonnes")
    print(f"After flotation:              {cu_content_tonnes:.2f} tonnes ({cu_content_tonnes/original_cu*100:.1f}%)")
    print(f"After leaching:               {final_cu_leached:.2f} tonnes ({final_cu_leached/original_cu*100:.1f}%)")
    print(f"After solvent extraction:     {cu_extracted_tonnes:.2f} tonnes ({cu_extracted_tonnes/original_cu*100:.1f}%)")
    print()
    print(f"Overall Cu recovery:          {cu_extracted_tonnes/original_cu*100:.2f}%")
