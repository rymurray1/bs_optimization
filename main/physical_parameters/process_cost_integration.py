"""
Process-to-Cost Integration Module - Grinding Stage

Starting with grinding to establish the pattern for connecting
physical operations to TEA costs.
"""

import numpy as np
from tea import (
    feed_size, product_size, p1, g,
    work_index, work,
    jaw_crusher_small, jaw_crusher_large,
    motor_small, motor_large,
    hours_per_op_year, electricity_cost,
    lang_factor, inflation_factor,
    equipment_cost,
    flotation_work, water_flow_rate_multiplier,
    residence_time_flotation, cfa_density, froth_factor_excess,
    seconds_per_op_year, residence_time_flotation_sec,
    water_cost, flotation_volume_required, flotation_rates,
    sl_ratio, acid_concentration, nitric_acid_density,
    residence_time_leaching, leaching_volume_required, leaching_rates
)
from grind_break import grind_break
from flotation_recovery import flot_rec
from leaching import leaching_copper_recovery, calculate_acid_consumption
from solvent_extraction import solvent_extraction
from electrowinning import cell_voltage


class GrindingStage:
    """
    Grinding stage: Connects particle size reduction to costs.

    Physical Parameters → Cost Calculation:
    - feed_size, product_size → work (kWh/t) → OPEX (electricity)
    - work × throughput → required_power_kw → CAPEX (equipment)
    """

    def __init__(self, feed_size_um=10000, product_size_um=100, p1_um=1000, g_val=2.5):
        """
        Args:
            feed_size_um: Feed particle size (microns)
            product_size_um: Product particle size (microns)
            p1_um: Reference size (microns)
            g_val: g value (g/rev)
        """
        self.feed_size = feed_size_um
        self.product_size = product_size_um
        self.p1 = p1_um
        self.g = g_val

        # Calculate work index and specific work using Bond's equation
        self.work_index = self._calculate_work_index()
        self.work = self._calculate_work()

    def _calculate_work_index(self):
        """Bond work index (kWh/t)."""
        work_index = 13
        return work_index
    
    def _calculate_work(self):
        """Specific grinding work (kWh/t)."""
        product_size_power = self.product_size ** 0.5
        feed_size_power = self.feed_size ** 0.5

        return 10 * self.work_index * ((1 / product_size_power) - (1 / feed_size_power))

    def calculate_capex(self, throughput_tpa):
        """
        CAPEX = equipment cost for jaw crusher + motor

        Steps:
        1. Calculate required power from throughput and work
        2. Select equipment basis (small vs large)
        3. Scale cost using power-law
        4. Apply inflation and lang factors
        """
        # Step 1: Required power
        required_power_kw = throughput_tpa * self.work / hours_per_op_year

        # Step 2: Select basis
        if required_power_kw <= jaw_crusher_small['power']:
            basis = 1
        else:
            basis = 2

        # Step 3: Equipment costs with power-law scaling
        crusher_cost = equipment_cost('jaw_crusher', required_power_kw, basis)
        motor_cost = equipment_cost('motor', required_power_kw, basis if basis == 1 else 2)

        # Step 4: Apply factors
        total_equipment_cost = (crusher_cost + motor_cost) * inflation_factor
        installed_cost = total_equipment_cost * lang_factor

        return installed_cost

    def calculate_opex(self, throughput_tpa):
        """
        OPEX = annual electricity cost

        Steps:
        1. Calculate annual energy from throughput and work
        2. Convert to cost using electricity price
        """
        # Annual energy consumption
        annual_power_kwh = throughput_tpa * self.work
        annual_power_mwh = annual_power_kwh / 1000

        # Annual cost
        return annual_power_mwh * electricity_cost

    def process(self, tonnage_input, top_size_input, bottom_size_input, selection_input, break_intensity):
        """
        Calculate grinding product using population balance model.

        Uses grind_break from equation_flow_official.py to calculate particle
        size distribution after grinding.

        Args:
            tonnage_input: Feed tonnage array by size class
            top_size_input: Top size for each size class (microns)
            bottom_size_input: Bottom size for each size class (microns)
            selection_input: Selection function values
            break_intensity: Breakage intensity (number of cycles)

        Returns:
            dict with keys:
                - 'ground_tonnage': Product tonnage distribution
                - 'total_feed': Total feed tonnage
                - 'total_product': Total product tonnage
                - 'p80_size': P80 particle size (microns)
        """
        # Run grind_break equation
        ground_tonnage = grind_break(tonnage_input, top_size_input, bottom_size_input,
                                     selection_input, break_intensity)

        total_feed = sum(tonnage_input)
        total_product = np.sum(ground_tonnage)

        # Calculate P80 particle size (80% passing size)
        cumulative_mass = 0
        p80_size = 0
        for tonnage, bottom, top in zip(ground_tonnage, bottom_size_input, top_size_input):
            cumulative_mass += tonnage
            if cumulative_mass >= 0.8 * total_product:
                p80_size = (top + bottom) / 2
                break

        return {
            'ground_tonnage': ground_tonnage,
            'total_feed': total_feed,
            'total_product': total_product,
            'p80_size': p80_size,
            'mass_balance_ok': abs(total_feed - total_product) < 0.01
        }


class FlotationStage:
    """
    Flotation stage: Separates valuable minerals using bubble attachment.

    Physical Parameters → Cost Calculation:
    - throughput → tank volume → CAPEX (vertical agitated tank)
    - throughput → power + water → OPEX (electricity + water)
    """

    def __init__(self, flotation_work_kwh_per_t=20):
        """
        Args:
            flotation_work_kwh_per_t: Specific flotation work (kWh/t)
        """
        self.flotation_work = flotation_work_kwh_per_t

    def calculate_capex(self, throughput_tpa):
        """
        CAPEX = vertical agitated tank cost

        Steps:
        1. Calculate required tank volume from throughput
        2. Scale cost using power-law
        3. Apply inflation and lang factors
        """
        # Step 1: Calculate required volume
        volume_required = flotation_volume_required(throughput_tpa)

        # Step 2: Equipment cost with power-law scaling
        tank_cost = equipment_cost('vertical_tank', volume_required)

        # Step 3: Apply factors
        installed_cost = tank_cost * inflation_factor * lang_factor

        return installed_cost

    def calculate_opex(self, throughput_tpa):
        """
        OPEX = annual power cost + water cost

        Steps:
        1. Calculate annual power consumption
        2. Calculate annual water consumption
        3. Sum costs
        """
        # Power cost
        annual_power_kwh = throughput_tpa * self.flotation_work
        annual_power_mwh = annual_power_kwh / 1000
        power_cost = annual_power_mwh * electricity_cost

        # Water cost
        rates = flotation_rates(throughput_tpa)
        annual_water_m3 = rates['water_rate'] * seconds_per_op_year
        water_cost_annual = annual_water_m3 * water_cost

        return {
            'power': power_cost,
            'water': water_cost_annual,
            'total': power_cost + water_cost_annual
        }

    def process(self, fitting_parameters, flotation_constants, optimizable_vars, inputs):
        """
        Calculate flotation recovery using physics-based model.

        Uses flot_rec from equation_flow_official.py to calculate recovery.
        This determines how much concentrate goes to downstream stages.

        Args:
            fitting_parameters: Dict with keys: b, alpha, coverage, bubble_f, detach_f, bulk_zone
            flotation_constants: Dict with keys: specific_gravity, grade, permitivity,
                                dielectric, pe, cell_area, bbl_ratio
            optimizable_vars: Dict with keys: particle_size, contact_angle, bubble_z_pot,
                             particle_z_pot, sp_power, sp_gas_rate, air_fraction, slurry_fraction,
                             num_cells, ret_time, cell_volume, froth_height, frother_conc
            inputs: Dict with keys: frother_type, water_or_particle

        Returns:
            dict with keys:
                - 'recovery': Flotation recovery (0-1)
                - 'feed_tonnage': Feed tonnage to flotation
                - 'concentrate_tonnage': Concentrate tonnage (feed × recovery)
                - 'tails_tonnage': Tailings tonnage (feed × (1-recovery))
                - 'particle_size': Particle size used (microns)
                - 'contact_angle': Contact angle used (degrees)
        """
        # Calculate recovery using physics-based model
        recovery = flot_rec(fitting_parameters, flotation_constants, optimizable_vars, inputs)

        # Get feed tonnage from flotation constants (if provided)
        feed_tonnage = flotation_constants.get('feed_tonnage', 0)

        # Calculate mass balance
        concentrate_tonnage = feed_tonnage * recovery
        tails_tonnage = feed_tonnage * (1 - recovery)

        return {
            'recovery': recovery,
            'feed_tonnage': feed_tonnage,
            'concentrate_tonnage': concentrate_tonnage,
            'tails_tonnage': tails_tonnage,
            'particle_size': optimizable_vars['particle_size'],
            'contact_angle': optimizable_vars['contact_angle']
        }


class LeachingStage:
    """
    Leaching stage: Dissolves copper from concentrate using acid.

    Physical Parameters → Cost Calculation:
    - concentrate_tonnage → tank volume → CAPEX (atmospheric tank)
    - copper_leached → acid consumption → OPEX (acid cost)

    Key dependency: Uses concentrate tonnage from flotation recovery.
    """

    def __init__(self, acid_price_per_tonne=100, acid_recovery_fraction=0.92):
        """
        Args:
            acid_price_per_tonne: Sulfuric acid price ($/tonne)
                                 Industrial H2SO4: $50-150/tonne
                                 Concentrated HNO3: $300-500/tonne
            acid_recovery_fraction: Fraction of acid recovered and recycled (0-1)
                                   Typical industrial: 0.90-0.95
                                   Net consumption = stoichiometric × (1 - recovery)
        """
        self.acid_price = acid_price_per_tonne
        self.acid_recovery = acid_recovery_fraction

    def calculate_capex(self, concentrate_tpa):
        """
        CAPEX = atmospheric tank cost

        Steps:
        1. Calculate required tank volume from concentrate throughput
        2. Scale cost using power-law
        3. Apply inflation and lang factors

        Args:
            concentrate_tpa: Concentrate throughput (NOT ore feed!)
        """
        # Calculate required volume
        volume_required = leaching_volume_required(concentrate_tpa)

        # Equipment cost with power-law scaling
        tank_cost = equipment_cost('atm_tank', volume_required)

        # Apply factors
        installed_cost = tank_cost * inflation_factor * lang_factor

        return installed_cost

    def calculate_opex(self, copper_leached_kg):
        """
        OPEX = annual acid cost (accounting for acid recovery/recycling)

        Steps:
        1. Calculate stoichiometric acid consumption
        2. Apply acid recovery factor to get net consumption
        3. Convert to annual cost

        Args:
            copper_leached_kg: Copper leached per year (kg)

        Returns:
            dict with acid cost breakdown
        """
        # Calculate stoichiometric acid consumption
        acid_result = calculate_acid_consumption(copper_leached_kg)

        # Stoichiometric acid requirement (convert kg to tonnes)
        stoichiometric_acid_kg = acid_result['acid_consumption_kg']
        stoichiometric_acid_tonnes = stoichiometric_acid_kg / 1000

        # Net consumption after recovery (only makeup acid needed)
        net_acid_tonnes = stoichiometric_acid_tonnes * (1 - self.acid_recovery)

        # Cost (only pay for makeup acid, not recycled acid)
        acid_cost = net_acid_tonnes * self.acid_price

        return {
            'acid_stoichiometric_tpa': stoichiometric_acid_tonnes,
            'acid_net_consumption_tpa': net_acid_tonnes,
            'acid_recovery_fraction': self.acid_recovery,
            'acid_per_kg_cu': acid_result['acid_per_kg_cu'],
            'acid_cost': acid_cost,
            'total': acid_cost
        }

    def process(self, concentrate_tpa, copper_in_concentrate_tpa, Ea=80, T=190,
                P_O2=12, n=0.75, H_plus=0.2, MFeS2=0.5, leach_time=12):
        """
        Calculate leaching recovery using kinetic model.

        Uses leaching_copper_recovery from equation_flow_official.py.

        Args:
            concentrate_tpa: Concentrate feed (tonnes/year)
            copper_in_concentrate_tpa: Copper in concentrate (tonnes/year)
            Ea: Activation energy (kJ/mol)
            T: Temperature (C)
            P_O2: Oxygen partial pressure (bar)
            n: Reaction order
            H_plus: H+ concentration
            MFeS2: FeS2 mass fraction
            leach_time: Leaching time (hours)

        Returns:
            dict with keys:
                - 'leach_recovery': Leaching recovery (0-1)
                - 'copper_leached_tpa': Copper leached (tonnes/year)
                - 'copper_leached_kg': Copper leached (kg/year)
                - 'solution': Leaching solution data
        """
        # Run leaching kinetics model
        X_Cu = 0  # Initial conversion
        leach_recovery, solution = leaching_copper_recovery(
            Ea, T, P_O2, n, H_plus, MFeS2, X_Cu, leach_time=leach_time
        )

        # Calculate copper leached
        copper_leached_tpa = copper_in_concentrate_tpa * leach_recovery
        copper_leached_kg = copper_leached_tpa * 1000

        return {
            'leach_recovery': leach_recovery,
            'concentrate_feed_tpa': concentrate_tpa,
            'copper_in_concentrate_tpa': copper_in_concentrate_tpa,
            'copper_leached_tpa': copper_leached_tpa,
            'copper_leached_kg': copper_leached_kg,
            'solution': solution
        }


class SolventExtractionStage:
    """
    Solvent Extraction (SX) stage: Purifies copper from leach solution using organic extractants.

    Physical Parameters → Cost Calculation:
    - solution_volume → mixer-settler volume → CAPEX (tank-like equipment)
    - copper_extracted → organic solvent consumption → OPEX (solvent makeup cost)

    Key dependency: Uses copper solution from leaching stage.
    """

    def __init__(self, solvent_price_per_litre=5, solvent_losses_fraction=0.01):
        """
        Args:
            solvent_price_per_litre: Organic extractant price ($/L)
                                    LIX reagents: $3-8/L
            solvent_losses_fraction: Fraction of solvent lost per cycle (0-1)
                                    Typical: 0.005-0.02 (0.5-2% losses)
        """
        self.solvent_price = solvent_price_per_litre
        self.solvent_losses = solvent_losses_fraction

    def calculate_capex(self, copper_leached_tpa):
        """
        CAPEX = mixer-settler equipment cost

        Uses atmospheric tank cost basis as proxy for mixer-settlers.

        Steps:
        1. Estimate volume from copper throughput
        2. Scale cost using power-law
        3. Apply inflation and lang factors

        Args:
            copper_leached_tpa: Copper from leaching (tonnes/year)
        """
        # Estimate volume needed (rough estimate: 0.1 m3 per tpa Cu)
        # This is a placeholder - ideally would calculate from residence time
        volume_required = copper_leached_tpa * 0.1

        # Use atmospheric tank basis as proxy
        tank_cost = equipment_cost('atm_tank', volume_required)

        # Apply factors
        installed_cost = tank_cost * inflation_factor * lang_factor

        return installed_cost

    def calculate_opex(self, copper_extracted_kg, solution_volume_m3):
        """
        OPEX = organic solvent makeup cost

        Steps:
        1. Calculate solvent inventory from solution volume
        2. Calculate annual solvent losses
        3. Convert to cost

        Args:
            copper_extracted_kg: Copper extracted per year (kg)
            solution_volume_m3: Solution volume processed (m3/year)

        Returns:
            dict with solvent cost breakdown
        """
        # Typical O/A ratio is 1:1, so organic volume ≈ aqueous volume
        organic_volume_m3 = solution_volume_m3
        organic_volume_litres = organic_volume_m3 * 1000

        # Annual solvent makeup (only replace losses)
        annual_solvent_losses_litres = organic_volume_litres * self.solvent_losses

        # Cost
        solvent_cost = annual_solvent_losses_litres * self.solvent_price

        return {
            'organic_volume_litres': organic_volume_litres,
            'solvent_losses_litres': annual_solvent_losses_litres,
            'solvent_losses_fraction': self.solvent_losses,
            'solvent_cost': solvent_cost,
            'total': solvent_cost
        }

    def process(self, copper_leached_kg, K_ex=10, RH=0.2, O_A=1.2, pH=1.1):
        """
        Calculate solvent extraction recovery using equilibrium model.

        Uses solvent_extraction from solvent_extraction.py.

        Args:
            copper_leached_kg: Copper from leaching (kg/year)
            K_ex: Extraction equilibrium constant (dimensionless)
            RH: Extractant concentration (molar)
            O_A: Organic/aqueous ratio (dimensionless)
            pH: Solution pH

        Returns:
            dict with keys:
                - 'sx_recovery': SX recovery (0-1)
                - 'copper_feed_kg': Copper feed from leaching (kg/year)
                - 'copper_extracted_kg': Copper extracted to organic (kg/year)
                - 'copper_raffinate_kg': Copper remaining in raffinate (kg/year)
                - 'copper_feed_moles': Copper feed (moles)
                - 'copper_extracted_moles': Copper extracted (moles)
        """
        # Convert copper to moles
        MW_Cu = 63.546  # g/mol
        copper_feed_moles = copper_leached_kg * 1000 / MW_Cu  # kg to g to mol

        # Run solvent extraction model
        copper_raffinate_moles = solvent_extraction(copper_feed_moles, K_ex, RH, O_A, pH)
        copper_extracted_moles = copper_feed_moles - copper_raffinate_moles

        # Calculate recovery
        sx_recovery = copper_extracted_moles / copper_feed_moles if copper_feed_moles > 0 else 0

        # Convert back to kg
        copper_extracted_kg = copper_extracted_moles * MW_Cu / 1000
        copper_raffinate_kg = copper_raffinate_moles * MW_Cu / 1000

        return {
            'sx_recovery': sx_recovery,
            'copper_feed_kg': copper_leached_kg,
            'copper_extracted_kg': copper_extracted_kg,
            'copper_raffinate_kg': copper_raffinate_kg,
            'copper_feed_moles': copper_feed_moles,
            'copper_extracted_moles': copper_extracted_moles,
            'K_ex': K_ex,
            'RH': RH,
            'O_A': O_A,
            'pH': pH
        }


class ElectrowinningStage:
    """
    Electrowinning (EW) stage: Plates pure copper metal from purified solution using electrolysis.

    Physical Parameters → Cost Calculation:
    - copper_plated → cell area/number → CAPEX (electrowinning cells)
    - copper_plated + cell_voltage → power consumption → OPEX (electricity cost)

    Key dependency: Uses purified copper solution from SX stage.
    """

    def __init__(self, current_density=250, current_efficiency=0.95, cell_cost_per_m2=2000):
        """
        Args:
            current_density: Current density (A/m2), typical: 200-300 A/m2
            current_efficiency: Current efficiency (0-1), typical: 0.90-0.98
            cell_cost_per_m2: Cost per unit cathode area ($/m2) - 1990s basis
                             Will be adjusted for inflation
                             Typical 1990s: $1500-2500/m2
                             Current (2020s): $3500-6000/m2
        """
        self.current_density = current_density
        self.current_efficiency = current_efficiency
        self.cell_cost_per_m2 = cell_cost_per_m2

    def calculate_capex(self, copper_plated_tpa):
        """
        CAPEX = electrowinning cell equipment cost

        Steps:
        1. Calculate required cathode area from production rate
        2. Estimate equipment cost per area
        3. Apply inflation and lang factors

        Args:
            copper_plated_tpa: Copper production rate (tonnes/year)

        Returns:
            Installed CAPEX ($)
        """
        # Faraday's law: kg Cu = (I * t * M_Cu * eta) / (n * F)
        # Rearrange for current: I = (kg Cu * n * F) / (t * M_Cu * eta)

        MW_Cu = 63.546  # g/mol
        F = 96485  # C/mol (Faraday constant)
        n = 2  # electrons per Cu atom

        # Annual copper production
        copper_kg_per_year = copper_plated_tpa * 1000

        # Operating time per year (seconds)
        operating_seconds = seconds_per_op_year

        # Required current (Amperes)
        required_current = (copper_kg_per_year * 1000 * n * F) / (operating_seconds * MW_Cu * self.current_efficiency)

        # Required cathode area (m2)
        # Current = current_density × area
        cathode_area = required_current / self.current_density

        # Equipment cost (cell cost per m2 is typically installed cost)
        equipment_cost = cathode_area * self.cell_cost_per_m2

        # Apply only inflation factor (cell_cost_per_m2 already includes installation)
        # Lang factor NOT applied since EW cell costs are quoted as installed costs
        installed_cost = equipment_cost * inflation_factor

        return installed_cost

    def calculate_opex(self, copper_plated_tpa, E_eq=0.34, eta_cath=0.1, eta_anode=0.3, IR=0.5):
        """
        OPEX = annual power cost for electrowinning

        Steps:
        1. Calculate required current from Faraday's law
        2. Calculate cell voltage
        3. Calculate power consumption
        4. Convert to annual cost

        Args:
            copper_plated_tpa: Copper production rate (tonnes/year)
            E_eq: Equilibrium potential (V), default 0.34 V for Cu
            eta_cath: Cathode overpotential (V), typical 0.1 V
            eta_anode: Anode overpotential (V), typical 0.3 V
            IR: Resistance losses (V), typical 0.5 V

        Returns:
            dict with power cost breakdown
        """
        MW_Cu = 63.546  # g/mol
        F = 96485  # C/mol
        n = 2  # electrons per Cu

        # Annual copper production
        copper_kg_per_year = copper_plated_tpa * 1000

        # Operating time per year (seconds)
        operating_seconds = seconds_per_op_year

        # Required current (Amperes)
        required_current = (copper_kg_per_year * 1000 * n * F) / (operating_seconds * MW_Cu * self.current_efficiency)

        # Cell voltage
        V_cell = cell_voltage(E_eq, eta_cath, eta_anode, IR)

        # Power (Watts)
        power_watts = required_current * V_cell
        power_kw = power_watts / 1000

        # Annual energy consumption
        annual_hours = hours_per_op_year
        annual_energy_kwh = power_kw * annual_hours
        annual_energy_mwh = annual_energy_kwh / 1000

        # Annual cost
        power_cost = annual_energy_mwh * electricity_cost

        return {
            'required_current': required_current,
            'cell_voltage': V_cell,
            'power_kw': power_kw,
            'annual_energy_kwh': annual_energy_kwh,
            'annual_energy_mwh': annual_energy_mwh,
            'power_cost': power_cost,
            'total': power_cost
        }

    def process(self, copper_extracted_kg, ew_recovery=0.95):
        """
        Calculate electrowinning production.

        Simple recovery model: assumes fixed recovery based on operating conditions.
        In practice, recovery depends on current efficiency and operating time.

        Args:
            copper_extracted_kg: Copper from SX (kg/year)
            ew_recovery: Electrowinning recovery (0-1), typical 0.93-0.98

        Returns:
            dict with keys:
                - 'ew_recovery': EW recovery (0-1)
                - 'copper_feed_kg': Copper feed from SX (kg/year)
                - 'copper_plated_kg': Copper plated (kg/year)
                - 'copper_plated_tpa': Copper plated (tonnes/year)
                - 'copper_remaining_kg': Copper remaining in solution (kg/year)
        """
        # Calculate copper plated
        copper_plated_kg = copper_extracted_kg * ew_recovery
        copper_remaining_kg = copper_extracted_kg * (1 - ew_recovery)
        copper_plated_tpa = copper_plated_kg / 1000

        return {
            'ew_recovery': ew_recovery,
            'copper_feed_kg': copper_extracted_kg,
            'copper_plated_kg': copper_plated_kg,
            'copper_plated_tpa': copper_plated_tpa,
            'copper_remaining_kg': copper_remaining_kg
        }


def test_grinding():
    """Test grinding stage with default parameters."""
    print("="*80)
    print("GRINDING STAGE - COST INTEGRATION TEST")
    print("="*80)

    # Create grinding stage with default parameters from tea.py
    grinding = GrindingStage()

    print(f"\nGrinding Parameters:")
    print(f"  Feed size: {grinding.feed_size} microns")
    print(f"  Product size: {grinding.product_size} microns")
    print(f"  Work index: {grinding.work_index:.2f} kWh/t")
    print(f"  Specific work: {grinding.work:.2f} kWh/t")

    # Test with example throughput
    throughput_tpa = 100000

    print(f"\nThroughput: {throughput_tpa} tpa")

    # Calculate costs
    capex = grinding.calculate_capex(throughput_tpa)
    opex = grinding.calculate_opex(throughput_tpa)

    print(f"\nCosts:")
    print(f"  CAPEX (installed): ${capex:,.0f}")
    print(f"  Annual OPEX: ${opex:,.0f}/year")

    # Show intermediate calculations for transparency
    required_power_kw = throughput_tpa * grinding.work / hours_per_op_year
    annual_power_kwh = throughput_tpa * grinding.work

    print(f"\nIntermediate Calculations:")
    print(f"  Required power: {required_power_kw:.2f} kW")
    print(f"  Annual energy: {annual_power_kwh:,.0f} kWh/year")
    print(f"  Annual energy: {annual_power_kwh/1000:,.2f} MWh/year")

    print("="*80)


def test_integrated_process():
    """Test integrated grinding to flotation process with costs."""
    print("\n" + "="*80)
    print("INTEGRATED PROCESS TEST: GRINDING -> FLOTATION -> COSTS")
    print("="*80)

    # Initialize stages
    grinding = GrindingStage()
    flotation = FlotationStage()

    throughput_tpa = 100000

    print(f"\n{'STEP 1: GRINDING':-^80}")
    print(f"Feed: {throughput_tpa:,} tpa")

    # Run grinding process
    tonnage_input = [10000, 30000, 40000, 15000, 5000]
    top_size_input = [600, 300, 150, 75, 37.5]
    bottom_size_input = [300, 150, 75, 37.5, 18.75]
    selection_input = [0.8, 0.6, 0.4, 0.2, 0.1]
    break_intensity = 5

    grinding_result = grinding.process(tonnage_input, top_size_input, bottom_size_input,
                                       selection_input, break_intensity)

    print(f"  Total feed: {grinding_result['total_feed']:,.0f} tonnes")
    print(f"  Total product: {grinding_result['total_product']:,.0f} tonnes")
    print(f"  P80 size: {grinding_result['p80_size']:.2f} microns")
    print(f"  Mass balance OK: {grinding_result['mass_balance_ok']}")

    # Grinding costs
    grinding_capex = grinding.calculate_capex(throughput_tpa)
    grinding_opex = grinding.calculate_opex(throughput_tpa)

    print(f"\nGrinding Costs:")
    print(f"  CAPEX: ${grinding_capex:,.0f}")
    print(f"  OPEX: ${grinding_opex:,.0f}/year")

    print(f"\n{'STEP 2: FLOTATION':-^80}")
    print(f"Feed from grinding: {grinding_result['total_product']:,.0f} tonnes")

    # Flotation costs (based on feed tonnage)
    flotation_capex = flotation.calculate_capex(throughput_tpa)
    flotation_opex = flotation.calculate_opex(throughput_tpa)

    print(f"\nFlotation Costs (for processing {throughput_tpa:,} tpa):")
    print(f"  CAPEX: ${flotation_capex:,.0f}")
    print(f"  OPEX (power): ${flotation_opex['power']:,.0f}/year")
    print(f"  OPEX (water): ${flotation_opex['water']:,.0f}/year")
    print(f"  OPEX (total): ${flotation_opex['total']:,.0f}/year")

    # Example: Run flotation recovery calculation
    print(f"\n{'FLOTATION RECOVERY CALCULATION':-^80}")
    print("(Using example parameters - would normally come from optimization)")

    # Example parameters for flotation
    fitting_params = {'b': 2, 'alpha': 0.10, 'coverage': 0.525,
                     'bubble_f': 0.825, 'detach_f': 0.5, 'bulk_zone': 0.5}
    flot_constants = {'specific_gravity': 4.2, 'grade': 0.28, 'permitivity': 8.854e-12,
                     'dielectric': 80, 'pe': 4, 'cell_area': 200, 'bbl_ratio': 10,
                     'feed_tonnage': throughput_tpa}
    opt_vars = {'particle_size': grinding_result['p80_size'], 'contact_angle': 52.3,
               'bubble_z_pot': -0.03, 'particle_z_pot': -50, 'sp_power': 1,
               'sp_gas_rate': 2, 'air_fraction': 0.05, 'slurry_fraction': 0.15,
               'num_cells': 1, 'ret_time': 23.67, 'cell_volume': 700,
               'froth_height': 0.165, 'frother_conc': 50}
    inputs_dict = {'frother_type': 4, 'water_or_particle': 'Particle'}

    flotation_result = flotation.process(fitting_params, flot_constants, opt_vars, inputs_dict)

    print(f"  Recovery: {flotation_result['recovery']*100:.2f}%")
    print(f"  Feed tonnage: {flotation_result['feed_tonnage']:,.0f} tpa")
    print(f"  Concentrate: {flotation_result['concentrate_tonnage']:,.0f} tpa")
    print(f"  Tailings: {flotation_result['tails_tonnage']:,.0f} tpa")

    print(f"\n{'DOWNSTREAM IMPLICATIONS':-^80}")
    print(f"Leaching/SX/EW will process: {flotation_result['concentrate_tonnage']:,.0f} tpa")
    print(f"  (NOT the full {throughput_tpa:,} tpa feed!)")
    print(f"\nThis is {flotation_result['concentrate_tonnage']/throughput_tpa*100:.1f}% of the feed tonnage.")

    print(f"\n{'STEP 3: LEACHING':-^80}")
    print(f"Feed from flotation: {flotation_result['concentrate_tonnage']:,.0f} tpa (concentrate)")

    # Initialize leaching stage with realistic sulfuric acid price
    # Industrial H2SO4 costs $50-150/tonne (vs $400 for concentrated HNO3)
    leaching = LeachingStage(acid_price_per_tonne=100)

    # Calculate copper in concentrate (assuming 28% Cu grade)
    concentrate_cu_grade = 0.28
    copper_in_concentrate_tpa = flotation_result['concentrate_tonnage'] * concentrate_cu_grade

    print(f"  Copper in concentrate: {copper_in_concentrate_tpa:,.0f} tpa Cu")

    # Run leaching process
    leaching_result = leaching.process(
        concentrate_tpa=flotation_result['concentrate_tonnage'],
        copper_in_concentrate_tpa=copper_in_concentrate_tpa
    )

    print(f"\n  Leaching recovery: {leaching_result['leach_recovery']*100:.2f}%")
    print(f"  Copper leached: {leaching_result['copper_leached_tpa']:,.0f} tpa Cu")
    print(f"  Copper leached: {leaching_result['copper_leached_kg']:,.0f} kg/year Cu")

    # Leaching costs (depends on concentrate tonnage and copper leached)
    leaching_capex = leaching.calculate_capex(flotation_result['concentrate_tonnage'])
    leaching_opex = leaching.calculate_opex(leaching_result['copper_leached_kg'])

    print(f"\nLeaching Costs:")
    print(f"  CAPEX: ${leaching_capex:,.0f}")
    print(f"  OPEX (acid): ${leaching_opex['acid_cost']:,.0f}/year")
    print(f"    - Stoichiometric acid: {leaching_opex['acid_stoichiometric_tpa']:,.0f} tpa")
    print(f"    - Acid recovery: {leaching_opex['acid_recovery_fraction']*100:.1f}%")
    print(f"    - Net acid consumption: {leaching_opex['acid_net_consumption_tpa']:,.0f} tpa")
    print(f"    - Acid price: ${leaching.acid_price}/tonne")

    print(f"\n{'DOWNSTREAM IMPLICATIONS':-^80}")
    print(f"SX/EW will process: {leaching_result['copper_leached_tpa']:,.0f} tpa Cu")
    print(f"  This is {leaching_result['leach_recovery']*100:.1f}% of the copper in concentrate")
    print(f"  Overall Cu recovery so far: {flotation_result['recovery']*leaching_result['leach_recovery']*100:.1f}%")

    print(f"\n{'STEP 4: SOLVENT EXTRACTION':-^80}")
    print(f"Feed from leaching: {leaching_result['copper_leached_tpa']:,.0f} tpa Cu")

    # Initialize SX stage
    sx = SolventExtractionStage(solvent_price_per_litre=5, solvent_losses_fraction=0.01)

    # Run SX process (using parameters from equation_flow_official.py)
    sx_result = sx.process(
        copper_leached_kg=leaching_result['copper_leached_kg'],
        K_ex=10,
        RH=0.2,
        O_A=1.2,
        pH=1.1
    )

    print(f"\n  SX Parameters:")
    print(f"    - K_ex: {sx_result['K_ex']}")
    print(f"    - RH (extractant conc): {sx_result['RH']} M")
    print(f"    - O/A ratio: {sx_result['O_A']}")
    print(f"    - pH: {sx_result['pH']}")

    print(f"\n  SX recovery: {sx_result['sx_recovery']*100:.2f}%")
    print(f"  Copper extracted: {sx_result['copper_extracted_kg']:,.0f} kg/year")
    print(f"  Copper in raffinate: {sx_result['copper_raffinate_kg']:,.0f} kg/year")

    # SX costs (depends on copper throughput and solution volume)
    sx_capex = sx.calculate_capex(leaching_result['copper_leached_tpa'])

    # Estimate solution volume (rough estimate: 1 m3 per kg Cu per year)
    solution_volume_m3 = leaching_result['copper_leached_kg'] * 0.001

    sx_opex = sx.calculate_opex(sx_result['copper_extracted_kg'], solution_volume_m3)

    print(f"\nSolvent Extraction Costs:")
    print(f"  CAPEX: ${sx_capex:,.0f}")
    print(f"  OPEX (solvent): ${sx_opex['solvent_cost']:,.0f}/year")
    print(f"    - Organic volume: {sx_opex['organic_volume_litres']:,.0f} L")
    print(f"    - Solvent losses: {sx_opex['solvent_losses_fraction']*100:.1f}%")
    print(f"    - Annual makeup: {sx_opex['solvent_losses_litres']:,.0f} L")
    print(f"    - Solvent price: ${sx.solvent_price}/L")

    print(f"\n{'DOWNSTREAM IMPLICATIONS':-^80}")
    print(f"EW will process: {sx_result['copper_extracted_kg']/1000:.0f} tpa Cu")
    print(f"  This is {sx_result['sx_recovery']*100:.1f}% of the copper from leaching")
    overall_recovery = flotation_result['recovery'] * leaching_result['leach_recovery'] * sx_result['sx_recovery']
    print(f"  Overall Cu recovery so far: {overall_recovery*100:.1f}%")

    print(f"\n{'STEP 5: ELECTROWINNING':-^80}")
    print(f"Feed from SX: {sx_result['copper_extracted_kg']/1000:.0f} tpa Cu")

    # Initialize EW stage (using 1990s cost basis, will be inflated)
    ew = ElectrowinningStage(current_density=250, current_efficiency=0.95, cell_cost_per_m2=2000)

    # Run EW process
    ew_result = ew.process(
        copper_extracted_kg=sx_result['copper_extracted_kg'],
        ew_recovery=0.95
    )

    print(f"\n  EW Parameters:")
    print(f"    - Current density: {ew.current_density} A/m2")
    print(f"    - Current efficiency: {ew.current_efficiency*100:.1f}%")

    print(f"\n  EW recovery: {ew_result['ew_recovery']*100:.2f}%")
    print(f"  Copper plated: {ew_result['copper_plated_kg']:,.0f} kg/year ({ew_result['copper_plated_tpa']:,.0f} tpa)")
    print(f"  Copper remaining: {ew_result['copper_remaining_kg']:,.0f} kg/year")

    # EW costs (depends on copper production)
    ew_capex = ew.calculate_capex(ew_result['copper_plated_tpa'])
    ew_opex = ew.calculate_opex(ew_result['copper_plated_tpa'])

    print(f"\nElectrowinning Costs:")
    print(f"  CAPEX: ${ew_capex:,.0f}")
    print(f"  OPEX (power): ${ew_opex['power_cost']:,.0f}/year")
    print(f"    - Required current: {ew_opex['required_current']:,.0f} A")
    print(f"    - Cell voltage: {ew_opex['cell_voltage']:.2f} V")
    print(f"    - Power: {ew_opex['power_kw']:,.0f} kW")
    print(f"    - Annual energy: {ew_opex['annual_energy_mwh']:,.0f} MWh")

    print(f"\n{'FINAL PRODUCT':-^80}")
    final_recovery = flotation_result['recovery'] * leaching_result['leach_recovery'] * sx_result['sx_recovery'] * ew_result['ew_recovery']
    print(f"Final copper product: {ew_result['copper_plated_tpa']:,.0f} tpa Cu (cathode)")
    print(f"Overall recovery: {final_recovery*100:.1f}%")
    print(f"  From {throughput_tpa:,} tpa ore feed")

    print(f"\n{'TOTAL COSTS (FULL PROCESS)':-^80}")
    total_capex = grinding_capex + flotation_capex + leaching_capex + sx_capex + ew_capex
    total_opex = grinding_opex + flotation_opex['total'] + leaching_opex['total'] + sx_opex['total'] + ew_opex['total']

    print(f"  Total CAPEX: ${total_capex:,.0f}")
    print(f"  Total OPEX: ${total_opex:,.0f}/year")
    print(f"\n  CAPEX Breakdown:")
    print(f"    - Grinding:       ${grinding_capex:,.0f} ({grinding_capex/total_capex*100:.1f}%)")
    print(f"    - Flotation:      ${flotation_capex:,.0f} ({flotation_capex/total_capex*100:.1f}%)")
    print(f"    - Leaching:       ${leaching_capex:,.0f} ({leaching_capex/total_capex*100:.1f}%)")
    print(f"    - SX:             ${sx_capex:,.0f} ({sx_capex/total_capex*100:.1f}%)")
    print(f"    - Electrowinning: ${ew_capex:,.0f} ({ew_capex/total_capex*100:.1f}%)")
    print(f"\n  OPEX Breakdown:")
    print(f"    - Grinding:       ${grinding_opex:,.0f}/year ({grinding_opex/total_opex*100:.1f}%)")
    print(f"    - Flotation:      ${flotation_opex['total']:,.0f}/year ({flotation_opex['total']/total_opex*100:.1f}%)")
    print(f"    - Leaching:       ${leaching_opex['total']:,.0f}/year ({leaching_opex['total']/total_opex*100:.1f}%)")
    print(f"    - SX:             ${sx_opex['total']:,.0f}/year ({sx_opex['total']/total_opex*100:.1f}%)")
    print(f"    - Electrowinning: ${ew_opex['total']:,.0f}/year ({ew_opex['total']/total_opex*100:.1f}%)")

    print("="*80)


if __name__ == "__main__":
    test_grinding()
    test_integrated_process()
