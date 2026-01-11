#need to know acid or base flow rate

#electrolyzer needs to know total amount of acid required from leaching and solvent exchange

variables = {
"sl_ratio": .05,
"acid_concentration": .5,
"kg_feed": 1000,
"acid_flow_rate": 6.17*1000, # total from leaching and solvent exchange
"molar_flow_rate": .5 * 6.17*1000,
"F": 96485,  # C/mol
"z": 2, # mol electron transfer for suplhuric acid
"current_required": .5 * 6.17*1000 * 96485 * 2,
"cell_potential": 2,
"sulph_acid_density": 1830, #kg/m^3
"recovery_rate": .85,
"grade": .03
}

#electrolyzer - energy cost = energy_required * $50mw/h

def electrolyzer(variables=variables):
    def energy_use(variables=variables, F = 96485, z=2):
        recovery_rate = .85
        grade = .03
        tons_feed_per_year = 100000 # tons/year
        kg_feed_per_second = (tons_feed_per_year * 1000) / (365 * 24 * 3600) # kg/s
        acid_flow_rate = 6.17*kg_feed_per_second*grade*recovery_rate # kg H2SO4/s
        molar_mass_H2SO4 = 98.08 # g/mol
        molar_flow_rate = (acid_flow_rate * 1000) / molar_mass_H2SO4 # mol/s
        F = 96485 # C/mol
        z = 2 # mol electron transfer for sulphuric acid
        area_required = molar_flow_rate * z * F
        current_required = molar_flow_rate * F * z # Amperes
        cell_potential = 2 # Volts
        power_required = (current_required * cell_potential)/1000 # kW
        operating_capacity = 0.9
        operating_hours_per_year = 365 * 24 * operating_capacity # hours/year
        energy_required = power_required * operating_hours_per_year # kWh/year

        return energy_required, area_required

    def water_use(tons_feed_year=100000):
        """
        Calculate water usage for electrolyzer based on stoichiometry.

        For H2SO4 electrolysis: H2SO4 -> H2 + SO2 + O2 + H2O
        Water is produced during electrolysis.

        However, water is also needed for the acid solution preparation.
        The relationship from the Excel model shows:
        - Acid solution flow rate (m3/s) is based on concentration and molar flow
        - Water flow rate = (Acid solution mass flow - Acid mass flow)

        Returns:
            annual_water_use (float): Water consumption in m³/year
        """
        recovery_rate = .85
        grade = .03
        tons_feed_year = tons_feed_year
        kg_feed_per_second = (tons_feed_year * 1000) / (365 * 24 * 3600)

        # Acid flow rate calculation
        acid_flow_rate = 6.17*kg_feed_per_second*grade*recovery_rate  # kg H2SO4/s

        # Convert to molar flow rate
        molar_mass_H2SO4 = 98.08  # g/mol
        molar_flow_rate = (acid_flow_rate * 1000) / molar_mass_H2SO4  # mol/s

        # Calculate acid solution flow rate based on concentration
        acid_concentration = 0.5  # M (mol/L) = 500 mol/m³
        acid_solution_volumetric_flow = molar_flow_rate / (acid_concentration * 1000)  # m³/s

        # Calculate mass flow rate of acid solution
        # Assuming density of 0.5 M H2SO4 solution ≈ 1030 kg/m³
        acid_solution_density = 1030  # kg/m³
        acid_solution_mass_flow = acid_solution_volumetric_flow * acid_solution_density  # kg/s

        # Water flow rate = total solution - acid
        water_flow_rate = acid_solution_mass_flow - acid_flow_rate  # kg/s

        # Convert to volumetric flow (water density = 1000 kg/m³)
        water_density = 1000  # kg/m³
        water_volumetric_flow = water_flow_rate / water_density  # m³/s

        # Calculate annual water use based on operating capacity
        operating_capacity = 0.9
        operating_seconds_per_year = 365 * 24 * 3600 * operating_capacity  # s/year
        annual_water_use = water_volumetric_flow * operating_seconds_per_year  # m³/year

        return annual_water_use, water_flow_rate


    def capex(area_required):
        capex_cost = 1500 * 50 * ((area_required / 1500) ** 0.85)
        return capex_cost

    def cost(energy_required, water_required, area_required, water_cost_per_m3=0.3):
        """
        Calculate costs for electrolyzer including CAPEX and OPEX.

        Args:
            energy_required (float): Energy consumption in kWh/year
            water_required (float): Water consumption in m³/year
            area_required (float): Electrolyzer area in m²
            water_cost_per_m3 (float): Water cost in $/m³ (default 0.3 from Excel model)

        Returns:
            tuple: (capex_cost, energy_cost, water_cost, opex_total, cost_total) all in $
        """
        capex_cost = capex(area_required)
        energy_cost_per_mwh = 50  # $/MWh
        energy_cost = (energy_required / 1000) * energy_cost_per_mwh  # $/year
        water_cost = water_required * water_cost_per_m3  # $/year
        opex_total = water_cost + energy_cost  # $/year
        cost_total = capex_cost + opex_total  # $
        return capex_cost, energy_cost, water_cost, opex_total, cost_total
    energy_required, area_required = energy_use()
    water_required, water_flow_rate = water_use()
    capex_cost, energy_cost, water_cost, opex_total, cost_total = cost(energy_required, water_required, area_required)
    

    # print(f"Energy Required: {energy_required:.2f} kWh/year")
    # print(f"Water Required: {water_flow_rate:.2f} m³/year")
    # print(f"Area Required: {area_required:.2f} m²")
    # print(f"CAPEX: ${capex_cost:,.2f}")
    # print(f"Energy Cost: ${energy_cost:.2f}/year")
    # print(f"Water Cost: ${water_cost:.2f}/year")
    # print(f"Total OPEX: ${opex_total:.2f}/year")
    # print(f"Total Cost (CAPEX + OPEX): ${cost_total:,.2f}")

    return energy_required, water_required, capex_cost, opex_total, cost_total, water_flow_rate



# TEST ELECTROLYZER WATER USAGE AND COST
def test_electrolyzer_water_cost():
    """
    Test function to evaluate electrolyzer water usage and costs.
    Based on Excel model relationships.
    """
    # print("="*80)
    # print("ELECTROLYZER WATER USAGE AND COST ANALYSIS")
    # print("="*80)
    # print()
    # Input parameters
    tons_feed_per_year = 100000
    recovery_rate = 0.85
    grade = 0.03
    operating_capacity = 0.9
    water_cost_per_m3 = 0.3
    
    energy_required, water_required, capex_cost, opex_total, cost_total, water_flow_rate= electrolyzer(variables)
    # print(f"\nINPUT PARAMETERS:")
    
    # print(f"  Feed rate: {tons_feed_per_year:,.0f} tons/year")
    print(f"water flow rate {water_flow_rate}")
    # print(f"  Grade: {grade*100:.1f}%")
    # print(f"  Recovery rate: {recovery_rate*100:.1f}%")
    # print(f"  Operating capacity: {operating_capacity*100:.1f}%")
    # print(f"  Water cost: ${water_cost_per_m3}/m³")

    # Calculate flow rates
    kg_feed_per_second = (tons_feed_per_year * 1000) / (365 * 24 * 3600)
    acid_flow_rate = 6.17 * kg_feed_per_second * grade * recovery_rate

    molar_mass_H2SO4 = 98.08
    molar_flow_rate = (acid_flow_rate * 1000) / molar_mass_H2SO4

    acid_concentration = 0.5
    acid_solution_volumetric_flow = molar_flow_rate / (acid_concentration * 1000)

    # Water flow using Excel relationship
    water_volumetric_flow = acid_solution_volumetric_flow * 1.2185
    operating_seconds_per_year = 365 * 24 * 3600 * operating_capacity
    annual_water_use = water_volumetric_flow * operating_seconds_per_year


    # print(f"\nFLOW CALCULATIONS:")
    # print(f"  Feed rate: {kg_feed_per_second:.6f} kg/s")
    # print(f"  Acid flow rate: {acid_flow_rate:.6f} kg H2SO4/s")
    # print(f"  Molar flow rate: {molar_flow_rate:.6f} mol/s")
    # print(f"  Acid solution flow: {acid_solution_volumetric_flow:.6f} m³/s")
    # print(f"  Water flow: {water_volumetric_flow:.6f} m³/s")

    # print(f"\nANNUAL WATER USAGE:")
    # print(f"  Operating seconds/year: {operating_seconds_per_year:,.0f} s")
    # print(f"  Total water usage: {annual_water_use:,.2f} m³/year")

    # Calculate water cost
    annual_water_cost = annual_water_use * water_cost_per_m3

    # print(f"\nANNUAL WATER COST:")
    # print(f"  Water cost: ${annual_water_cost:,.2f}/year")
    # print(f"  Water cost per ton feed: ${annual_water_cost/tons_feed_per_year:.2f}/ton")

    # Energy calculation for comparison
    F = 96485
    z = 2
    current_required = molar_flow_rate * F * z
    cell_potential = 2
    power_required = (current_required * cell_potential) / 1000
    operating_hours_per_year = 365 * 24 * operating_capacity
    energy_required = power_required * operating_hours_per_year
    energy_cost_per_mwh = 50
    energy_cost = (energy_required / 1000) * energy_cost_per_mwh

    # print(f"\nENERGY COST:")
    # print(f"  Energy required: {energy_required:,.2f} kWh/year")
    # print(f"  Energy cost: ${energy_cost:,.2f}/year")

    # CAPEX calculation
    current_density = 4000  # A/m² from Excel
    area_required = current_required / current_density  # m²
    capex_per_m2 = 1500  # $/m² from Excel
    basis_area = 50  # m²
    exponent = 0.85
    capex_cost = capex_per_m2 * basis_area * ((area_required / capex_per_m2) ** exponent)

    # print(f"\nCAPEX COST:")
    # print(f"  Current required: {current_required:,.2f} A")
    # print(f"  Current density: {current_density:,.0f} A/m²")
    # print(f"  Area required: {area_required:,.2f} m²")
    # print(f"  CAPEX: ${capex_cost:,.2f}")

    # print(f"\nTOTAL OPERATING COST (OPEX):")
    total_opex = annual_water_cost + energy_cost
    # print(f"  Water + Energy: ${total_opex:,.2f}/year")
    

    # print(f"\nTOTAL COST (CAPEX + Annual OPEX):")
    # print(f"  CAPEX + OPEX: ${capex_cost + total_opex:,.2f}")
    # print(f"  OPEX per ton feed: ${total_opex/tons_feed_per_year:.2f}/ton")
    # print(f"  CAPEX per ton feed: ${capex_cost/tons_feed_per_year:.2f}/ton")
    # print(f"  OPEX per ton feed: ${(capex_cost+total_opex)/tons_feed_per_year:.2f}/ton")


    # print("="*80)

    return annual_water_use, annual_water_cost, capex_cost, total_opex

if __name__ == "__main__":
    test_electrolyzer_water_cost()