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
    def energy_use(variables=variables):
        recovery_rate = .85
        grade = .03
        tonnes_feed_per_year = 100000 # tonnes/year
        kg_feed_per_second = (tonnes_feed_per_year * 1000) / (365 * 24 * 3600) # kg/s
        acid_flow_rate = 6.17*kg_feed_per_second*grade*recovery_rate # kg H2SO4/s
        molar_mass_H2SO4 = 98.08 # g/mol
        molar_flow_rate = (acid_flow_rate * 1000) / molar_mass_H2SO4 # mol/s
        F = 96485 # C/mol
        z = 2 # mol electron transfer for sulphuric acid
        current_required = molar_flow_rate * F * z # Amperes
        cell_potential = 2 # Volts
        power_required = (current_required * cell_potential)/1000 # kW
        operating_capacity = 0.9
        operating_hours_per_year = 365 * 24 * operating_capacity # hours/year
        energy_required = power_required * operating_hours_per_year # kWh/year

        return energy_required

    def water_use():
        recovery_rate = .85
        grade = .03
        tonnes_feed_per_year = 100000
        kg_feed_per_second = (tonnes_feed_per_year * 1000) / (365 * 24 * 3600)
        acid_flow_rate = 6.17*kg_feed_per_second*grade*recovery_rate
        molar_mass_H2SO4 = 98.08
        molar_flow_rate = (acid_flow_rate * 1000) / molar_mass_H2SO4
        molar_mass_H2O = 18.015
        water_flow_rate = (molar_flow_rate * molar_mass_H2O) / 1000
        water_density = 1000
        water_volumetric_flow = water_flow_rate / water_density
        operating_capacity = 0.9
        operating_seconds_per_year = 365 * 24 * 3600 * operating_capacity
        annual_water_use = water_volumetric_flow * operating_seconds_per_year

        return annual_water_use
    
    def capex(area_required):
        capex_cost = 1500 * 50 * ((area_required / 1500) ** 0.85)
        return capex_cost

    def cost(energy_required, water_required):
        energy_cost_per_mwh = 50
        energy_cost = (energy_required / 1000) * energy_cost_per_mwh
        water_cost = water_required * 50
        cost_total = water_cost + energy_cost
        return energy_cost, water_cost, cost_total
    energy_required = energy_use()
    water_required = water_use()
    energy_cost, water_cost, cost_total = cost(energy_required, water_required)

    print(f"Energy Required: {energy_required:.2f} kWh/year")
    print(f"Water Required: {water_required:.2f} m³/year")
    print(f"Energy Cost: ${energy_cost:.2f}/year")
    print(f"Water Cost: ${water_cost:.2f}/year")
    print(f"Total Cost: ${cost_total:.2f}/year")

    return energy_required, water_required, cost_total



#how to calculate acid flow rate? we have total acid consumption and acid consumption by kg
# know:
    # production rate of copper
    # production rate of copper * 6.17
#what should we be using for sulphuric acid density?
# search for sulphuric acid density, acid concentration value
# 


#water flow rate = acid solution rate - (mass flow rate / 1000) 
# convert to m^3/s -> water flow rate / 1000