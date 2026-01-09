import math
import scipy.integrate

variables = {
    "Ea": (70, 90),  # kJ/mol
    "T": (180, 202),  # °C
    "P_O2": (10, 15),  # bar
    "n": (0.5, 1),  # dimensionless
    "Phi": 10,  # kg/L
    "H_plus": (.1-.3),  # -log(pH)
    "P(O2_eff)": (10-15)  # effective oxygen pressure
}

def leaching_copper_recovery(Ea, T, P_O2, n, H_plus, MFeS2, X_Cu, Phi=.1, k0=2500, R=8.314, leach_time=12):
    """
    Calculate total copper recovery after a given leaching time.

    Args:
        Ea: Activation energy (kJ/mol)
        T: Temperature (°C)
        P_O2: Partial pressure of oxygen (bar)
        n: Stoichiometry parameter (dimensionless)
        Phi: Mass of solids per solution volume (kg/L)
        H_plus: H+ ion concentration, typically expressed as -log(pH)
        MFeS2: Mass fraction of FeS2 (dimensionless)
        X_Cu: Initial copper conversion (not used, kept for compatibility)
        k0: Pre-exponential factor (1/s), default=1.0
        R: Universal gas constant (J/mol·K), default=8.314
        leach_time: Leaching time in hours

    Returns:
        X_Cu_final: Final copper recovery fraction (0-1)
        solution: ODE solution object
    """
    from scipy.integrate import solve_ivp

    MFeS2 = .86 * 1000 * .58
    T_K = T + 273.15
    Ea_J = Ea * 1000
    P_O2_eff = P_O2

    def copper_leaching_rate(t, y):
        """
        ODE function for solve_ivp.

        Equation: dX_Cu/dt = 3*k0*exp(-Ea/RT)*(P_O2_eff)^n*(1-X_Cu)^(2/3)

        Args:
            t: time (seconds)
            y: state vector [X_Cu]

        Returns:
            [dX_Cu/dt]
        """
        X_Cu_current = y[0]

        # Calculate the leaching rate
        exponential_term = math.exp(-Ea_J / (R * T_K))
        pressure_term = P_O2_eff ** n
        conversion_term = (1 - X_Cu_current) ** (2/3)

        dX_Cu_dt = 3 * k0 * exponential_term * pressure_term * conversion_term
        return [dX_Cu_dt]

    X_Cu_initial = [0.0]
    t_span = (0, leach_time * 3600)

    solution = solve_ivp(copper_leaching_rate, t_span, X_Cu_initial, dense_output=True)

    X_Cu_final = solution.y[0][-1]

    return X_Cu_final, solution


def calculate_acid_consumption(Cu_recovered_kg):
    """
    Calculate acid consumption for copper leaching based on chalcopyrite stoichiometry.

    Stoichiometry for chalcopyrite (CuFeS2):
    CuFeS2 + 4H2SO4 + O2 -> CuSO4 + FeSO4 + 4S + 4H2O

    Args:
        Cu_recovered_kg: Mass of copper recovered (kg)

    Returns:
        dict with acid consumption data
    """
    MW_Cu = 63.546  # g/mol
    MW_H2SO4 = 98.079  # g/mol

    # Stoichiometric ratio: 4 mol H2SO4 per 1 mol Cu
    acid_stoich_ratio = 4 * MW_H2SO4 / MW_Cu  # kg H2SO4 / kg Cu = 6.17

    # Total acid consumption
    acid_consumption_kg = Cu_recovered_kg * acid_stoich_ratio

    return {
        "acid_consumption_kg": acid_consumption_kg,
        "acid_per_kg_cu": acid_stoich_ratio,
        "Cu_recovered_kg": Cu_recovered_kg
    }

if __name__ == "__main__":
    print("Testing leaching_copper_recovery:")
    X_Cu_final, solution = leaching_copper_recovery(Ea=80, T=190, P_O2=12, n=0.75, H_plus=0.2, MFeS2=0.5, X_Cu=0)
    print(f"Final copper recovery: {X_Cu_final:.4f} ({X_Cu_final*100:.2f}%)")

    print("\nTesting calculate_acid_consumption:")
    result = calculate_acid_consumption(Cu_recovered_kg=1000)
    print(f"Copper recovered: {result['Cu_recovered_kg']} kg")
    print(f"Acid consumption: {result['acid_consumption_kg']:.2f} kg")
    print(f"Acid per kg Cu: {result['acid_per_kg_cu']:.2f}")


"""
total acid consumption for the system should be 6.17 * feed * 

"""