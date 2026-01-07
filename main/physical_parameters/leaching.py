import math

def copper_leaching_rate(Ea, T, P_O2, n, Phi, H_plus, MFeS2, X_Cu, k0=1.0, R=8.314):
    """
    Calculate the copper leaching rate based on shrinking core kinetics.

    Equation: dX_Cu/dt = 3*k0*exp(-Ea/RT)*(P_O2_eff)^n*(1-X_Cu)^(2/3)

    Args:
        Ea: Activation energy (kJ/mol)
        T: Temperature (°C)
        P_O2: Partial pressure of oxygen (bar)
        n: Stoichiometry parameter (dimensionless)
        Phi: Mass of solids per solution volume (kg/L)
        H_plus: H+ ion concentration, typically expressed as -log(pH)
        MFeS2: Mass fraction of FeS2 (dimensionless)
        X_Cu: Molar ratio of Cu converted (fraction, 0-1)
        k0: Pre-exponential factor (1/s), default=1.0
        R: Universal gas constant (J/mol·K), default=8.314

    Returns:
        float: Rate of copper conversion (dX_Cu/dt) in 1/s

    Base case parameters (from image):
        Ea: 70-90 kJ/mol
        T: 180-202°C
        P_O2: 10-15 bar
        n: 0.5-1
        Phi: 10 kg/L
        H_plus: -log(pH)
        MFeS2: varies
    """
    # Convert temperature from Celsius to Kelvin
    T_K = T + 273.15

    # Convert Ea from kJ/mol to J/mol for consistency with R
    Ea_J = Ea * 1000

    # Calculate effective oxygen pressure (P_O2_eff = P_O2 / (1 + alpha*MFeS2))
    # For now, assuming P_O2_eff = P_O2 based on the equation shown
    # This can be modified if the relationship P_O2/(1+alpha*MFeS2) needs to be applied
    P_O2_eff = P_O2

    # Calculate the leaching rate
    # dX_Cu/dt = 3*k0*exp(-Ea/RT)*(P_O2_eff)^n*(1-X_Cu)^(2/3)
    exponential_term = math.exp(-Ea_J / (R * T_K))
    pressure_term = P_O2_eff ** n
    conversion_term = (1 - X_Cu) ** (2/3)

    dX_Cu_dt = 3 * k0 * exponential_term * pressure_term * conversion_term

    return dX_Cu_dt


# Test leaching
if __name__ == "__main__":
    print("COPPER LEACHING TEST")
    print("=" * 80)

    # Base case parameters
    Ea = 80  # kJ/mol
    T = 190  # °C
    P_O2 = 12  # bar
    n = 0.75
    Phi = 10  # kg/L
    H_plus = 1.0  # -log(pH)
    MFeS2 = 0.1
    X_Cu = 0.0  # Starting conversion
    k0 = 1.0

    print("Leaching conditions:")
    print(f"  Temperature: {T}°C")
    print(f"  Oxygen pressure: {P_O2} bar")
    print(f"  Activation energy: {Ea} kJ/mol")
    print(f"  Stoichiometry parameter (n): {n}")
    print(f"  Initial conversion: {X_Cu*100}%")
    print()

    # Calculate initial leaching rate
    rate = copper_leaching_rate(Ea, T, P_O2, n, Phi, H_plus, MFeS2, X_Cu, k0)

    print(f"Initial leaching rate: {rate:.6e} s^-1")
    print(f"Initial leaching rate: {rate*3600:.6e} hr^-1")
    print()

    # Simulate leaching over time
    print("Leaching progression:")
    print(f"{'Time (hr)':>10} | {'Conversion':>12} | {'Rate (1/hr)':>15}")
    print("-" * 80)

    time_hours = [0, 1, 2, 4, 8, 16, 24, 48]
    for t in time_hours:
        # Simple approximation: X_Cu increases over time
        # In reality, this would require numerical integration
        # For now, just show how rate changes with conversion
        X_Cu_t = min(0.95, t / 50)  # Simplified conversion over time
        rate_t = copper_leaching_rate(Ea, T, P_O2, n, Phi, H_plus, MFeS2, X_Cu_t, k0)

        print(f"{t:10.1f} | {X_Cu_t*100:11.2f}% | {rate_t*3600:15.6e}")

    print()
    print("NOTE: This shows how leaching rate decreases as conversion increases")
    print("due to the (1-X_Cu)^(2/3) term (shrinking core model)")
