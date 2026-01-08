import math
import numpy as np

def electrowinning_concentration(C_Cu_in, I, eta_CE=.98, n=2, F=96485):
    """
    Calculate output copper concentration from electrowinning.

    Equation: C_Cu_out = C_Cu_in - (eta_CE * I) / (n * F * F)

    Args:
        C_Cu_in: Input copper concentration (g/mol or appropriate units)
        eta_CE: Current efficiency (dimensionless, 0-1)
        I: Current (Amperes)
        n: Number of electrons per Cu (2 e- per Cu)
        F: Faraday constant (96485 C/mol)
        M_Cu: Molecular weight of Cu (63.55 g/mol)

    Returns:
        float: Output copper concentration (C_Cu_out)
    """
    # Calculate concentration change due to electrowinning

    C_Cu_out = C_Cu_in - (eta_CE * I) / (n * F)

    return C_Cu_out

def cell_voltage(E_eq, eta_cath = .1, eta_anode = .3, IR = .5):
    """
    Calculate cell voltage for electrowinning.

    Equation: V_cell = E_eq + eta_cath + eta_anode + IR

    Args:
        E_eq: Nernst equilibrium potential (V)
        eta_cath: Overpotentials at cathode (V)
        eta_anode: Overpotentials at anode (V)
        I: Current (A)
        R: Resistance (Ohms)

    Returns:
        float: Cell voltage (V)
    """
    V_cell = E_eq + eta_cath + eta_anode + IR

    return V_cell

def nernst_potential(E_Cu_0=.34, a_Cu2=1, R=8.314, T=298.15, F=96485, n=2):
    """
    Calculate Nernst equilibrium potential.

    Equation: E_eq = E_Cu^0 + (RT/2F) * ln(a_Cu2+)

    Args:
        E_Cu_0: Standard potential for Cu (V), typically 0.34 V
        a_Cu2: Activity of Cu2+ ions (dimensionless)
        R: Universal gas constant (J/mol·K), default=8.314
        T: Temperature (K), default=298.15
        F: Faraday constant (C/mol), default=96485
        n: Number of electrons (2 for Cu), default=2

    Returns:
        float: Nernst equilibrium potential (V)
    """
    E_eq = E_Cu_0 + (R * T / (n * F)) * math.log(a_Cu2)

    return E_eq

def overpotential(a, j, b):
    """
    Calculate overpotential as function of current density.

    Equation: eta = a * ln(j) + b

    Args:
        a: Tafel slope coefficient (V)
        j: Current density (I/A)
        b: Constant (V)

    Returns:
        float: Overpotential (V)
    """
    eta = a * math.log(j) + b

    return eta

if __name__ == "__main__":
    # Example usage
    C_Cu_in = 10.0  # g/L
    eta_CE = 0.95
    I = 500.0  # Amperes

    C_Cu_out = electrowinning_concentration(C_Cu_in, eta_CE, I)
    print(f"Output Copper Concentration: {C_Cu_out:.2f} g/L")

    E_eq = nernst_potential(a_Cu2=0.1)
    V_cell = cell_voltage(E_eq)
    print(f"Cell Voltage: {V_cell:.2f} V")