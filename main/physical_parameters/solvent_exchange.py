import math

upper_lower_bounds = {
    "K_ex": (100,300),  # Extraction equilibrium constant (dimensionless)
    "RH": (0.1,.35),  # Extractant concentration (molar)
    "H_plus": (0.01,0.05),  # Hydrogen ion concentration (molar)
    "O_A": (.05,1.5)  # Organic/aqueous solution ratio (dimensionless)
}

def solvent_extraction(C_Cu_in_aq, K_ex=10, RH=0.1, O_A=1, pH=1.1):
    """
    Calculate output copper concentration from solvent extraction.

    Equation: C_Cu_out^aq = C_Cu_in^aq / (1 + K_ex * [RH]^2 / [H+]^2 * (O/A))

    Args:
        C_Cu_in_aq: Amount of Cu in feed from leaching
        K_ex: Extraction equilibrium constant (dimensionless)
        RH: Extractant concentration (molar)
        H_plus: Hydrogen ion concentration (molar)
        O_A: Organic/aqueous solution ratio (dimensionless)

    Returns:
        float: Output copper concentration (C_Cu_out^aq) in moles
    """

    H_plus = 10**(-pH)  # convert pH to H+ concentration: [H+] = 10^(-pH)

    # Calculate the denominator term
    denominator = 1 + K_ex * (RH ** 2) / (H_plus ** 2) * O_A

    # Calculate output concentration
    C_Cu_out_aq = C_Cu_in_aq / denominator

    return C_Cu_out_aq

# Test example
if __name__ == "__main__":
    test_input = 1000
    result = solvent_extraction(test_input)
    extraction_rate = ((test_input - result) / test_input) * 100
    print(f'Input: {test_input} moles Cu')
    print(f'Output (raffinate): {result:.2f} moles Cu')
    print(f'Extracted: {test_input - result:.2f} moles Cu')
    print(f'Extraction efficiency: {extraction_rate:.2f}%')
