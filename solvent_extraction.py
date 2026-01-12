import math

upper_lower_bounds = {
    "K_ex": (2, 6),
    "RH": (0.05, 0.10),
    "H_plus": (0.01, 0.05),
    "O_A": (0.08, 0.18),
    "pH": (1.8, 2.5)
}

def solvent_extraction(C_Cu_in_aq, K_ex=10, RH=0.5, O_A=.5, pH=1.1):
    """
    Calculate output copper concentration from solvent extraction.

    Equation: C_Cu_out^aq = C_Cu_in^aq / (1 + K_ex * [RH]^2 / [H+]^2 * (O/A))

    Args:
        C_Cu_in_aq: Amount of Cu in feed from leaching
        K_ex: Extraction equilibrium constant
        RH: Extractant concentration (M)
        O_A: Organic/aqueous ratio
        pH: pH of aqueous solution

    Returns:
        float: Copper concentration in raffinate (C_Cu_out^aq)
    """

    H_plus = 10**(-pH)
    denominator = 1 + K_ex * (RH ** 2) / (H_plus ** 2) * O_A
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
