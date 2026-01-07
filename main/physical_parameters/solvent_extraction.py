def solvent_extraction(C_Cu_in_aq, K_ex, RH, H_plus, O_A):
    """
    Calculate output copper concentration from solvent extraction.

    Equation: C_Cu_out^aq = C_Cu_in^aq / (1 + K_ex * [RH]^2 / [H+]^2 * (O/A))

    Args:
        C_Cu_in_aq: Amount of Cu in feed solution to SX extractant (moles of Cu)
        K_ex: Extraction equilibrium constant (dimensionless)
        RH: Extractant concentration (molar)
        H_plus: Hydrogen ion concentration (molar)
        O_A: Organic/aqueous solution ratio (dimensionless)

    Returns:
        float: Output copper concentration (C_Cu_out^aq) in moles
    """
    # Calculate the denominator term
    denominator = 1 + K_ex * (RH ** 2) / (H_plus ** 2) * O_A

    # Calculate output concentration
    C_Cu_out_aq = C_Cu_in_aq / denominator

    return C_Cu_out_aq