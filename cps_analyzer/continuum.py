"""
cps_analyzer/continuum.py
-------------------------
Voigt-Reuss-Hill (VRH) continuum elastic homogenization and Poisson's ratio calculations.
"""

import numpy as np
from typing import Dict, List
from .affinity import derive_affinity_and_moduli


def calculate_vrh_elasticity(composition: Dict[str, float], use_derived_moduli: bool = False) -> dict:
    """
    Computes effective bulk modulus (K_VRH), shear modulus (G_VRH),
    Poisson's ratio (nu), and composite molar volume (V_m).
    
    Args:
        composition: Dict of mole fractions, e.g., {'Zr': 0.05, 'Ta': 0.55, 'Ni': 0.35, 'Ti': 0.05}
        use_derived_moduli: If True, uses Rose UBER moduli; if False, uses tabulated reference moduli.
    """
    elements = list(composition.keys())
    x = np.array([composition[el] for el in elements], dtype=float)
    x /= np.sum(x)

    k_vals, g_vals, v_vals, r_vals, z_vals, a_vals = [], [], [], [], [], []
    for el in elements:
        props = derive_affinity_and_moduli(el)
        r_vals.append(props["r_pm"])
        z_vals.append(props["Z_val"])
        a_vals.append(props["A_i_eV"])
        v_vals.append(props["V_m_cm3"] * 1e-6)  # Convert to m^3/mol

        if use_derived_moduli:
            k_vals.append(props["K_der_GPa"])
            g_vals.append(props["G_der_GPa"])
        else:
            k_vals.append(props["K_tab_GPa"])
            g_vals.append(props["G_tab_GPa"])

    k_vals = np.array(k_vals)
    g_vals = np.array(g_vals)
    v_vals = np.array(v_vals)
    r_vals = np.array(r_vals)
    z_vals = np.array(z_vals)
    a_vals = np.array(a_vals)

    # Voigt and Reuss averaging
    k_voigt = np.sum(x * k_vals)
    g_voigt = np.sum(x * g_vals)

    k_reuss = 1.0 / np.sum(x / k_vals)
    g_reuss = 1.0 / np.sum(x / g_vals)

    # Voigt-Reuss-Hill effective moduli
    k_vrh = 0.5 * (k_voigt + k_reuss)
    g_vrh = 0.5 * (g_voigt + g_reuss)

    # Poisson's ratio (Eq. 6)
    nu = (3.0 * k_vrh - 2.0 * g_vrh) / (2.0 * (3.0 * k_vrh + g_vrh))

    # Concentration-weighted averages
    r_bar = np.sum(x * r_vals)
    z_bar = np.sum(x * z_vals)
    v_m = np.sum(x * v_vals)
    a_bar = np.sum(x * a_vals)

    return {
        "K_VRH": k_vrh,
        "G_VRH": g_vrh,
        "nu": nu,
        "V_m": v_m,
        "r_bar": r_bar,
        "z_bar": z_bar,
        "A_bar": a_bar,
        "elements": elements,
        "x": x,
        "r_vals": r_vals,
        "z_vals": z_vals
    }