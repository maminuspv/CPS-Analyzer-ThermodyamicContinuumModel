"""
cps_analyzer/thermodynamics.py
------------------------------
Miedema multi-component mixing enthalpy calculations (analytical & curated table).
"""

import numpy as np
from typing import Dict
from .affinity import PERIODIC_METALS_DB

# Curated Takeuchi & Inoue (2005) Benchmark Lookup Table (kJ/mol)
MIEDEMA_CURATED_PAIRS = {
    ("Zr", "Cu"): -23.0, ("Zr", "Al"): -44.0, ("Zr", "Ni"): -49.0,
    ("Zr", "Ti"):   0.0, ("Zr", "Sn"): -43.0, ("Zr", "Ta"):   3.0,
    ("Cu", "Al"):  -1.0, ("Cu", "Ni"):   0.0, ("Cu", "Ti"):  -9.0,
    ("Cu", "Sn"):   7.0, ("Al", "Ni"): -22.0, ("Al", "Ti"): -30.0,
    ("Al", "Sn"):   4.0, ("Ni", "Ti"): -35.0, ("Ni", "Sn"): -21.0,
    ("Ti", "Sn"): -21.0, ("Fe", "Si"): -35.0, ("Fe", "B"):  -26.0,
    ("Si", "B"):  -14.0, ("Fe", "Al"): -11.0, ("Pd", "Ni"):   0.0,
    ("Pd", "P"):  -37.0, ("Ni", "P"):  -34.0, ("Al", "La"): -38.0,
    ("Ni", "La"): -27.0, ("Cu", "La"): -21.0, ("Ta", "Ni"): -29.0,
    ("Ta", "Ti"):   1.0, ("Ta", "Cu"):   2.0, ("Ta", "Al"): -19.0,
    ("W",  "Ni"):  -3.0, ("W",  "Zr"):  -9.0, ("W",  "Ti"):  -6.0,
    ("W",  "Ta"):  -7.0, ("Mo", "Ni"):  -7.0, ("Mo", "Zr"):  -6.0,
    ("Mo", "Ti"):  -4.0, ("Mo", "Ta"):  -5.0, ("Nb", "Ni"): -30.0,
    ("Nb", "Zr"):   4.0, ("Nb", "Ti"):   2.0, ("Nb", "Ta"):   0.0,
}

METALLOID_SET = {"B", "C", "Si", "Ge", "P", "As", "Sb", "Te"}
TRANSITION_SET = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Th", "U"
}


def get_miedema_pairwise_enthalpy(el1: str, el2: str) -> float:
    """Returns curated Miedema binary enthalpy, or computes via macroscopic-atom equations."""
    if el1 == el2:
        return 0.0
    if (el1, el2) in MIEDEMA_CURATED_PAIRS:
        return MIEDEMA_CURATED_PAIRS[(el1, el2)]
    if (el2, el1) in MIEDEMA_CURATED_PAIRS:
        return MIEDEMA_CURATED_PAIRS[(el2, el1)]

    # Analytical derivation for arbitrary periodic elements
    vm1, phi1, nws1 = PERIODIC_METALS_DB[el1][6], PERIODIC_METALS_DB[el1][8], PERIODIC_METALS_DB[el1][9]
    vm2, phi2, nws2 = PERIODIC_METALS_DB[el2][6], PERIODIC_METALS_DB[el2][8], PERIODIC_METALS_DB[el2][9]

    d_phi = phi1 - phi2
    d_nws = nws1 - nws2
    nws_inv_mean = 0.5 * ((1.0 / nws1) + (1.0 / nws2))
    v_surf_mean = 0.5 * ((vm1 ** (2.0 / 3.0)) + (vm2 ** (2.0 / 3.0)))

    both_tm = (el1 in TRANSITION_SET) and (el2 in TRANSITION_SET)
    p_const = 14.2 if both_tm else 12.35
    q_over_p = 9.4

    has_tm = (el1 in TRANSITION_SET) or (el2 in TRANSITION_SET)
    has_met = (el1 in METALLOID_SET) or (el2 in METALLOID_SET)
    r_over_p = 1.65 if (has_tm and has_met) else 0.0

    f_interface = (2.0 * v_surf_mean) / nws_inv_mean
    h_ij = f_interface * p_const * (-(d_phi ** 2) + q_over_p * (d_nws ** 2) - r_over_p)
    return float(np.clip(h_ij, -80.0, 50.0))


def calculate_mixing_enthalpy(composition: Dict[str, float]) -> float:
    """Computes total enthalpy of mixing Delta_H_mix for an N-component alloy (Eq. 8)."""
    elements = list(composition.keys())
    x = np.array([composition[el] for el in elements], dtype=float)
    x /= np.sum(x)

    delta_h_mix = 0.0
    n = len(elements)
    for i in range(n):
        for j in range(i + 1, n):
            h_ij = get_miedema_pairwise_enthalpy(elements[i], elements[j])
            delta_h_mix += 4.0 * h_ij * x[i] * x[j]

    return delta_h_mix