"""
cps_analyzer/kinetics.py
------------------------
Johnson-Samwer Cooperative Shear Model (CSM), Egami size polydispersity (delta),
and mismatch potential perturbations for glass transition (Tg) predictions.
"""

import numpy as np
from typing import Dict
from .continuum import calculate_vrh_elasticity
from .thermodynamics import calculate_mixing_enthalpy

# Constants
R_GAS = 8.314462618              # J/(mol*K)
C_CSM = 0.01135                  # STZ cooperative scaling constant
H_REF = 250.0                    # Reference normalization enthalpy (kJ/mol)
GFA_DELTA_THRESHOLD = 8.0        # Egami topological vitrification threshold (%)


def predict_glass_properties(
    composition: Dict[str, float],
    alpha: float = 0.000,
    beta: float = 0.384,
    use_derived_moduli: bool = False,
    enforce_gfa_gate: bool = False
) -> dict:
    """
    Evaluates topological mismatch, mismatch potential, and predicted Tg.
    """
    cont = calculate_vrh_elasticity(composition, use_derived_moduli=use_derived_moduli)
    delta_h_mix = calculate_mixing_enthalpy(composition)

    x = cont["x"]
    elements = cont["elements"]
    r_vals = cont["r_vals"]
    z_vals = cont["z_vals"]
    r_bar = cont["r_bar"]
    z_bar = cont["z_bar"]

    # Egami topological size polydispersity delta (%)
    delta_pct = np.sqrt(np.sum(x * ((1.0 - r_vals / r_bar) ** 2))) * 100.0

    # Structural and electronic mismatch perturbation Phi (Eq. 10)
    size_mis = 0.0
    charge_mis = 0.0
    n = len(elements)
    for i in range(n):
        for j in range(i + 1, n):
            size_mis += x[i] * x[j] * ((r_vals[i] - r_vals[j]) / r_bar) ** 2
            charge_mis += x[i] * x[j] * ((z_vals[i] - z_vals[j]) / z_bar) ** 2

    phi_mismatch = (alpha * size_mis) + (beta * charge_mis)

    # Johnson-Samwer CSM Glass Transition (Eq. 9)
    elastic_energy_scale = C_CSM * (cont["G_VRH"] * 1e9) * cont["V_m"] / R_GAS
    tg_base = elastic_energy_scale * (1.0 - (delta_h_mix / H_REF))
    tg_cal = elastic_energy_scale * (1.0 - (delta_h_mix / H_REF) + phi_mismatch)

    vitrifies = bool(delta_pct >= GFA_DELTA_THRESHOLD)
    tg_gated = tg_cal if (vitrifies or not enforce_gfa_gate) else None

    return {
        "Tg_Calibrated_K": round(tg_cal, 1) if tg_cal is not None else None,
        "Tg_Gated_K": round(tg_gated, 1) if tg_gated is not None else None,
        "Tg_Base_K": round(tg_base, 1),
        "G_VRH_GPa": round(cont["G_VRH"], 1),
        "K_VRH_GPa": round(cont["K_VRH"], 1),
        "Poisson_nu": round(cont["nu"], 3),
        "Delta_H_mix_kJ": round(delta_h_mix, 2),
        "delta_pct": round(delta_pct, 2),
        "Phi_mismatch": round(phi_mismatch, 4),
        "Vitrifies_GFA": vitrifies,
        "Mean_Affinity_eV": round(cont["A_bar"], 2)
    }