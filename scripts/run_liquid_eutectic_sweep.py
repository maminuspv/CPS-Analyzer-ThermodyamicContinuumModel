"""
scripts/run_liquid_eutectic_sweep.py
-----------------------------------
High-throughput screening script for liquid metal CPS coolants.
Identifies low-melting eutectics (Li-Sn-Bi-Pb-In-Ga) and evaluates interfacial 
chemical compatibility against Zr5Ta55Ni35Ti5 refractory armor.
"""

import itertools
import csv
import numpy as np
from typing import Dict, List, Tuple
from cps_analyzer.thermodynamics import get_miedema_pairwise_enthalpy
from cps_analyzer.affinity import derive_affinity_and_moduli

# Gas Constant (J / mol*K)
R_GAS = 8.314462618

# Elemental reference melting points (Kelvin) and atomic numbers
LIQUID_ELEMENTS_DB = {
    "Li": {"Tm_K": 453.65, "Z": 3,  "M_g_mol": 6.94},
    "Sn": {"Tm_K": 505.08, "Z": 50, "M_g_mol": 118.71},
    "Bi": {"Tm_K": 544.55, "Z": 83, "M_g_mol": 208.98},
    "Pb": {"Tm_K": 600.61, "Z": 82, "M_g_mol": 207.20},
    "In": {"Tm_K": 429.75, "Z": 49, "M_g_mol": 114.82},
    "Ga": {"Tm_K": 302.91, "Z": 31, "M_g_mol": 69.72},
}

# Optimal solid refractory armor baseline from Minus (2026)
REFRACTORY_ARMOR_COMP = {"Zr": 0.05, "Ta": 0.55, "Ni": 0.35, "Ti": 0.05}


def calculate_substrate_mean_affinity(substrate: Dict[str, float]) -> float:
    """Calculates the concentration-weighted mean ion affinity of the solid armor."""
    total_x = sum(substrate.values())
    a_bar = 0.0
    for el, frac in substrate.items():
        props = derive_affinity_and_moduli(el)
        a_bar += (frac / total_x) * props["A_i_eV"]
    return a_bar


def evaluate_liquid_coolant(comp: Dict[str, float], armor_a_bar: float) -> dict:
    """
    Computes thermodynamic liquidus depression, effective plasma charge Z_eff,
    and interfacial work-of-adhesion mismatch against the refractory substrate.
    """
    elements = list(comp.keys())
    x = np.array([comp[el] for el in elements], dtype=float)
    x /= np.sum(x)
    n = len(elements)

    # 1. Ideal and regular solution thermodynamic mixing
    ideal_tm = sum(x[i] * LIQUID_ELEMENTS_DB[elements[i]]["Tm_K"] for i in range(n))
    
    delta_h_mix = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            h_ij = get_miedema_pairwise_enthalpy(elements[i], elements[j])
            delta_h_mix += 4.0 * h_ij * x[i] * x[j]  # kJ/mol

    # Ideal entropy of mixing (kJ / mol*K)
    delta_s_ideal = -R_GAS * np.sum(x * np.log(x + 1e-12)) / 1000.0

    # Approximate eutectic depression via thermodynamic excess free energy:
    # Strong exothermic mixing (Delta_H_mix < 0) deepens liquid stability troughs
    delta_t_eutectic = (delta_h_mix / (delta_s_ideal + 0.008)) * 0.45
    t_liq_pred_k = ideal_tm + delta_t_eutectic
    t_liq_pred_c = t_liq_pred_k - 273.15

    # 2. Core plasma radiation metric: Effective atomic number Z_eff
    z_eff = sum(x[i] * LIQUID_ELEMENTS_DB[elements[i]]["Z"] for i in range(n))

    # 3. Interfacial affinity disparity (Capillary spreading driving force)
    # Delta_A_inter = |A_substrate - A_liquid|
    liquid_a_bar = 0.0
    for i in range(n):
        props = derive_affinity_and_moduli(elements[i])
        liquid_a_bar += x[i] * props["A_i_eV"]

    delta_a_interface = abs(armor_a_bar - liquid_a_bar)

    # 4. Armor dissolution risk: Check cross-mixing enthalpy into solid Ni and Zr
    # Highly exothermic negative cross-enthalpies indicate aggressive dissolution
    cross_dissolution_penalty = 0.0
    for i in range(n):
        cross_dissolution_penalty += x[i] * get_miedema_pairwise_enthalpy(elements[i], "Ni")

    return {
        "T_liq_K": round(t_liq_pred_k, 1),
        "T_liq_C": round(t_liq_pred_c, 1),
        "Delta_H_mix_kJ": round(delta_h_mix, 2),
        "Z_eff": round(z_eff, 2),
        "Delta_A_inter_eV": round(delta_a_interface, 2),
        "Armor_Corrosion_Index": round(cross_dissolution_penalty, 2),
    }


def generate_simplex_grid(elements: List[str], step_size: float = 0.05, min_pct: float = 0.05):
    """Generates compositional grid points constrained to a sum of 1.0."""
    n = len(elements)
    step_int = int(round(step_size * 100))
    min_int = int(round(min_pct * 100))
    valid_range = range(min_int, 100 + step_int, step_int)
    for combo in itertools.product(valid_range, repeat=n - 1):
        rem = 100 - sum(combo)
        if rem >= min_int and rem % step_int == 0:
            fracs = [c / 100.0 for c in combo] + [rem / 100.0]
            yield dict(zip(elements, fracs))


def run_eutectic_sweep(
    candidate_elements: List[str],
    step_size: float = 0.05,
    max_t_liq_c: float = 250.0,
    max_z_eff: float = 40.0
) -> List[dict]:
    """Sweeps candidate liquid systems and applies fusion CPS coolant guardrails."""
    armor_a_bar = calculate_substrate_mean_affinity(REFRACTORY_ARMOR_COMP)
    print(f"--> Target Substrate: Zr5Ta55Ni35Ti5 (Mean Affinity = {armor_a_bar:.2f} eV)")
    print(f"--> Sweeping liquid system {candidate_elements} at {step_size*100:.0f} mol% resolution...")

    grid = list(generate_simplex_grid(candidate_elements, step_size=step_size, min_pct=0.05))
    print(f"--> Permutations generated: {len(grid)}")

    passed = []
    for comp in grid:
        metrics = evaluate_liquid_coolant(comp, armor_a_bar)
        
        # Enforce liquidus ceiling (e.g. <= 250 deg C) and core plasma Z limit
        if metrics["T_liq_C"] <= max_t_liq_c and metrics["Z_eff"] <= max_z_eff:
            formula = "".join(f"{el}{int(round(comp[el]*100))}" for el in candidate_elements if comp[el] > 0)
            row = {"Alloy": formula}
            for el in candidate_elements:
                row[f"x_{el}"] = round(comp[el], 2)
            row.update(metrics)
            passed.append(row)

    # Sort primarily by lowest liquidus temperature
    passed.sort(key=lambda r: r["T_liq_C"])
    return passed


if __name__ == "__main__":
    # Screen ternary Li-Sn-In systems (Low melting, intermediate-Z plasma compatibility)
    target_elements = ["Li", "Sn", "In"]
    results = run_eutectic_sweep(
        target_elements,
        step_size=0.05,
        max_t_liq_c=250.0,
        max_z_eff=35.0
    )

    print("\n" + "=" * 92)
    print(f" TOP LOW-MELTING LIQUID CPS CANDIDATES (Found: {len(results)})")
    print("=" * 92)

    header = f"{'Alloy':<16} | {'T_liq (°C)':<10} | {'T_liq (K)':<10} | {'Z_eff':<7} | {'ΔH_mix':<8} | {'ΔA_inter':<8}"
    print(header)
    print("-" * len(header))

    for r in results[:15]:
        print(
            f"{r['Alloy']:<16} | {r['T_liq_C']:<10.1f} | {r['T_liq_K']:<10.1f} | "
            f"{r['Z_eff']:<7.1f} | {r['Delta_H_mix_kJ']:<8.2f} | {r['Delta_A_inter_eV']:<8.2f}"
        )

    # Save outputs to CSV
    output_filename = "data/liquid_eutectic_cps_candidates.csv"
    if results:
        fieldnames = list(results[0].keys())
        with open(output_filename, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        print(f"\nSaved {len(results)} viable coolants to '{output_filename}'.")