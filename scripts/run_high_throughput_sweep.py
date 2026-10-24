"""
scripts/run_high_throughput_sweep.py
------------------------------------
Executes high-throughput combinatorial screening across the quaternary 
Zr-Ta-Ni-Ti Capillary Porous System (CPS) alloy space using first-principles 
derived ion affinities (A_i), Rose UBER elastic moduli, and Egami vitrification gates.
"""

import sys
from pathlib import Path
import csv
import itertools
import numpy as np
import pandas as pd

# Automatically append project root to sys.path so cps_analyzer is discoverable
project_root = Path(__file__).resolve().parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from cps_analyzer.affinity import PERIODIC_METALS_DB, derive_affinity_and_moduli
from cps_analyzer.kinetics import predict_glass_properties


def generate_simplex_grid(elements: list, step_size: float = 0.05, min_pct: float = 0.05):
    """Generates exact stoichiometric combinations summing to 1.0 across N elements."""
    n = len(elements)
    step_int = int(round(step_size * 100))
    min_int = int(round(min_pct * 100))
    valid_range = range(min_int, 100 + step_int, step_int)
    
    for combo in itertools.product(valid_range, repeat=n - 1):
        rem = 100 - sum(combo)
        if rem >= min_int and rem % step_int == 0:
            fracs = [c / 100.0 for c in combo] + [rem / 100.0]
            yield dict(zip(elements, fracs))


def run_high_throughput_sweep(
    target_elements: list = ["Zr", "Ta", "Ni", "Ti"],
    step_size: float = 0.05,
    min_tg: float = 700.0,
    min_g: float = 45.0,
    min_delta: float = 8.0
):
    print("=" * 78)
    print("       HIGH-THROUGHPUT QUATERNARY CPS ALLOY SCREENING ENGINE         ")
    print("=" * 78)
    print(f"Target System : {'-'.join(target_elements)}")
    print(f"Grid Step     : {step_size * 100:.0f} mol%")
    print(f"Constraints   : Tg >= {min_tg} K | G >= {min_g} GPa | delta >= {min_delta}% (Egami Gate)\n")

    grid_permutations = list(generate_simplex_grid(target_elements, step_size=step_size, min_pct=0.05))
    print(f"Total permutations generated in compositional simplex: {len(grid_permutations)}")

    passed_candidates = []

    for comp in grid_permutations:
        # Build formula string for readout (e.g., Zr5Ta55Ni35Ti5)
        formula = "".join([f"{el}{int(round(frac*100))}" for el, frac in comp.items()])

        try:
            # Predict properties using the unified first-principles analytical engine
            res = predict_glass_properties(
                comp, 
                alpha=0.000, 
                beta=0.384, 
                use_derived_moduli=True, 
                enforce_gfa_gate=True
            )

            tg_val = res["Tg_Gated_K"]
            g_val = res["G_VRH_GPa"]
            delta_val = res["delta_pct"]

            # Apply fusion CPS engineering guardrails
            if tg_val is not None and tg_val >= min_tg and g_val >= min_g and delta_val >= min_delta:
                row = {
                    "Composition": formula,
                    **{f"x_{el}": round(frac, 3) for el, frac in comp.items()},
                    "Tg_K": round(tg_val, 1),
                    "G_VRH_GPa": round(g_val, 1),
                    "K_VRH_GPa": round(res["K_VRH_GPa"], 1),
                    "Poisson_nu": round(res["Poisson_nu"], 3),
                    "Delta_H_mix_kJ": round(res["Delta_H_mix_kJ"], 2),
                    "delta_pct": round(delta_val, 2),
                }
                passed_candidates.append(row)

        except Exception as e:
            # Skip invalid configurations gracefully
            continue

    df_results = pd.DataFrame(passed_candidates)

    if not df_results.empty:
        df_results = df_results.sort_values(by="Tg_K", ascending=False).reset_index(drop=True)
        print(f"\nScreening Complete! Found {len(df_results)} qualifying candidates meeting all constraints.")
        print("\nTop 5 Optimal Refractory CPS Candidates:")
        print(df_results[["Composition", "Tg_K", "G_VRH_GPa", "K_VRH_GPa", "Poisson_nu", "Delta_H_mix_kJ", "delta_pct"]].head(5).to_string(index=False))

        # Save results to CSV
        output_dir = project_root / "data"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "quaternary_cps_screening_results.csv"
        df_results.to_csv(output_path, index=False)
        print(f"\nFull screening matrix successfully saved to: {output_path}")
    else:
        print("\nNo candidates met the strict engineering thresholds. Try relaxing min_tg or min_g bounds.")

    return df_results


if __name__ == "__main__":
    run_high_throughput_sweep()