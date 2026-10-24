"""
cli/terminal_frontend.py
-------------------------
Interactive Terminal Command Line Interface (CLI) for the 76-Element
Universal CPS Alloy Screener & Glass Transition Predictor.
"""

import sys
from pathlib import Path

# Automatically append project root to sys.path if run directly
project_root = Path(__file__).resolve().parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from cps_analyzer.affinity import PERIODIC_METALS_DB
from cps_analyzer.kinetics import predict_glass_properties


def print_banner():
    print("=" * 78)
    print("       76-ELEMENT UNIVERSAL CPS ALLOY & GLASS SCREENING TERMINAL CLI      ")
    print("=" * 78)
    print(f"Loaded Elements in Database: {len(PERIODIC_METALS_DB)}")
    print("Instructions:")
    print("  - Enter element symbols and atomic percentages (e.g., Zr 5, Ta 55, Ni 35, Ti 5)")
    print("  - Type 'preset' to load standard benchmarks (Vitreloy, Refractory Armor, etc.)")
    print("  - Type 'list' to view all available 76 elements in the database")
    print("  - Type 'exit' to quit.\n")


def run_cli():
    presets = {
        "1": ("Zr55Cu30Al10Ni5 (Vitreloy 1 BMG)", {"Zr": 55, "Cu": 30, "Al": 10, "Ni": 5}),
        "2": ("Zr5Ta55Ni35Ti5 (Fusion CPS Armor)", {"Zr": 5, "Ta": 55, "Ni": 35, "Ti": 5}),
        "3": ("Fe78Si9B13 (Amorphous Ribbon)", {"Fe": 78, "Si": 9, "B": 13}),
        "4": ("Pd40Ni40P20 (High Stability BMG)", {"Pd": 40, "Ni": 40, "P": 20}),
    }

    print_banner()

    while True:
        current_comp = {}
        print("-" * 78)
        user_input = input("Enter command [element pct / preset / list / exit]: ").strip()

        if not user_input:
            continue

        cmd = user_input.lower()

        if cmd == 'exit':
            print("Exiting terminal screening tool. Goodbye!")
            sys.exit(0)

        elif cmd == 'list':
            elements = sorted(PERIODIC_METALS_DB.keys())
            print(f"\nAvailable Elements ({len(elements)} total):")
            for i in range(0, len(elements), 12):
                print("  " + ", ".join(elements[i:i+12]))
            print()
            continue

        elif cmd == 'preset':
            print("\nSelect Preset Alloy System:")
            for k, v in presets.items():
                print(f"  [{k}] {v[0]}")
            p_choice = input("Enter preset number: ").strip()
            if p_choice in presets:
                name, current_comp = presets[p_choice]
                print(f"\nLoaded Preset: {name}")
            else:
                print("Invalid preset selection.")
                continue

        else:
            # Parse custom element-percentage pairs
            try:
                tokens = user_input.replace(',', ' ').split()
                if len(tokens) % 2 != 0:
                    print("Error: Provide alternating element and percentage pairs (e.g., Zr 55 Cu 30).")
                    continue
                for i in range(0, len(tokens), 2):
                    el = tokens[i].capitalize()
                    pct = float(tokens[i+1])
                    if el not in PERIODIC_METALS_DB:
                        print(f"Warning: Element '{el}' not found in 76-element database. Skipping.")
                        continue
                    current_comp[el] = pct
            except Exception as e:
                print(f"Parsing error: {e}")
                continue

        if not current_comp:
            print("Error: Composition is empty.")
            continue

        # Normalize percentages to fractional moles for the analyzer
        total_pct = sum(current_comp.values())
        fraction_comp = {el: val / total_pct for el, val in current_comp.items()}
        formula = "".join([f"{el}{pct:g}" for el, pct in current_comp.items()])

        try:
            res = predict_glass_properties(
                fraction_comp, 
                alpha=0.000, 
                beta=0.384, 
                use_derived_moduli=True,
                enforce_gfa_gate=True
            )

            print("\n" + "=" * 54)
            print(f" RESULTS FOR ALLOY: {formula}")
            print("=" * 54)
            print(f"  Bulk Modulus (K_VRH):       {res['K_VRH_GPa']:.2f} GPa")
            print(f"  Shear Modulus (G_VRH):      {res['G_VRH_GPa']:.2f} GPa")
            print(f"  Poisson's Ratio (ν):        {res['Poisson_nu']:.3f}")
            print(f"  Enthalpy of Mixing (ΔH):    {res['Delta_H_mix_kJ']:.2f} kJ/mol")
            print(f"  Size Mismatch (δ):          {res['delta_pct']:.2f} %")
            print(f"  Vitrification Gate (δ≥8%):  {'PASS' if res['Vitrifies_GFA'] else 'CAUTION (<8%)'}")
            
            if res['Tg_Gated_K'] is not None:
                print(f"  Predicted Glass Temp (Tg):  {res['Tg_Gated_K']:.1f} K")
            else:
                print(f"  Predicted Glass Temp (Tg):  N/A (Crystallizes / δ < 8%)")
            print("=" * 54)
        except Exception as e:
            print(f"Calculation error: {e}")


if __name__ == "__main__":
    run_cli()