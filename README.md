# High-Throughput Refractory BMG Screening Pipeline

An open-source Python analytical engine and high-throughput combinatorial screening framework developed for rapid discovery of high-temperature bulk metallic glasses (BMGs) tailored for nuclear fusion **Capillary Porous Systems (CPS)** divertor armor.

---

## 🚀 Theoretical Foundations

This framework couples first-principles electronic parameters directly to macroscopic thermodynamic and structural flow barriers without costly density functional theory (DFT) or molecular dynamics (MD) computations:
* **Screened Nuclear Ion Affinity ($\mathcal{A}_i$):** Derives elemental cohesive strength and electrostatic well-depths directly from core effective charge ($Z_{\text{eff}}$) and Wigner-Seitz metallic radius ($r_i$) modulated by Thomas-Fermi conduction screening ($f_{\text{screen}} \approx 0.145$).
* **Rose Universal Binding Energy Relation (UBER):** Computes intrinsic elemental bulk ($K_i$) and shear ($G_i$) moduli from volume derivatives of atomic binding potentials.
* **Voigt–Reuss–Hill (VRH) Homogenization:** Evaluates effective composite elastic matrices ($K_{\text{VRH}}, G_{\text{VRH}}, \nu$).
* **Topological Frustration ($\delta$):** Monitors structural packaging confusion and vitrification preconditions via Egami size polydispersity ($\delta \ge 8.0\%$).
* **Cooperative Shear Model (CSM):** Predicts glass transition temperatures ($T_g$) via Johnson-Samwer shear arrest scaling, incorporating optimized core valence disparity corrections ($\beta = 0.384$).

---

## 🛠️ Performance Features
* Evaluates complete physical property profiles in **< 0.001 ms per composition**.
* Scans complete multicomponent spaces (e.g., 839 candidate structures inside the quaternary $\text{Zr--Ta--Ni--Ti}$ grid matrix) in **< 1 second**.

---

## 📁 Repository Structure

```text
cps-refractory-screening/
├── cps_analyzer/                 # Core modular physics and thermodynamic engine
│   ├── __init__.py
│   ├── affinity.py             # Thomas-Fermi screening & Rose UBER derivations (A_i)
│   ├── continuum.py            # Voigt-Reuss-Hill homogenization (K_VRH, G_VRH, nu)
│   ├── thermodynamics.py       # Miedema pairwise Delta_H_mix summation
│   └── kinetics.py             # Johnson-Samwer CSM & Egami delta polydispersity (Tg)
├── data/                         # Compiled reference, calibration, and database files
│   ├── bmg_calibration_experimental_data.csv  # Experimental BMG training data
│   └── Table_S1_76_Elements_Database.csv      # Complete 76-element first-principles table
├── cli/                          # Interactive front-end tools
│   └── terminal_frontend.py     # Command-line interface for custom alloy screening
├── scripts/                      # Production execution and simulation scripts
│   ├── run_calibration.py        # SciPy optimizer for alpha and beta fitting
│   └── run_high_throughput_sweep.py # Quaternary Zr-Ta-Ni-Ti grid sweep
├── manuscript/                   # LaTeX source files, figures, and preprint PDF
│   ├── manuscript_cps_refractory.tex
│   ├── references.bib
│   └── chemrxiv.15010896_v1.pdf
├── tests/                        # Automated unit testing (regression invariants)
│   └── test_physics_engine.py    # Asserts pure element & benchmark limits
├── LICENSE                       # Apache 2.0 Software License
├── CITATION.cff                  # Machine-readable citation metadata linking to preprint
└── README.md                     # Project overview and quick-start guide