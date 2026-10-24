"""
tests/test_physics_engine.py
----------------------------
Automated unit tests for the CPS Refractory Screening physics engine.
Asserts pure element limits, Voigt-Reuss-Hill bounds, Miedema enthalpies,
and BMG glass transition (Tg) regression invariants.
"""

import pytest
import numpy as np
from cps_analyzer.affinity import derive_affinity_and_moduli, PERIODIC_METALS_DB
from cps_analyzer.continuum import calculate_vrh_elasticity
from cps_analyzer.thermodynamics import calculate_mixing_enthalpy, get_miedema_pairwise_enthalpy
from cps_analyzer.kinetics import predict_glass_properties


def test_pure_element_database():
    """Asserts that essential elements exist in the 76-element database."""
    assert "Zr" in PERIODIC_METALS_DB
    assert "Ta" in PERIODIC_METALS_DB
    assert "Ti" in PERIODIC_METALS_DB
    assert "Ni" in PERIODIC_METALS_DB


def test_affinity_and_moduli_derivation():
    """Asserts that Rose UBER derivations yield physically realistic bulk moduli for Ti."""
    props = derive_affinity_and_moduli("Ti")
    assert props["r_pm"] == 147.0
    assert props["K_der_GPa"] > 50.0  # Should reflect refractory stiffness scale
    assert props["G_der_GPa"] > 20.0


def test_miedema_pairwise_enthalpy():
    """Asserts that known exothermic binary pairs return negative mixing enthalpies."""
    # Zr-Ni has a deep exothermic mixing enthalpy (-49.0 kJ/mol)
    h_zr_ni = get_miedema_pairwise_enthalpy("Zr", "Ni")
    assert h_zr_ni < 0.0
    assert h_zr_ni == -49.0


def test_vrh_homogenization_pure_limit():
    """Asserts that a 100% single-element composition matches its reference moduli."""
    comp = {"Ti": 1.0}
    res = calculate_vrh_elasticity(comp, use_derived_moduli=False)
    assert pytest.approx(res["K_VRH"], rel=0.01) == PERIODIC_METALS_DB["Ti"][4]
    assert pytest.approx(res["G_VRH"], rel=0.01) == PERIODIC_METALS_DB["Ti"][5]


def test_vitreloy_bmg_regression():
    """Asserts that Vitreloy (Zr55Cu30Al10Ni5) passes GFA thresholds and predicts realistic Tg."""
    composition = {"Zr": 0.55, "Cu": 0.30, "Al": 0.10, "Ni": 0.05}
    res = predict_glass_properties(composition, alpha=0.000, beta=0.384, enforce_gfa_gate=True)
    
    # Vitreloy should vitrify (delta >= 8.0%)
    assert res["Vitrifies_GFA"] is True
    assert res["delta_pct"] >= 8.0
    
    # Predicted Tg should align closely with experimental Vitreloy (~685 K)
    assert res["Tg_Calibrated_K"] is not None
    assert 600.0 <= res["Tg_Calibrated_K"] <= 780.0


def test_quaternary_fusion_armor():
    """Asserts that the optimal fusion armor candidate (Zr5Ta55Ni35Ti5) meets mechanical floors."""
    composition = {"Zr": 0.05, "Ta": 0.55, "Ni": 0.35, "Ti": 0.05}
    res = predict_glass_properties(composition, alpha=0.000, beta=0.384)
    
    # Must meet fusion divertor thresholds: G >= 45 GPa, Tg >= 700 K
    assert res["G_VRH_GPa"] >= 45.0
    assert res["Tg_Calibrated_K"] >= 700.0
    assert res["Poisson_nu"] > 0.31  # Ductile threshold check