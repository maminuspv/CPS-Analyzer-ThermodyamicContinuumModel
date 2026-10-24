"""
cps_analyzer
------------
First-Principles Continuum-Thermodynamic Screening Engine for Refractory Bulk Metallic Glasses.
"""

from .affinity import derive_affinity_and_moduli, PERIODIC_METALS_DB
from .continuum import calculate_vrh_elasticity
from .thermodynamics import calculate_mixing_enthalpy, get_miedema_pairwise_enthalpy
from .kinetics import predict_glass_properties, GFA_DELTA_THRESHOLD

__version__ = "1.0.0"
__author__ = "Matthew B. Minus"

__all__ = [
    "derive_affinity_and_moduli",
    "calculate_vrh_elasticity",
    "calculate_mixing_enthalpy",
    "get_miedema_pairwise_enthalpy",
    "predict_glass_properties",
    "PERIODIC_METALS_DB",
    "GFA_DELTA_THRESHOLD"
]