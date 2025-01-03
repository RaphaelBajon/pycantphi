from dataclasses import dataclass
from typing import List, Dict
import numpy as np


@dataclass(frozen=True)
class ChemicalConstants:
    ATMOSPHERIC_PRESSURE: float = 1.0
    PREINDUSTRIAL_CO2_MOLAR_FRACTION: float = 278.0
    K1K2: int = 4  # previously 10: Lueker, 2000. according best practices guide in 2007 - NOTE: Previously it was option 4, MEHRBACH refit BY DICKSON AND MILLERO.
    KSO4: int = 1  # Dickson's KSO4
    RN: float = 9.0  # N:C ratio
    RP: float = 135.0  # P:C ratio
    MAUNA_LOA_OFFSET_YEARS: int = 5

@dataclass(frozen=True)
class ProcessingConstants:
    OUTPUT_YEARS: List[int] = (1994, 2000, 2004, 2005, 2010)
    base_params_to_return = [
        'cAntPhiCt0ML', 'cAntPhiCt0ML1994', 'cAntPhiCt0ML2004',
        'cAntTroca', 'cAntTrocaML1994', 'cAntTrocaML2004',
        'cAntSatML', 'cAntSatML1994', 'cAntSatML2004'
    ]
    additional_params = [
        'x1', 'x2', 'x3', 'x4', 'x5', 'x6',
        'percentMW', 'salinityCorrectedByMW', 'thetaCorrectedByMW',
        'noConservativeCorrectedByMW', 'poConservativeCorrectedByMW',
        'alkalinity_preformed_parametric', 'alkalinity_preformed_parametric_revised',
        'deltaCDiseq', 'deltaCa', 'silicate_teorical', 'cStar'
    ]

@dataclass(frozen=True)
class MethodConstants:
    REQUIRED_PARAMETERS: List[str] = (
        'longitude', 'latitude', 'pressure', 'theta',
        'salinity', 'oxygen', 'silicate', 'nitrate', 'phosphate',
        'carbon', 'alkalinity', 'year'
    )
