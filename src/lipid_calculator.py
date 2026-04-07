"""
Lipid Calculator Module
=======================

Calculates lipid envelope requirements for HMPV virions based on
virion geometry and paramyxovirus lipid composition data.

Key Functions:
- calculate_lipid_stoichiometry: Calculate envelope lipid requirements

Author: Syed Mushahid Hussain
"""

import logging
from typing import Dict, Optional

from .config import LIPID_PACKING_DENSITY_NM2, LIPID_FRACTIONS
from .protein_analyzer import calculate_virion_surface_area

logger = logging.getLogger(__name__)


def calculate_lipid_stoichiometry(
    virion_diameter_nm: float = 209,
    include_lipids: bool = True,
    morphology: str = "spherical",
    virion_length_nm: Optional[float] = None,
) -> Dict[str, float]:
    """
    Calculate lipid requirements for the viral envelope.

    For HMPV (enveloped paramyxovirus):
    - Envelope derived from host plasma membrane
    - Typical composition similar to other paramyxoviruses

    Parameters:
    -----------
    virion_diameter_nm : float
        Average virion diameter in nanometers (default: 150 nm for HMPV)
    include_lipids : bool
        Whether to include lipid requirements (default: True)
    morphology : str
        Virion morphology ("spherical" or "filamentous")
    virion_length_nm : float, optional
        Virion length for filamentous morphology

    Returns:
    --------
    dict : Lipid stoichiometry (negative = consumption)

    Note:
    -----
    HMPV can be pleomorphic (variable shape), so this is an approximation
    based on spherical geometry.

    Lipid composition based on paramyxovirus literature:
    - Phosphatidylcholine (PC): ~45%
    - Phosphatidylethanolamine (PE): ~25%
    - Phosphatidylserine (PS): ~5%
    - Sphingomyelin (SM): ~15%
    - Cholesterol: ~10%

    Sources:
    --------
    1. Jinjun Shan. (2018) "High-resolution lipidomics reveals dysregulation of lipid metabolism"
       https://pubs.rsc.org/en/content/articlelanding/2018/ra/c8ra05640d

    2. Brügger et al. (2006) "The HIV lipidome"
       https://pmc.ncbi.nlm.nih.gov/articles/PMC1413831/

    3. Nagle & Tristram-Nagle (2000) "Structure of lipid bilayers"
       Biochim Biophys Acta 1469(3):159-195. DOI: 10.1016/S0304-4157(00)00016-2
       https://www.sciencedirect.com/science/article/pii/S0304415700000162
       - Lipid packing density: ~0.65 nm² per lipid molecule

    4. Aller et al. (2018) "Integrated human-virus metabolic stoichiometric modelling"
       https://pmc.ncbi.nlm.nih.gov/articles/PMC6170780/
       - VBOF methodology for enveloped viruses

    5. Characterization of the phospholipid and fatty acid composition of Sendai virus.
	Barnes, J A and Pehowich, D J and Allen, T M. (1987) Journal of Lipid Research 28(2):130-137. DOI: 10.1016/S0022-2275(20)38714-9
	https://www.sciencedirect.com/science/article/pii/S0022227520387149
	- Sendai virus lipid composition
	- Phosphatidylcholine: 18.65%
	- Phosphatidylethanolamine: 13.40%
	- Phosphatidylserine: 6.00%
	- Sphingomyelin: 4.40%
	- Cardiolipin: 4.50%
	- Phosphatidylinositol: 2.75%
	- Lysophosphatidylethanolamine: 0.80%
	- Cholesterol: 50.00%
    """
    if not include_lipids:
        logger.info("Lipid stoichiometry skipped (include_lipids=False)")
        return {}

    # Surface area — morphology-aware (spherical or filamentous capsule)
    surface_area_nm2 = calculate_virion_surface_area(
        morphology=morphology,
        diameter_nm=virion_diameter_nm,
        length_nm=virion_length_nm,
    )

    # Lipid packing density — controlled via config.LIPID_PACKING_DENSITY_NM2
    lipid_area_nm2 = LIPID_PACKING_DENSITY_NM2

    # Total lipids in bilayer (outer + inner leaflet)
    total_lipids = int((surface_area_nm2 / lipid_area_nm2) * 2)

    # Lipid composition fractions — controlled via config.LIPID_FRACTIONS
    lipid_fractions = LIPID_FRACTIONS

    stoichiometry = {}
    for lipid_id, fraction in lipid_fractions.items():
        count = int(total_lipids * fraction)
        if count > 0:
            stoichiometry[lipid_id] = -count

    logger.info(f"Lipid stoichiometry calculated ({morphology}):")
    logger.info(f"  Virion diameter: {virion_diameter_nm} nm"
                + (f"  length: {virion_length_nm} nm" if virion_length_nm else ""))
    logger.info(f"  Surface area: {surface_area_nm2:.0f} nm²")
    logger.info(f"  Total lipids: {total_lipids}")
    for lipid_id, count in stoichiometry.items():
        logger.info(f"  {lipid_id}: {count}")

    return stoichiometry
