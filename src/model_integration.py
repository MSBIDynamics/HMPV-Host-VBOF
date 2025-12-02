"""
Model Integration Module
========================

This module handles the integration of the HMPV VBOF into host genome-scale 
metabolic models (GEMs) using COBRApy.

Key Functions:
- load_host_model: Load and validate the host metabolic model
- map_metabolites: Map VBOF metabolite IDs to host model IDs
- create_vbof_reaction: Create COBRApy reaction from VBOF stoichiometry
- integrate_vbof: Add VBOF reaction to host model
- validate_integrated_model: Validate the model after integration

Supported Host Models:
- iHsaEC21.xml (Human epithelial cell - primary target)
- Recon3D, Human1, iHsa (alternative models)

=============================================================================
SOURCES AND REFERENCES:
=============================================================================

1. COBRApy Documentation:
   https://cobrapy.readthedocs.io/

2. BiGG Models Database:
   http://bigg.ucsd.edu/
   King et al. (2016) Nucleic Acids Res. DOI: 10.1093/nar/gkv1049

3. VBOF Integration Methodology:
   - Aller et al. (2018) J R Soc Interface 15:20180442
   - Thiele et al. (2020) Bioinformatics

=============================================================================

Author: Syed Mushahid Hussain
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass

import cobra
from cobra import Model, Reaction, Metabolite
from cobra.io import read_sbml_model, write_sbml_model

from .exceptions import (
    ModelIntegrationError,
    IDMappingError,
    HMPVModelError
)
from .config import (
    MODEL_DIR,
    DEFAULT_HOST_MODEL,
    OUTPUT_DIR,
    VBOF_REACTION_ID,
    VBOF_REACTION_NAME,
    VBOF_SUBSYSTEM
)

logger = logging.getLogger(__name__)


@dataclass
class MetaboliteMapping:
    """Container for metabolite mapping results."""
    vbof_id: str
    model_id: str
    found: bool
    coefficient: float
    suggestions: List[str] = None


# =============================================================================
# METABOLITE ID MAPPINGS
# =============================================================================

# Common BiGG ID variations between models
# Format: {generic_id: [possible_model_ids]}
# NOTE: iHsaEC21 model uses specific conventions:
#   - Some use uppercase: GLN_c, ILE_c, TYR_c, TRP_c
#   - Some use _DASH_DASH_DASH_ instead of __: arg_DASH_DASH_DASH_L_c
#   - Some use _DASH_ for -: pchol_DASH_hs_c
METABOLITE_ID_VARIANTS = {
    # Nucleotides
    'atp_c': ['atp_c', 'ATP_c', 'M_atp_c'],
    'utp_c': ['utp_c', 'UTP_c', 'M_utp_c'],
    'gtp_c': ['gtp_c', 'GTP_c', 'M_gtp_c'],
    'ctp_c': ['ctp_c', 'CTP_c', 'M_ctp_c'],
    'adp_c': ['adp_c', 'ADP_c', 'M_adp_c'],
    'gdp_c': ['gdp_c', 'GDP_c', 'M_gdp_c'],
    'amp_c': ['amp_c', 'AMP_c', 'M_amp_c'],
    
    # Phosphate
    'pi_c': ['pi_c', 'Pi_c', 'M_pi_c', 'PI_c'],
    'ppi_c': ['ppi_c', 'PPi_c', 'M_ppi_c', 'PPI_c'],
    
    # Water and protons
    'h2o_c': ['h2o_c', 'H2O_c', 'M_h2o_c', 'WATER_c', 'h2o_m', 'h2o_n', 'h2o_r', 'h2o_l'],
    'h_c': ['h_c', 'H_c', 'M_h_c', 'PROTON_c'],
    
    # Amino acids (L-forms) - iHsaEC21 specific mappings
    'ala__L_c': ['ala__L_c', 'ala_L_c', 'ALA_c', 'L_45_ALPHA_45_ALANINE_c', 'ala_DASH_DASH_DASH_L_c'],
    'arg__L_c': ['arg__L_c', 'arg_L_c', 'ARG_c', 'arg_DASH_DASH_DASH_L_c', 'L_45_ARGININE_c'],
    'asn__L_c': ['asn__L_c', 'asn_L_c', 'ASN_c', 'asn_DASH_DASH_DASH_L_c', 'L_45_ASPARAGINE_c'],
    'asp__L_c': ['asp__L_c', 'asp_L_c', 'ASP_c', 'asp_DASH_DASH_DASH_L_c', 'L_45_ASPARTATE_c'],
    'cys__L_c': ['cys__L_c', 'cys_L_c', 'CYS_c', 'cys_DASH_DASH_DASH_L_c', 'L_45_CYSTEINE_c'],
    'glu__L_c': ['glu__L_c', 'glu_L_c', 'GLU_c', 'glu_DASH_DASH_DASH_L_c', 'L_45_GLUTAMATE_c'],
    'gln__L_c': ['gln__L_c', 'gln_L_c', 'GLN_c', 'gln_DASH_DASH_DASH_L_c', 'L_45_GLUTAMINE_c'],
    'gly_c': ['gly_c', 'GLY_c', 'Gly_c', 'GLYCINE_c'],
    'his__L_c': ['his__L_c', 'his_L_c', 'HIS_c', 'his_DASH_DASH_DASH_L_c', 'L_45_HISTIDINE_c'],
    'ile__L_c': ['ile__L_c', 'ile_L_c', 'ILE_c', 'ile_DASH_DASH_DASH_L_c', 'L_45_ISOLEUCINE_c'],
    'leu__L_c': ['leu__L_c', 'leu_L_c', 'LEU_c', 'leu_DASH_DASH_DASH_L_c', 'L_45_LEUCINE_c'],
    'lys__L_c': ['lys__L_c', 'lys_L_c', 'LYS_c', 'lys_DASH_DASH_DASH_L_c', 'L_45_LYSINE_c'],
    'met__L_c': ['met__L_c', 'met_L_c', 'MET_c', 'met_DASH_DASH_DASH_L_c', 'L_45_METHIONINE_c'],
    'phe__L_c': ['phe__L_c', 'phe_L_c', 'PHE_c', 'phe_DASH_DASH_DASH_L_c', 'L_45_PHENYLALANINE_c'],
    'pro__L_c': ['pro__L_c', 'pro_L_c', 'PRO_c', 'pro_DASH_DASH_DASH_L_c', 'L_45_PROLINE_c'],
    'ser__L_c': ['ser__L_c', 'ser_L_c', 'SER_c', 'ser_DASH_DASH_DASH_L_c', 'L_45_SERINE_c'],
    'thr__L_c': ['thr__L_c', 'thr_L_c', 'THR_c', 'thr_DASH_DASH_DASH_L_c', 'L_45_THREONINE_c'],
    'trp__L_c': ['trp__L_c', 'trp_L_c', 'TRP_c', 'trp_DASH_DASH_DASH_L_c', 'L_45_TRYPTOPHAN_c'],
    'tyr__L_c': ['tyr__L_c', 'tyr_L_c', 'TYR_c', 'tyr_DASH_DASH_DASH_L_c', 'L_45_TYROSINE_c'],
    'val__L_c': ['val__L_c', 'val_L_c', 'VAL_c', 'val_DASH_DASH_DASH_L_c', 'L_45_VALINE_c'],
    
    # Lipids - iHsaEC21 uses _DASH_ for hyphens
    'pc_hs_c': ['pc_hs_c', 'pchol_DASH_hs_c', 'pchol_hs_c', 'PHOSPHATIDYLCHOLINE_c'],
    'pe_hs_c': ['pe_hs_c', 'pe_DASH_hs_c', 'PHOSPHATIDYLETHANOLAMINE_c', 'pail_hs_c'],
    'ps_hs_c': ['ps_hs_c', 'ps_DASH_hs_c', 'PHOSPHATIDYLSERINE_c'],
    'sphmyln_hs_c': ['sphmyln_hs_c', 'sphmyln_DASH_hs_c', 'sphings_c', 'SPHINGOMYELIN_c'],
    'chsterol_c': ['chsterol_c', 'chol_c', 'cholesterol_c', 'CHOLESTEROL_c', 'CHSTEROL_c'],
    
    # Glycan precursors (nucleotide sugars)
    'uacgam_c': ['uacgam_c', 'udpacgam_c', 'udpglcnac_c', 'UDP_45_N_45_ACETYL_45_D_45_GLUCOSAMINE_c'],
    'gdpmann_c': ['gdpmann_c', 'gdpman_c', 'GDP_45_MANNOSE_c', 'GDP_45_D_45_MANNOSE_c'],
    'udpgal_c': ['udpgal_c', 'udpg_c', 'UDP_45_GALACTOSE_c', 'UDP_45_D_45_GALACTOSE_c'],
    'cmpacna_c': ['cmpacna_c', 'cmpneu5ac_c', 'CMP_45_N_45_ACETYLNEURAMINATE_c'],
    'gdpfuc_c': ['gdpfuc_c', 'gdpfuca_c', 'GDP_45_FUCOSE_c', 'GDP_45_L_45_FUCOSE_c'],
    'udpgalfur_c': ['udpgalfur_c', 'udpgalnac_c', 'UDP_45_N_45_ACETYL_45_D_45_GALACTOSAMINE_c'],
}


def load_host_model(model_path: str = None) -> Model:
    """
    Load the host metabolic model from SBML file.
    
    Parameters:
    -----------
    model_path : str, optional
        Path to the SBML model file. If not provided, uses default.
    
    Returns:
    --------
    cobra.Model : Loaded metabolic model
    
    Raises:
    -------
    ModelIntegrationError : If model cannot be loaded
    """
    if model_path is None:
        model_path = MODEL_DIR / DEFAULT_HOST_MODEL
    
    model_path = Path(model_path)
    
    if not model_path.exists():
        raise ModelIntegrationError(f"Model file not found: {model_path}")
    
    logger.info(f"Loading host model from: {model_path}")
    
    try:
        model = read_sbml_model(str(model_path))
        
        logger.info(f"Model loaded successfully:")
        logger.info(f"  Model ID: {model.id}")
        logger.info(f"  Reactions: {len(model.reactions)}")
        logger.info(f"  Metabolites: {len(model.metabolites)}")
        logger.info(f"  Genes: {len(model.genes)}")
        
        return model
    
    except Exception as e:
        raise ModelIntegrationError(f"Failed to load model: {e}")


def get_model_metabolite_ids(model: Model) -> Set[str]:
    """
    Get all metabolite IDs in the model.
    
    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    
    Returns:
    --------
    set : Set of metabolite IDs
    """
    return {met.id for met in model.metabolites}


def find_metabolite_in_model(
    model: Model,
    vbof_id: str,
    model_metabolite_ids: Set[str] = None
) -> Tuple[str, bool, List[str]]:
    """
    Find the corresponding metabolite ID in the model.
    
    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    vbof_id : str
        The VBOF metabolite ID
    model_metabolite_ids : set, optional
        Pre-computed set of model metabolite IDs
    
    Returns:
    --------
    tuple : (model_id, found, suggestions)
    """
    if model_metabolite_ids is None:
        model_metabolite_ids = get_model_metabolite_ids(model)
    
    # Direct match
    if vbof_id in model_metabolite_ids:
        return vbof_id, True, []
    
    # Check variants
    variants = METABOLITE_ID_VARIANTS.get(vbof_id, [])
    for variant in variants:
        if variant in model_metabolite_ids:
            return variant, True, []
    
    # Try with M_ prefix (common in SBML)
    m_prefix_id = f"M_{vbof_id}"
    if m_prefix_id in model_metabolite_ids:
        return m_prefix_id, True, []
    
    # Find similar IDs for suggestions
    base_id = vbof_id.split('_')[0]
    suggestions = [
        mid for mid in model_metabolite_ids
        if base_id.lower() in mid.lower()
    ][:5]
    
    return vbof_id, False, suggestions


def map_vbof_metabolites(
    model: Model,
    vbof_stoichiometry: Dict[str, float]
) -> Tuple[Dict[str, MetaboliteMapping], List[str]]:
    """
    Map all VBOF metabolites to model metabolites.
    
    Parameters:
    -----------
    model : cobra.Model
        The host metabolic model
    vbof_stoichiometry : dict
        VBOF stoichiometric coefficients
    
    Returns:
    --------
    tuple : (mappings_dict, unmapped_ids_list)
    """
    logger.info("Mapping VBOF metabolites to host model...")
    
    model_met_ids = get_model_metabolite_ids(model)
    mappings = {}
    unmapped = []
    
    for vbof_id, coef in vbof_stoichiometry.items():
        model_id, found, suggestions = find_metabolite_in_model(
            model, vbof_id, model_met_ids
        )
        
        mapping = MetaboliteMapping(
            vbof_id=vbof_id,
            model_id=model_id,
            found=found,
            coefficient=coef,
            suggestions=suggestions
        )
        mappings[vbof_id] = mapping
        
        if not found:
            unmapped.append(vbof_id)
            logger.warning(f"Metabolite not found: {vbof_id}")
            if suggestions:
                logger.warning(f"  Suggestions: {suggestions}")
    
    # Summary
    found_count = sum(1 for m in mappings.values() if m.found)
    logger.info(f"Mapping complete: {found_count}/{len(mappings)} metabolites found")
    
    if unmapped:
        logger.warning(f"Unmapped metabolites: {unmapped}")
    
    return mappings, unmapped


def create_vbof_reaction(
    model: Model,
    vbof_stoichiometry: Dict[str, float],
    reaction_id: str = VBOF_REACTION_ID,
    reaction_name: str = VBOF_REACTION_NAME,
    subsystem: str = VBOF_SUBSYSTEM,
    skip_unmapped: bool = False
) -> Tuple[Reaction, List[str]]:
    """
    Create a COBRApy reaction for the VBOF.
    
    Parameters:
    -----------
    model : cobra.Model
        The host metabolic model
    vbof_stoichiometry : dict
        VBOF stoichiometric coefficients
    reaction_id : str
        Reaction ID (default: 'HMPV_VBOF')
    reaction_name : str
        Reaction name
    subsystem : str
        Subsystem for the reaction
    skip_unmapped : bool
        If True, skip unmapped metabolites. If False, raise error.
    
    Returns:
    --------
    tuple : (Reaction, list of unmapped IDs)
    
    Raises:
    -------
    IDMappingError : If metabolites cannot be mapped and skip_unmapped=False
    """
    logger.info(f"Creating VBOF reaction: {reaction_id}")
    
    # Map metabolites
    mappings, unmapped = map_vbof_metabolites(model, vbof_stoichiometry)
    
    if unmapped and not skip_unmapped:
        raise IDMappingError(
            f"Cannot create VBOF: {len(unmapped)} metabolites not found in model",
            suggestions=unmapped
        )
    
    # Create reaction
    reaction = Reaction(reaction_id)
    reaction.name = reaction_name
    reaction.subsystem = subsystem
    reaction.lower_bound = 0  # Irreversible (production only)
    reaction.upper_bound = 1000
    
    # Add metabolites
    metabolites_to_add = {}
    for vbof_id, mapping in mappings.items():
        if mapping.found:
            met = model.metabolites.get_by_id(mapping.model_id)
            metabolites_to_add[met] = mapping.coefficient
        else:
            if skip_unmapped:
                logger.warning(f"Skipping unmapped metabolite: {vbof_id}")
    
    reaction.add_metabolites(metabolites_to_add)
    
    logger.info(f"VBOF reaction created:")
    logger.info(f"  Metabolites: {len(metabolites_to_add)}")
    logger.info(f"  Skipped (unmapped): {len(unmapped)}")
    
    return reaction, unmapped


def integrate_vbof(
    model: Model,
    vbof_stoichiometry: Dict[str, float],
    reaction_id: str = VBOF_REACTION_ID,
    skip_unmapped: bool = True
) -> Tuple[Model, Reaction, List[str]]:
    """
    Integrate the VBOF reaction into the host model.
    
    Parameters:
    -----------
    model : cobra.Model
        The host metabolic model
    vbof_stoichiometry : dict
        VBOF stoichiometric coefficients
    reaction_id : str
        Reaction ID for the VBOF
    skip_unmapped : bool
        If True, skip metabolites not found in model
    
    Returns:
    --------
    tuple : (modified_model, vbof_reaction, unmapped_metabolites)
    """
    logger.info("Integrating VBOF into host model...")
    
    # Check if reaction already exists
    if reaction_id in [r.id for r in model.reactions]:
        logger.warning(f"Reaction {reaction_id} already exists. Removing old version.")
        model.remove_reactions([reaction_id])
    
    # Create VBOF reaction
    vbof_reaction, unmapped = create_vbof_reaction(
        model=model,
        vbof_stoichiometry=vbof_stoichiometry,
        reaction_id=reaction_id,
        skip_unmapped=skip_unmapped
    )
    
    # Add reaction to model
    model.add_reactions([vbof_reaction])
    
    logger.info(f"VBOF integrated successfully!")
    logger.info(f"  Model now has {len(model.reactions)} reactions")
    
    return model, vbof_reaction, unmapped


def validate_integrated_model(
    model: Model,
    vbof_reaction_id: str = VBOF_REACTION_ID
) -> Dict:
    """
    Validate the integrated model.
    
    Parameters:
    -----------
    model : cobra.Model
        The integrated model
    vbof_reaction_id : str
        ID of the VBOF reaction
    
    Returns:
    --------
    dict : Validation results
    """
    logger.info("Validating integrated model...")
    
    results = {
        'valid': True,
        'errors': [],
        'warnings': [],
        'vbof_info': {},
        'fba_result': None
    }
    
    # Check VBOF reaction exists
    if vbof_reaction_id not in [r.id for r in model.reactions]:
        results['valid'] = False
        results['errors'].append(f"VBOF reaction {vbof_reaction_id} not found")
        return results
    
    vbof = model.reactions.get_by_id(vbof_reaction_id)
    results['vbof_info'] = {
        'id': vbof.id,
        'name': vbof.name,
        'num_metabolites': len(vbof.metabolites),
        'bounds': (vbof.lower_bound, vbof.upper_bound),
        'reaction_string': vbof.reaction[:200] + '...' if len(vbof.reaction) > 200 else vbof.reaction
    }
    
    # Try FBA with VBOF as objective
    try:
        with model:
            model.objective = vbof_reaction_id
            solution = model.optimize()
            
            results['fba_result'] = {
                'status': solution.status,
                'objective_value': solution.objective_value if solution.status == 'optimal' else None,
                'feasible': solution.status == 'optimal'
            }
            
            if solution.status != 'optimal':
                results['warnings'].append(f"FBA with VBOF objective: {solution.status}")
            else:
                logger.info(f"FBA successful! VBOF flux: {solution.objective_value:.6f}")
                
    except Exception as e:
        results['errors'].append(f"FBA failed: {str(e)}")
        results['valid'] = False
    
    # Check for unbounded metabolites in VBOF
    for met, coef in vbof.metabolites.items():
        # Check if metabolite has producing/consuming reactions
        producing = [r for r in met.reactions if r.id != vbof_reaction_id and r.metabolites[met] > 0]
        consuming = [r for r in met.reactions if r.id != vbof_reaction_id and r.metabolites[met] < 0]
        
        if coef < 0 and not producing:
            results['warnings'].append(f"VBOF consumes {met.id} but no producing reactions exist")
        if coef > 0 and not consuming:
            results['warnings'].append(f"VBOF produces {met.id} but no consuming reactions exist")
    
    logger.info(f"Validation complete:")
    logger.info(f"  Valid: {results['valid']}")
    logger.info(f"  Errors: {len(results['errors'])}")
    logger.info(f"  Warnings: {len(results['warnings'])}")
    
    return results


def save_integrated_model(
    model: Model,
    output_path: str = None,
    format: str = 'sbml'
) -> str:
    """
    Save the integrated model to file.
    
    Parameters:
    -----------
    model : cobra.Model
        The integrated model
    output_path : str, optional
        Output file path. If not provided, uses default.
    format : str
        Output format ('sbml' or 'json')
    
    Returns:
    --------
    str : Path to saved file
    """
    if output_path is None:
        output_dir = Path(OUTPUT_DIR)
        output_dir.mkdir(exist_ok=True)
        
        if format == 'sbml':
            output_path = output_dir / f"{model.id}_with_HMPV_VBOF.xml"
        else:
            output_path = output_dir / f"{model.id}_with_HMPV_VBOF.json"
    
    output_path = Path(output_path)
    
    logger.info(f"Saving integrated model to: {output_path}")
    
    if format == 'sbml':
        write_sbml_model(model, str(output_path))
    else:
        cobra.io.save_json_model(model, str(output_path))
    
    logger.info(f"Model saved successfully!")
    
    return str(output_path)


def get_integration_summary(
    model: Model,
    vbof_reaction_id: str,
    unmapped: List[str],
    validation: Dict
) -> str:
    """
    Generate a summary of the integration.
    
    Parameters:
    -----------
    model : cobra.Model
        The integrated model
    vbof_reaction_id : str
        ID of the VBOF reaction
    unmapped : list
        List of unmapped metabolites
    validation : dict
        Validation results
    
    Returns:
    --------
    str : Summary text
    """
    vbof = model.reactions.get_by_id(vbof_reaction_id)
    
    summary = f"""
================================================================================
HMPV VBOF INTEGRATION SUMMARY
================================================================================

MODEL INFORMATION:
------------------
Model ID: {model.id}
Total Reactions: {len(model.reactions)}
Total Metabolites: {len(model.metabolites)}
Total Genes: {len(model.genes)}

VBOF REACTION:
--------------
Reaction ID: {vbof.id}
Reaction Name: {vbof.name}
Subsystem: {vbof.subsystem}
Metabolites in VBOF: {len(vbof.metabolites)}
Bounds: [{vbof.lower_bound}, {vbof.upper_bound}]

METABOLITE MAPPING:
-------------------
Mapped metabolites: {len(vbof.metabolites)}
Unmapped metabolites: {len(unmapped)}
"""
    
    if unmapped:
        summary += f"\nUnmapped IDs: {', '.join(unmapped)}\n"
    
    summary += f"""
VALIDATION:
-----------
Valid: {validation['valid']}
Errors: {len(validation['errors'])}
Warnings: {len(validation['warnings'])}
"""
    
    if validation['fba_result']:
        summary += f"""
FBA RESULT:
-----------
Status: {validation['fba_result']['status']}
Objective Value: {validation['fba_result']['objective_value']}
Feasible: {validation['fba_result']['feasible']}
"""
    
    summary += "================================================================================\n"
    
    return summary

