"""
Protein Analyzer Module
=======================

This module provides functions to parse and analyze HMPV protein sequences,
calculating amino acid composition for VBOF construction.

Key Functions:
- load_proteins: Load protein sequences from FASTA file
- count_amino_acids: Count amino acids in a protein sequence
- calculate_protein_stoichiometry: Calculate amino acid requirements for VBOF

Author: Syed Mushahid Hussain
"""

import logging
import math
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, field

from .exceptions import MissingProteinDataError, MissingCopyNumberError
from .config import AMINO_ACID_BIGG_IDS, METABOLITE_MOLECULAR_WEIGHTS

logger = logging.getLogger(__name__)


# =============================================================================
# AMINO ACID PROPERTIES
# =============================================================================
#
# Standard amino acid molecular weights and BiGG metabolite IDs
#
# SOURCES AND REFERENCES:
#
#
# 1. BIGG DATABASE IDs:
#    - King et al. (2016) "BiGG Models: A platform for integrating genome-scale models"
#      Nucleic Acids Res 44(D1):D515-22. DOI: 10.1093/nar/gkv1049
#      http://bigg.ucsd.edu/universal/metabolites
#
# 2. AMINO ACID BIOCHEMISTRY:
#    - Berg, Tymoczko, Stryer (2015) "Biochemistry" 8th Edition
#      Chapter 2: Protein Composition and Structure
#      https://dn790008.ca.archive.org/0/items/JeremyM.BergJohnL.TymoczkoGregoryJ.GattoJr.LubertStryerBiochemistry_201708/Jeremy%20M.%20Berg%2C%20John%20L.%20Tymoczko%2C%20Gregory%20J.%20Gatto%20Jr.%2C%20Lubert%20Stryer%20Biochemistry.pdf
#
# BiGG IDs and MWs sourced from config — edit there to propagate everywhere
_AA_META = [
    ('A', 'Alanine',       'Ala'),
    ('R', 'Arginine',      'Arg'),
    ('N', 'Asparagine',    'Asn'),
    ('D', 'Aspartate',     'Asp'),
    ('C', 'Cysteine',      'Cys'),
    ('E', 'Glutamate',     'Glu'),
    ('Q', 'Glutamine',     'Gln'),
    ('G', 'Glycine',       'Gly'),
    ('H', 'Histidine',     'His'),
    ('I', 'Isoleucine',    'Ile'),
    ('L', 'Leucine',       'Leu'),
    ('K', 'Lysine',        'Lys'),
    ('M', 'Methionine',    'Met'),
    ('F', 'Phenylalanine', 'Phe'),
    ('P', 'Proline',       'Pro'),
    ('S', 'Serine',        'Ser'),
    ('T', 'Threonine',     'Thr'),
    ('W', 'Tryptophan',    'Trp'),
    ('Y', 'Tyrosine',      'Tyr'),
    ('V', 'Valine',        'Val'),
]
AMINO_ACID_INFO = {
    code: {
        'name':    name,
        'abbrev':  abbrev,
        'bigg_id': AMINO_ACID_BIGG_IDS[code],
        'mw':      METABOLITE_MOLECULAR_WEIGHTS[AMINO_ACID_BIGG_IDS[code]],
    }
    for code, name, abbrev in _AA_META
}


@dataclass
class ProteinInfo:
    """Container for protein information."""
    protein_id: str
    gene_name: str
    description: str
    sequence: str
    length: int
    amino_acid_counts: Dict[str, int] = field(default_factory=dict)
    molecular_weight: float = 0.0
    copy_number: Optional[int] = None
    
    def __post_init__(self):
        """Calculate derived properties after initialization."""
        if not self.amino_acid_counts:
            self.amino_acid_counts = count_amino_acids(self.sequence)
        if self.molecular_weight == 0.0:
            self.molecular_weight = calculate_molecular_weight(self.sequence)


def load_proteins(fasta_path: str) -> Dict[str, ProteinInfo]:
    """
    Load and parse protein sequences from FASTA file.
    
    Parameters:
    -----------
    fasta_path : str
        Path to the protein FASTA file (.faa)
    
    Returns:
    --------
    dict : Dictionary mapping gene names to ProteinInfo objects
    
    Raises:
    -------
    MissingProteinDataError : If file not found or invalid format
    
    Example:
    --------
    >>> proteins = load_proteins("Data/protein/GCF_002815375.1_ASM281537v1_protein.faa")
    >>> print(f"Loaded {len(proteins)} proteins")
    """
    filepath = Path(fasta_path)
    
    if not filepath.exists():
        raise MissingProteinDataError(f"Protein file not found: {fasta_path}")
    
    logger.info(f"Loading proteins from: {fasta_path}")
    
    proteins = {}
    current_header = None
    current_sequence = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('>'):
                # Save previous protein if exists
                if current_header is not None:
                    protein_info = _parse_protein_entry(current_header, ''.join(current_sequence))
                    proteins[protein_info.gene_name] = protein_info
                
                # Start new protein
                current_header = line[1:]  # Remove '>'
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # Don't forget the last protein
        if current_header is not None:
            protein_info = _parse_protein_entry(current_header, ''.join(current_sequence))
            proteins[protein_info.gene_name] = protein_info
    
    if not proteins:
        raise MissingProteinDataError("No proteins found in file")
    
    logger.info(f"Loaded {len(proteins)} proteins:")
    for gene_name, protein in proteins.items():
        logger.info(f"  {gene_name}: {protein.length} aa, MW={protein.molecular_weight:.2f} Da")
    
    return proteins


def _parse_protein_entry(header: str, sequence: str) -> ProteinInfo:
    """
    Parse a single protein entry from FASTA.
    
    Parameters:
    -----------
    header : str
        FASTA header line (without '>')
    sequence : str
        Amino acid sequence
    
    Returns:
    --------
    ProteinInfo : Parsed protein information
    """
    # Parse header: >YP_009513265.1 nucleoprotein [human metapneumovirus]
    parts = header.split(' ', 1)
    protein_id = parts[0]
    description = parts[1] if len(parts) > 1 else ""
    
    # Extract gene name from description
    # Common patterns: "nucleoprotein", "phosphoprotein", "matrix protein", etc.
    gene_name = _extract_gene_name(description, protein_id)
    
    sequence = sequence.upper().replace('*', '')  # Remove stop codon marker
    
    return ProteinInfo(
        protein_id=protein_id,
        gene_name=gene_name,
        description=description,
        sequence=sequence,
        length=len(sequence)
    )


def _extract_gene_name(description: str, protein_id: str) -> str:
    """
    Extract standardized gene name from protein description.
    
    Parameters:
    -----------
    description : str
        Protein description from FASTA header
    protein_id : str
        Protein accession ID
    
    Returns:
    --------
    str : Standardized gene name (N, P, M, F, M2-1, M2-2, SH, G, L)
    """
    description_lower = description.lower()
    
    # HMPV protein name mapping
    if 'nucleoprotein' in description_lower:
        return 'N'
    elif 'phosphoprotein' in description_lower:
        return 'P'
    elif 'matrix protein 2-1' in description_lower or 'matrix protein m2-1' in description_lower:
        return 'M2-1'
    elif 'matrix protein 2-2' in description_lower or 'matrix protein m2-2' in description_lower:
        return 'M2-2'
    elif 'matrix protein' in description_lower or 'matrix' in description_lower:
        return 'M'
    elif 'fusion' in description_lower:
        return 'F'
    elif 'small hydrophobic' in description_lower:
        return 'SH'
    elif 'attachment' in description_lower or 'glycoprotein g' in description_lower:
        return 'G'
    elif 'polymerase' in description_lower or 'rna-dependent' in description_lower:
        return 'L'
    else:
        # Return protein_id if can't determine gene name
        return protein_id


def count_amino_acids(sequence: str) -> Dict[str, int]:
    """
    Count amino acids in a protein sequence.
    
    Parameters:
    -----------
    sequence : str
        Amino acid sequence (single-letter code)
    
    Returns:
    --------
    dict : Count for each amino acid
    """
    sequence = sequence.upper()
    
    counts = {aa: 0 for aa in AMINO_ACID_INFO.keys()}
    
    for aa in sequence:
        if aa in counts:
            counts[aa] += 1
        elif aa not in ['X', '*', '-']:  # Ignore unknown, stop, gap
            logger.warning(f"Unknown amino acid: {aa}")
    
    return counts


def calculate_molecular_weight(sequence: str) -> float:
    """
    Calculate molecular weight of a protein.
    
    Parameters:
    -----------
    sequence : str
        Amino acid sequence
    
    Returns:
    --------
    float : Molecular weight in Daltons
    
    Note:
    -----
    Uses average molecular weights. Subtracts water for peptide bonds.
    """
    sequence = sequence.upper()
    water_mw = 18.015  # Water released per peptide bond
    
    mw = 0.0
    for aa in sequence:
        if aa in AMINO_ACID_INFO:
            mw += AMINO_ACID_INFO[aa]['mw']
    
    # Subtract water for peptide bonds (n-1 bonds for n amino acids)
    if len(sequence) > 1:
        mw -= water_mw * (len(sequence) - 1)
    
    return mw


def set_copy_numbers(
    proteins: Dict[str, ProteinInfo],
    copy_numbers: Dict[str, int]
) -> Dict[str, ProteinInfo]:
    """
    Set copy numbers for each protein.
    
    Parameters:
    -----------
    proteins : dict
        Dictionary of ProteinInfo objects
    copy_numbers : dict
        Dictionary mapping gene names to copy numbers per virion
    
    Returns:
    --------
    dict : Updated proteins dictionary
    
    Raises:
    -------
    MissingCopyNumberError : If copy number missing for any protein
    """
    for gene_name, protein in proteins.items():
        if gene_name in copy_numbers:
            protein.copy_number = copy_numbers[gene_name]
            logger.info(f"Set copy number for {gene_name}: {protein.copy_number}")
        else:
            raise MissingCopyNumberError(
                f"Copy number not provided for protein",
                protein_name=gene_name
            )
    
    return proteins


def calculate_protein_stoichiometry(
    proteins: Dict[str, ProteinInfo],
    copy_numbers: Optional[Dict[str, int]] = None
) -> Dict[str, float]:
    """
    Calculate total amino acid requirements for all viral proteins.
    
    Parameters:
    -----------
    proteins : dict
        Dictionary of ProteinInfo objects
    copy_numbers : dict, optional
        Copy numbers per virion. If not provided, uses values from ProteinInfo.
    
    Returns:
    --------
    dict : Stoichiometric coefficients for amino acids (BiGG IDs)
           Negative values indicate consumption
    
    Note:
    -----
    If copy numbers are not provided, raises MissingCopyNumberError.
    """
    # Initialize total counts
    total_aa_counts = {aa: 0 for aa in AMINO_ACID_INFO.keys()}
    
    for gene_name, protein in proteins.items():
        # Get copy number
        if copy_numbers and gene_name in copy_numbers:
            copy_num = copy_numbers[gene_name]
        elif protein.copy_number is not None:
            copy_num = protein.copy_number
        else:
            raise MissingCopyNumberError(
                "Copy number required for stoichiometry calculation",
                protein_name=gene_name
            )
        
        # Add weighted amino acid counts
        for aa, count in protein.amino_acid_counts.items():
            total_aa_counts[aa] += count * copy_num
    
    # Convert to BiGG IDs with negative coefficients (consumption)
    stoichiometry = {}
    for aa, count in total_aa_counts.items():
        if count > 0:
            bigg_id = AMINO_ACID_INFO[aa]['bigg_id']
            stoichiometry[bigg_id] = -count
    
    # Log summary
    total_aa = sum(total_aa_counts.values())
    logger.info(f"Protein stoichiometry calculated:")
    logger.info(f"  Total amino acids required: {total_aa}")
    logger.info(f"  Most abundant: {max(total_aa_counts, key=total_aa_counts.get)} "
                f"({max(total_aa_counts.values())})")
    
    return stoichiometry


def get_protein_summary(proteins: Dict[str, ProteinInfo]) -> Dict:
    """
    Get a summary of all proteins.
    
    Parameters:
    -----------
    proteins : dict
        Dictionary of ProteinInfo objects
    
    Returns:
    --------
    dict : Summary statistics for each protein
    """
    summary = {
        'total_proteins': len(proteins),
        'proteins': {}
    }
    
    total_aa = 0
    total_mw = 0.0
    
    for gene_name, protein in proteins.items():
        summary['proteins'][gene_name] = {
            'protein_id': protein.protein_id,
            'length': protein.length,
            'molecular_weight': round(protein.molecular_weight, 2),
            'copy_number': protein.copy_number,
            'description': protein.description
        }
        total_aa += protein.length
        total_mw += protein.molecular_weight
    
    summary['total_amino_acids'] = total_aa
    summary['total_molecular_weight'] = round(total_mw, 2)

    return summary


# =============================================================================
# BIOPHYSICAL COPY NUMBER CALCULATIONS
# =============================================================================

def calculate_virion_surface_area(
    morphology: str,
    diameter_nm: float,
    length_nm: Optional[float] = None,
) -> float:
    """
    Calculate virion outer surface area based on particle morphology.

    Spherical
    ---------
    Models the particle as a perfect sphere:
        A = 4π r²
    where r = diameter_nm / 2.

    Filamentous
    -----------
    Models the particle as a capsule (straight cylinder capped at both ends
    with hemispheres of the same radius):
        A = 2π r (L - 2r) + 4π r²   →   simplifies to   A = 2π r L
    where r = diameter_nm / 2  and  L = total tip-to-tip length_nm.

    Derivation of the simplification:
      Lateral cylinder area  = 2π r (L - 2r)
      Two hemispherical caps = 4π r²  (= one full sphere)
      Total = 2πrL - 4πr² + 4πr² = 2πrL  

    Note: when L = 2r (degenerate case, a sphere), A = 2πr·2r = 4πr² 

    Parameters
    ----------
    morphology : str
        "spherical" or "filamentous".
    diameter_nm : float
        Sphere diameter for spherical; tube cross-section diameter for
        filamentous (nm).
    length_nm : float | None
        Total tip-to-tip length for filamentous particles (nm).
        Must be >= diameter_nm.  Ignored for spherical.

    Returns
    -------
    float  Surface area in nm².

   
    """
    radius = diameter_nm / 2.0
    if morphology == "spherical":
        return 4.0 * math.pi * radius ** 2
    elif morphology == "filamentous":
        if length_nm is None:
            raise ValueError(
                "length_nm must be provided for filamentous morphology."
            )
        if length_nm < diameter_nm:
            raise ValueError(
                f"length_nm ({length_nm} nm) must be >= diameter_nm ({diameter_nm} nm) "
                "for a valid capsule shape."
            )
        return 2.0 * math.pi * radius * length_nm
    else:
        raise ValueError(
            f"Unknown morphology {morphology!r}. Expected 'spherical' or 'filamentous'."
        )


def calculate_n_copy_number(
    genome_length: int,
    nt_per_protomer: int = 7,
    genome_copies: int = 1
) -> int:
    """
    Calculate N protein copy number from genome length and binding stoichiometry.

    Formula: N_copies = (L_genome * genome_copies) / nt_per_protomer

    Parameters
    ----------
    genome_length : int
        Total nucleotides in the viral genome.
    nt_per_protomer : int
        Nucleotides bound per N protomer (default 7, from Tawar et al. 2009).
    genome_copies : int
        Number of genome copies per virion (default 1 for spherical particles).

    Returns
    -------
    int : Estimated N protein copy number.

    References
    ----------
    Tawar et al. (2009) Science 326:1279-83  — RSV N-RNA crystal structure, 7 nt/N.
    """
    copies = (genome_length * genome_copies) / nt_per_protomer
    result = int(round(copies))
    logger.info(f"N protein copy number: {genome_length} nt / {nt_per_protomer} nt per protomer "
                f"× {genome_copies} genome(s) = {result}")
    return result


def calculate_m_copy_number(
    virion_diameter_nm: float = 209.0,
    coverage_fraction: float = 0.24,
    dimer_footprint_nm2: float = 35.8,
    morphology: str = "spherical",
    virion_length_nm: Optional[float] = None,
) -> int:
    """
    Estimate M protein copy number from virion surface area and M-layer geometry.

    Formula (shape-independent):
        A_total    = calculate_virion_surface_area(morphology, diameter, length)
        A_M        = A_total × coverage_fraction
        N_dimers   = A_M / dimer_footprint_nm2
        N_monomers = N_dimers × 2

    Shape-specific surface area:
        Spherical   : A = 4π r²
        Filamentous : A = 2π r L   (capsule formula)

    Parameters
    ----------
    virion_diameter_nm : float
        Sphere diameter for spherical; tube cross-section diameter for
        filamentous (nm).
    coverage_fraction : float
        Fraction of inner membrane covered by M layer
        (default 0.24, from RSV cryo-ET, Kiss et al. 2014).
    dimer_footprint_nm2 : float
        Membrane-facing convex-hull area of the M dimer in nm²
        (default 35.8 nm², YZ-plane projection of PDB 4LP7, Leyrat et al. 2014).
    morphology : str
        "spherical" (default) or "filamentous".
    virion_length_nm : float | None
        Total tip-to-tip length for filamentous particles (nm).
        Required when morphology="filamentous".

    Returns
    -------
    int : Estimated M protein monomer copy number.
    """
    surface_area = calculate_virion_surface_area(morphology, virion_diameter_nm, virion_length_nm)
    m_covered_area = surface_area * coverage_fraction
    n_dimers = m_covered_area / dimer_footprint_nm2
    n_monomers = int(round(n_dimers * 2))
    logger.info(
        f"M protein copy number ({morphology}): surface={surface_area:.0f} nm², "
        f"M-covered={m_covered_area:.0f} nm², "
        f"dimers={n_dimers:.1f} → monomers={n_monomers}"
    )
    return n_monomers


def calculate_glycoprotein_copy_numbers(
    virion_diameter_nm: float = 209.0,
    spike_width_nm: float = 11.5,
    spike_spacing_nm: float = 8.0,
    ratio: Dict[str, int] = None,
    oligomeric_states: Dict[str, int] = None,
    morphology: str = "spherical",
    virion_length_nm: Optional[float] = None,
) -> Dict[str, int]:
    """
    Estimate F, G, and SH copy numbers using a hexagonal spike-packing model.

    Formula (shape-independent):
        A_virion       = calculate_virion_surface_area(morphology, diameter, length)
        d              = spike_width + spike_spacing      (centre-to-centre distance)
        A_spike        = (√3 / 2) × d²                   (Voronoi cell, hex packing)
        N_spikes       = A_virion / A_spike
        copies_i       = (ratio_i / Σratio) × N_spikes × oligomeric_state_i

    Shape-specific surface area:
        Spherical   : A = 4π r²
        Filamentous : A = 2π r L   (capsule formula)

    Parameters
    ----------
    virion_diameter_nm : float
        Sphere diameter for spherical; tube cross-section diameter for
        filamentous (nm).
    spike_width_nm : float
        Average glycoprotein spike width in nm (default 11.5, Walsh et al. 2015).
    spike_spacing_nm : float
        Inter-spike gap in nm (default 8.0 = average spacing).
    ratio : dict
        F:G:SH stoichiometric ratio (default {'F': 5, 'G': 3, 'SH': 1}).
    oligomeric_states : dict
        Subunits per spike: F=3 (trimer), G=1 (monomer), SH=5 (pentamer).
    morphology : str
        "spherical" (default) or "filamentous".
    virion_length_nm : float | None
        Total tip-to-tip length for filamentous particles (nm).
        Required when morphology="filamentous".

    Returns
    -------
    dict : {'F': int, 'G': int, 'SH': int} copy numbers per virion.

    References
    ----------
    Walsh et al. (2015)     — RSV spike width/spacing by EM.
    Conley et al. (2022)    — RSV cryo-ET hexagonal spike arrays.
    McLellan (2013), Battles (2017) — F homotrimer.
    Thammawat (2008)        — G monomer.
    Gan et al. (2012)       — RSV SH pentamer (applied by analogy to HMPV SH).
    Schildgen et al. (2011) — transcription gradient supporting F > G > SH.
    """
    if ratio is None:
        ratio = {'F': 5, 'G': 3, 'SH': 1}
    if oligomeric_states is None:
        oligomeric_states = {'F': 3, 'G': 1, 'SH': 5}

    surface_area = calculate_virion_surface_area(morphology, virion_diameter_nm, virion_length_nm)

    # Hexagonal packing
    d = spike_width_nm + spike_spacing_nm
    spike_area = (math.sqrt(3) / 2.0) * d ** 2

    n_spikes = surface_area / spike_area
    ratio_total = sum(ratio.values())

    copy_numbers = {}
    for protein, r in ratio.items():
        spikes_for_protein = (r / ratio_total) * n_spikes
        oligo = oligomeric_states.get(protein, 1)
        copy_numbers[protein] = int(round(spikes_for_protein * oligo))

    logger.info(
        f"Glycoprotein copy numbers ({morphology}, d={d:.1f} nm, {n_spikes:.0f} total spikes): "
        + ", ".join(f"{p}={c}" for p, c in copy_numbers.items())
    )
    return copy_numbers


def get_copy_numbers(
    genome_length: int,
    default_copy_numbers: Dict[str, int],
    use_calculated: bool = False,
    virion_diameter_nm: float = 209.0,
    nt_per_protomer: int = 7,
    genome_copies: int = 1,
    coverage_fraction: float = 0.24,
    dimer_footprint_nm2: float = 35.8,
    spike_width_nm: float = 11.5,
    spike_spacing_mode: str = 'average',
    spike_spacing_min_nm: float = 6.0,
    spike_spacing_avg_nm: float = 8.0,
    spike_spacing_max_nm: float = 10.0,
    ratio: Dict[str, int] = None,
    oligomeric_states: Dict[str, int] = None,
    morphology: str = "spherical",
    virion_length_nm: Optional[float] = None,
) -> Dict[str, int]:
    """
    Return protein copy numbers, either default or biophysically calculated.

    For proteins not covered by the calculation (P, L, M2-1, M2-2), the
    default values from default_copy_numbers are always used.

    Parameters
    ----------
    genome_length : int
        Viral genome length in nucleotides (needed for N protein calculation).
    default_copy_numbers : dict
        Fallback copy numbers (typically HMPV_COPY_NUMBERS from config).
    use_calculated : bool
        If True, calculate N, M, F, G, SH from biophysical models.
        If False (default), return default_copy_numbers unchanged.
    (remaining parameters are forwarded to the individual calc functions)

    Returns
    -------
    dict : Copy numbers for all 9 HMPV proteins.
    """
    if not use_calculated:
        logger.info("Using default copy numbers from config (use_calculated=False).")
        return dict(default_copy_numbers)

    logger.info("Calculating copy numbers from biophysical models (use_calculated=True).")

    # Start from defaults so P, L, M2-1, M2-2 are always present
    copy_numbers = dict(default_copy_numbers)

    logger.info(f"Morphology: {morphology}" + (f"  (length={virion_length_nm} nm)" if virion_length_nm else ""))

    # N protein  — genome-stoichiometry based; not surface-area dependent
    copy_numbers['N'] = calculate_n_copy_number(genome_length, nt_per_protomer, genome_copies)

    # M protein  — surface-area dependent
    copy_numbers['M'] = calculate_m_copy_number(
        virion_diameter_nm=virion_diameter_nm,
        coverage_fraction=coverage_fraction,
        dimer_footprint_nm2=dimer_footprint_nm2,
        morphology=morphology,
        virion_length_nm=virion_length_nm,
    )

    # Select spike spacing based on mode
    spacing_map = {
        'min': spike_spacing_min_nm,
        'average': spike_spacing_avg_nm,
        'max': spike_spacing_max_nm,
    }
    spike_spacing_nm = spacing_map.get(spike_spacing_mode, spike_spacing_avg_nm)
    logger.info(f"Spike spacing mode: '{spike_spacing_mode}' → {spike_spacing_nm} nm")

    # F, G, SH proteins  — surface-area dependent
    glyco = calculate_glycoprotein_copy_numbers(
        virion_diameter_nm=virion_diameter_nm,
        spike_width_nm=spike_width_nm,
        spike_spacing_nm=spike_spacing_nm,
        ratio=ratio,
        oligomeric_states=oligomeric_states,
        morphology=morphology,
        virion_length_nm=virion_length_nm,
    )
    copy_numbers.update(glyco)

    logger.info("Final copy numbers:")
    for protein, count in copy_numbers.items():
        logger.info(f"  {protein}: {count}")

    return copy_numbers
