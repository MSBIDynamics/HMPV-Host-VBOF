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
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, field

from .exceptions import MissingProteinDataError, MissingCopyNumberError

logger = logging.getLogger(__name__)


# =============================================================================
# AMINO ACID PROPERTIES
# =============================================================================
#
# Standard amino acid molecular weights and BiGG metabolite IDs
#
# SOURCES AND REFERENCES:
#
# 1. MOLECULAR WEIGHTS:
#    - IUPAC-IUB Joint Commission on Biochemical Nomenclature (1984)
#      "Nomenclature and Symbolism for Amino Acids"
#      Pure Appl Chem 56(5):595-624. DOI: 10.1351/pac198456050595
#      https://iupac.org/what-we-do/nomenclature/
#
# 2. BIGG DATABASE IDs:
#    - King et al. (2016) "BiGG Models: A platform for integrating genome-scale models"
#      Nucleic Acids Res 44(D1):D515-22. DOI: 10.1093/nar/gkv1049
#      http://bigg.ucsd.edu/universal/metabolites
#
# 3. AMINO ACID BIOCHEMISTRY:
#    - Berg, Tymoczko, Stryer (2015) "Biochemistry" 8th Edition
#      Chapter 2: Protein Composition and Structure
#      https://www.ncbi.nlm.nih.gov/books/NBK22364/
#
AMINO_ACID_INFO = {
    'A': {'name': 'Alanine', 'abbrev': 'Ala', 'mw': 89.09, 'bigg_id': 'ala__L_c'},
    'R': {'name': 'Arginine', 'abbrev': 'Arg', 'mw': 174.20, 'bigg_id': 'arg__L_c'},
    'N': {'name': 'Asparagine', 'abbrev': 'Asn', 'mw': 132.12, 'bigg_id': 'asn__L_c'},
    'D': {'name': 'Aspartate', 'abbrev': 'Asp', 'mw': 133.10, 'bigg_id': 'asp__L_c'},
    'C': {'name': 'Cysteine', 'abbrev': 'Cys', 'mw': 121.16, 'bigg_id': 'cys__L_c'},
    'E': {'name': 'Glutamate', 'abbrev': 'Glu', 'mw': 147.13, 'bigg_id': 'glu__L_c'},
    'Q': {'name': 'Glutamine', 'abbrev': 'Gln', 'mw': 146.15, 'bigg_id': 'gln__L_c'},
    'G': {'name': 'Glycine', 'abbrev': 'Gly', 'mw': 75.07, 'bigg_id': 'gly_c'},
    'H': {'name': 'Histidine', 'abbrev': 'His', 'mw': 155.16, 'bigg_id': 'his__L_c'},
    'I': {'name': 'Isoleucine', 'abbrev': 'Ile', 'mw': 131.17, 'bigg_id': 'ile__L_c'},
    'L': {'name': 'Leucine', 'abbrev': 'Leu', 'mw': 131.17, 'bigg_id': 'leu__L_c'},
    'K': {'name': 'Lysine', 'abbrev': 'Lys', 'mw': 146.19, 'bigg_id': 'lys__L_c'},
    'M': {'name': 'Methionine', 'abbrev': 'Met', 'mw': 149.21, 'bigg_id': 'met__L_c'},
    'F': {'name': 'Phenylalanine', 'abbrev': 'Phe', 'mw': 165.19, 'bigg_id': 'phe__L_c'},
    'P': {'name': 'Proline', 'abbrev': 'Pro', 'mw': 115.13, 'bigg_id': 'pro__L_c'},
    'S': {'name': 'Serine', 'abbrev': 'Ser', 'mw': 105.09, 'bigg_id': 'ser__L_c'},
    'T': {'name': 'Threonine', 'abbrev': 'Thr', 'mw': 119.12, 'bigg_id': 'thr__L_c'},
    'W': {'name': 'Tryptophan', 'abbrev': 'Trp', 'mw': 204.23, 'bigg_id': 'trp__L_c'},
    'Y': {'name': 'Tyrosine', 'abbrev': 'Tyr', 'mw': 181.19, 'bigg_id': 'tyr__L_c'},
    'V': {'name': 'Valine', 'abbrev': 'Val', 'mw': 117.15, 'bigg_id': 'val__L_c'},
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

