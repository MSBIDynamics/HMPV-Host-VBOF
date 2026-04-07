"""
Genome Analyzer Module
======================

This module provides functions to parse and analyze the HMPV genome sequence,
extracting nucleotide composition for VBOF construction.

Key Functions:
- load_genome: Load genome from FASTA file
- count_nucleotides: Count A, U, G, C in genome
- parse_gff_annotations: Parse gene annotations from GFF file
- calculate_genome_stoichiometry: Calculate nucleotide requirements for VBOF

Author: Syed Mushahid Hussain
"""

import logging
from pathlib import Path
from typing import Dict, Tuple, List, Optional
from dataclasses import dataclass

from .exceptions import MissingGenomeError

logger = logging.getLogger(__name__)
@dataclass
class GenomeInfo:
    """Container for genome information."""
    accession: str
    description: str
    sequence: str
    length: int
    nucleotide_counts: Dict[str, int]
    gc_content: float


@dataclass
class GeneAnnotation:
    """Container for gene annotation from GFF file."""
    gene_name: str
    start: int
    end: int
    strand: str
    product: str
    protein_id: str
    locus_tag: str




def load_genome(fasta_path: str) -> GenomeInfo:
    """
    Load and parse HMPV genome from FASTA file.
    
    Parameters:
    -----------
    fasta_path : str
        Path to the genome FASTA file (.fna)
    
    Returns:
    --------
    GenomeInfo : Dataclass containing genome information
    
    Raises:
    -------
    MissingGenomeError : If file not found or invalid format
    
    Example:
    --------
    >>> genome = load_genome("Data/genomic/GCF_002815375.1_ASM281537v1_genomic.fna")
    >>> print(f"Genome length: {genome.length} nt")
    """
    filepath = Path(fasta_path)
    
    if not filepath.exists():
        raise MissingGenomeError(f"Genome file not found: {fasta_path}", filepath=str(filepath))
    
    logger.info(f"Loading genome from: {fasta_path}")
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        raise MissingGenomeError(f"Failed to read genome file: {e}", filepath=str(filepath))
    
    # Parse FASTA format
    if not lines or not lines[0].startswith('>'):
        raise MissingGenomeError("Invalid FASTA format: missing header line", filepath=str(filepath))
    
    # Parse header
    header = lines[0].strip()[1:]  # Remove '>'
    parts = header.split(' ', 1)
    accession = parts[0]
    description = parts[1] if len(parts) > 1 else ""
    
    # Parse sequence (join all non-header lines)
    sequence = ''.join(line.strip() for line in lines[1:] if not line.startswith('>'))
    sequence = sequence.upper()
    
    if not sequence:
        raise MissingGenomeError("Empty genome sequence", filepath=str(filepath))
    
    # Count nucleotides
    nucleotide_counts = count_nucleotides(sequence)
    
    # Calculate GC content
    gc_content = (nucleotide_counts.get('G', 0) + nucleotide_counts.get('C', 0)) / len(sequence) * 100
    
    genome_info = GenomeInfo(
        accession=accession,
        description=description,
        sequence=sequence,
        length=len(sequence),
        nucleotide_counts=nucleotide_counts,
        gc_content=gc_content
    )
    
    logger.info(f"Loaded genome: {accession}")
    logger.info(f"  Length: {genome_info.length} nt")
    logger.info(f"  GC content: {genome_info.gc_content:.2f}%")
    logger.info(f"  Nucleotide counts: A={nucleotide_counts['A']}, T/U={nucleotide_counts['T']}, "
                f"G={nucleotide_counts['G']}, C={nucleotide_counts['C']}")
    
    return genome_info


def count_nucleotides(sequence: str) -> Dict[str, int]:
    """
    Count nucleotides in a RNA sequence.
    
    Parameters:
    -----------
    sequence : str
        Nucleotide sequence
    
    Returns:
    --------
    dict : Counts for each nucleotide {A, T/U, G, C}
    
    Note:
    -----
    For HMPV (negative-sense RNA virus), the genome is stored as DNA in databases.
    T counts will be converted to U for VBOF calculations.
    """
    sequence = sequence.upper()
    
    counts = {
        'A': sequence.count('A'),
        'T': sequence.count('T'),  # Will be U in RNA
        'G': sequence.count('G'),
        'C': sequence.count('C'),
        'U': sequence.count('U'),  # In case sequence is already RNA
    }
    
    # Combine T and U (for sequences that might have either)
    counts['T'] = counts['T'] + counts['U']
    
    return counts
def parse_gff_annotations(gff_path: str) -> List[GeneAnnotation]:
    """
    Parse gene annotations from GFF file.
    
    Parameters:
    -----------
    gff_path : str
        Path to the GFF annotation file
    
    Returns:
    --------
    list : List of GeneAnnotation objects
    
    Raises:
    -------
    MissingGenomeError : If file not found or invalid format
    """
    filepath = Path(gff_path)
    
    if not filepath.exists():
        raise MissingGenomeError(f"GFF file not found: {gff_path}", filepath=str(filepath))
    
    logger.info(f"Parsing GFF annotations from: {gff_path}")
    
    annotations = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip comments and empty lines
            if line.startswith('#') or not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            
            # Only process CDS features (protein-coding)
            feature_type = parts[2]
            if feature_type != 'CDS':
                continue
            
            # Parse fields
            seqid = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
            
            gene_name = attr_dict.get('gene', 'Unknown')
            product = attr_dict.get('product', 'Unknown')
            protein_id = attr_dict.get('protein_id', attr_dict.get('Name', 'Unknown'))
            locus_tag = attr_dict.get('locus_tag', 'Unknown')
            
            annotation = GeneAnnotation(
                gene_name=gene_name,
                start=start,
                end=end,
                strand=strand,
                product=product,
                protein_id=protein_id,
                locus_tag=locus_tag
            )
            
            annotations.append(annotation)
            logger.debug(f"  Found gene: {gene_name} ({product}) at {start}-{end}")
    
    logger.info(f"Parsed {len(annotations)} gene annotations")
    
    return annotations


def calculate_genome_stoichiometry(
    genome_info: GenomeInfo,
    copies_per_virion: int = 1
) -> Dict[str, float]:
    """
    Calculate nucleotide stoichiometry for VBOF.
    
    For HMPV (negative-sense ssRNA virus):
    - The virion contains one copy of the (-)ssRNA genome
    - During replication, both (+) and (-) strands are synthesized
    - We calculate requirements for the packaged genome
    
    Parameters:
    -----------
    genome_info : GenomeInfo
        Parsed genome information
    copies_per_virion : int
        Number of genome copies per virion (default: 1)
    
    Returns:
    --------
    dict : Stoichiometric coefficients for nucleotides
           Negative values indicate consumption
    
    Note:
    -----
    Uses BiGG metabolite IDs for compatibility with Recon models.
    T is converted to U for RNA.
    """
    nt_counts = genome_info.nucleotide_counts
    
    # Convert DNA counts to RNA (T -> U)
    # For (-)ssRNA genome: A in DNA = A in RNA, T in DNA = U in RNA
    rna_counts = {
        'A': nt_counts['A'],
        'U': nt_counts['T'],  # T becomes U in RNA
        'G': nt_counts['G'],
        'C': nt_counts['C'],
    }
    
    # BiGG metabolite IDs for NTPs (cytosolic)
    # Negative coefficients = consumption
    stoichiometry = {
        'atp_c': -rna_counts['A'] * copies_per_virion,
        'utp_c': -rna_counts['U'] * copies_per_virion,
        'gtp_c': -rna_counts['G'] * copies_per_virion,
        'ctp_c': -rna_counts['C'] * copies_per_virion,
    }
    
    # RNA polymerization releases pyrophosphate (PPi)
    # NTP -> NMP (in RNA) + PPi
    total_nt = sum(rna_counts.values()) * copies_per_virion
    stoichiometry['ppi_c'] = total_nt  # Positive = production
    
    logger.info(f"Genome stoichiometry calculated:")
    logger.info(f"  ATP: {stoichiometry['atp_c']}")
    logger.info(f"  UTP: {stoichiometry['utp_c']}")
    logger.info(f"  GTP: {stoichiometry['gtp_c']}")
    logger.info(f"  CTP: {stoichiometry['ctp_c']}")
    logger.info(f"  PPi produced: {stoichiometry['ppi_c']}")
    
    return stoichiometry


def get_genome_summary(genome_info: GenomeInfo) -> Dict:
    """
    Get a summary dictionary of genome information.
    
    Parameters:
    -----------
    genome_info : GenomeInfo
        Parsed genome information
    
    Returns:
    --------
    dict : Summary statistics
    """
    return {
        'accession': genome_info.accession,
        'length': genome_info.length,
        'gc_content': round(genome_info.gc_content, 2),
        'nucleotide_counts': genome_info.nucleotide_counts,
        'a_percent': round(genome_info.nucleotide_counts['A'] / genome_info.length * 100, 2),
        't_u_percent': round(genome_info.nucleotide_counts['T'] / genome_info.length * 100, 2),
        'g_percent': round(genome_info.nucleotide_counts['G'] / genome_info.length * 100, 2),
        'c_percent': round(genome_info.nucleotide_counts['C'] / genome_info.length * 100, 2),
    }

