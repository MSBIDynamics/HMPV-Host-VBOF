"""
Custom Exceptions for HMPV VBOF Construction
============================================

This module defines custom exceptions for error handling throughout the VBOF
construction pipeline. Following project rules: descriptive, clear messages.

Author: Syed Mushahid Hussain
"""

import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger(__name__)


class HMPVModelError(Exception):
    """Base exception for HMPV metabolic modeling errors."""
    pass


class MissingGenomeError(HMPVModelError):
    """
    Raised when genome data is missing or invalid.
    
    Examples:
        - Genome file not found
        - Invalid FASTA format
        - Empty sequence
    """
    def __init__(self, message: str, filepath: str = None):
        self.filepath = filepath
        super().__init__(f"Genome Error: {message}" + 
                        (f" (File: {filepath})" if filepath else ""))


class MissingProteinDataError(HMPVModelError):
    """
    Raised when protein data is incomplete or missing.
    
    Examples:
        - Protein file not found
        - Missing protein sequences
        - Invalid amino acid sequences
    """
    def __init__(self, message: str, protein_name: str = None):
        self.protein_name = protein_name
        super().__init__(f"Protein Data Error: {message}" +
                        (f" (Protein: {protein_name})" if protein_name else ""))


class MissingCopyNumberError(HMPVModelError):
    """
    Raised when protein copy number data is missing.
    
    Examples:
        - Copy number not provided for a protein
        - Copy number is zero or negative
    """
    def __init__(self, message: str, protein_name: str = None):
        self.protein_name = protein_name
        super().__init__(f"Copy Number Error: {message}" +
                        (f" (Protein: {protein_name})" if protein_name else ""))


class EnergyDataMissingError(HMPVModelError):
    """
    Raised when energy calculation data is missing.
    
    Examples:
        - Invalid energy coefficients
        - Missing ATP/GTP cost data
    """
    def __init__(self, message: str):
        super().__init__(f"Energy Data Error: {message}")


class IDMappingError(HMPVModelError):
    """
    Raised when metabolite ID cannot be mapped to model.
    
    Examples:
        - Metabolite ID not found in host model
        - Compartment mismatch
        - Invalid BiGG identifier
    """
    def __init__(self, message: str, metabolite_id: str = None, suggestions: list = None):
        self.metabolite_id = metabolite_id
        self.suggestions = suggestions or []
        error_msg = f"ID Mapping Error: {message}"
        if metabolite_id:
            error_msg += f" (Metabolite: {metabolite_id})"
        if suggestions:
            error_msg += f" Suggestions: {suggestions[:5]}"
        super().__init__(error_msg)


class ModelIntegrationError(HMPVModelError):
    """
    Raised when VBOF cannot be integrated into host model.
    
    Examples:
        - Model validation failed
        - Reaction already exists
        - Infeasible model after integration
    """
    def __init__(self, message: str):
        super().__init__(f"Model Integration Error: {message}")


class VBOFConstructionError(HMPVModelError):
    """
    Raised when VBOF construction fails.
    
    Examples:
        - Invalid stoichiometric coefficients
        - Missing required components
        - Mass balance violation
    """
    def __init__(self, message: str):
        super().__init__(f"VBOF Construction Error: {message}")

