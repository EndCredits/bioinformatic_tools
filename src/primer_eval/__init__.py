"""Primer Evaluation Tool - Evaluate DNA primer pairs for PCR experiments."""

from primer_eval.validator import Primer3Validator, PrimerAnalysis, PrimerPairAnalysis, SpecificityResult

__version__ = "1.0.0"
__all__ = ["Primer3Validator", "PrimerAnalysis", "PrimerPairAnalysis", "SpecificityResult"]
