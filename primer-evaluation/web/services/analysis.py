"""Analysis service wrapping the primer evaluation core library."""

import sys
from pathlib import Path
from typing import Any, Dict, Optional

# Add src/ to path for the installed package
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from primer_eval.validator import (
    Primer3Validator,
    PrimerAnalysis,
    PrimerPairAnalysis,
    SpecificityResult,
)
from web.config import config


class AnalysisService:
    def __init__(self):
        self.validator = Primer3Validator()
        self.max_mismatches = config.MAX_MISMATCHES
        self.allow_3prime_mismatches = config.ALLOW_3PRIME_MISMATCHES

    def analyze(
        self,
        forward: str,
        reverse: str,
        template: Optional[str] = None,
        max_mismatches: Optional[int] = None,
        allow_3prime_mismatches: Optional[int] = None,
    ) -> Dict[str, Any]:
        forward = forward.upper().strip()
        reverse = reverse.upper().strip()
        template = template.upper().strip() if template else None

        mm = max_mismatches if max_mismatches is not None else self.max_mismatches
        tp = allow_3prime_mismatches if allow_3prime_mismatches is not None else self.allow_3prime_mismatches

        if template and len(template) > config.MAX_TEMPLATE_LENGTH:
            raise ValueError(
                f"Template too long: {len(template)} bp. "
                f"Maximum: {config.MAX_TEMPLATE_LENGTH} bp"
            )

        if template:
            result = self.validator.analyze_primer_pair_with_template(
                forward, reverse, template, mm, tp
            )
        else:
            result = self.validator.analyze_primer_pair(forward, reverse)

        return result.to_dict()

    def generate_text_report(self, result: Dict[str, Any]) -> str:
        forward = PrimerAnalysis(**result["forward"])
        reverse = PrimerAnalysis(**result["reverse"])
        pair = result["pair"]
        specificity = result.get("specificity")

        if specificity:
            spec_obj = SpecificityResult(
                forward_matches=[(m["start"], m["end"], m["sequence"], m["mismatches"])
                                 for m in specificity["forward_matches"]],
                reverse_matches=[(m["start"], m["end"], m["sequence"], m["mismatches"])
                                 for m in specificity["reverse_matches"]],
                forward_specific=specificity["forward_specific"],
                reverse_specific=specificity["reverse_specific"],
                potential_products=[(p["fw_start"], p["rv_start"], p["length"])
                                     for p in specificity["potential_products"]],
                specificity_warnings=specificity["specificity_warnings"],
            )
        else:
            spec_obj = None

        analysis = PrimerPairAnalysis(
            forward=forward,
            reverse=reverse,
            tm_difference=pair["tm_difference"],
            heterodimer_dg=pair["heterodimer_dg"],
            three_prime_risk=pair["three_prime_risk"],
            three_prime_details=pair["three_prime_details"],
            warnings=result["warnings"],
            passed=result["passed"],
            specificity=spec_obj,
        )
        return self.validator.generate_report(analysis)


_analysis_service: Optional[AnalysisService] = None


def get_analysis_service() -> AnalysisService:
    global _analysis_service
    if _analysis_service is None:
        _analysis_service = AnalysisService()
    return _analysis_service
