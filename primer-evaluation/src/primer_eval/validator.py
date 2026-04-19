#!/usr/bin/env python3
"""Core primer validation module using primer3-py."""

import logging
import os
import re
from dataclasses import dataclass, asdict
from typing import List, Tuple, Optional

from primer3 import (
    calc_tm, calc_hairpin, calc_homodimer,
    calc_heterodimer, calc_end_stability
)

logger = logging.getLogger(__name__)

# 1 kcal = 1000 cal
CAL_TO_KCAL = 1000.0


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class PrimerAnalysis:
    """Single primer analysis result."""
    sequence: str
    length: int
    gc_content: float
    tm: float
    hairpin_dg: float  # kcal/mol
    homodimer_dg: float  # kcal/mol

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class SpecificityResult:
    """Primer-template specificity check result."""
    forward_matches: List[Tuple[int, int, str, int]]  # (start, end, seq, mismatches)
    reverse_matches: List[Tuple[int, int, str, int]]
    forward_specific: bool
    reverse_specific: bool
    potential_products: List[Tuple[int, int, int]]  # (fw_start, rv_start, length)
    specificity_warnings: List[str]


@dataclass
class PrimerPairAnalysis:
    """Primer pair analysis result."""
    forward: PrimerAnalysis
    reverse: PrimerAnalysis
    tm_difference: float
    heterodimer_dg: float  # kcal/mol
    three_prime_risk: bool
    three_prime_details: List[str]
    warnings: List[str]
    passed: bool
    specificity: Optional[SpecificityResult] = None

    def to_dict(self) -> dict:
        result = {
            "forward": self.forward.to_dict(),
            "reverse": self.reverse.to_dict(),
            "pair": {
                "tm_difference": self.tm_difference,
                "heterodimer_dg": self.heterodimer_dg,
                "three_prime_risk": self.three_prime_risk,
                "three_prime_details": self.three_prime_details,
            },
            "warnings": self.warnings,
            "passed": self.passed,
        }
        if self.specificity:
            result["specificity"] = {
                "forward_matches": [
                    {"start": s, "end": e, "sequence": seq, "mismatches": m}
                    for s, e, seq, m in self.specificity.forward_matches
                ],
                "reverse_matches": [
                    {"start": s, "end": e, "sequence": seq, "mismatches": m}
                    for s, e, seq, m in self.specificity.reverse_matches
                ],
                "forward_specific": self.specificity.forward_specific,
                "reverse_specific": self.specificity.reverse_specific,
                "potential_products": [
                    {"fw_start": fw, "rv_start": rv, "length": length}
                    for fw, rv, length in self.specificity.potential_products
                ],
                "specificity_warnings": self.specificity.specificity_warnings,
            }
        return result


# ---------------------------------------------------------------------------
# SequenceMatcher
# ---------------------------------------------------------------------------

class SequenceMatcher:
    """Sequence matcher for primer-template specificity checking."""

    def __init__(self, max_mismatches: int = 3, allow_3prime_mismatches: int = 1):
        self.max_mismatches = max_mismatches
        self.allow_3prime_mismatches = allow_3prime_mismatches

    @staticmethod
    def reverse_complement(seq: str) -> str:
        """Return reverse complement of a DNA sequence."""
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        return "".join(complement[b] for b in reversed(seq.upper()))

    def count_mismatches(self, primer: str, template: str, start: int) -> Tuple[int, int]:
        """Count total and 3'-end mismatches."""
        primer = primer.upper()
        template = template.upper()
        end_pos = start + len(primer)

        if end_pos > len(template) or len(template[start:end_pos]) != len(primer):
            return float("inf"), float("inf")

        total = sum(1 for i in range(len(primer)) if primer[i] != template[start + i])

        end_len = min(3, len(primer))
        end_mismatches = sum(
            1 for i in range(len(primer) - end_len, len(primer))
            if primer[i] != template[start + i]
        )
        return total, end_mismatches

    def find_matches(self, primer: str, template: str, is_reverse: bool = False) -> List[Tuple[int, int, str, int]]:
        """Find all binding sites for a primer in the template."""
        matches = []
        primer_len = len(primer)
        template_len = len(template)

        if primer_len > template_len:
            return matches

        search_seq = self.reverse_complement(primer) if is_reverse else primer

        for i in range(template_len - primer_len + 1):
            mismatches, end_mismatches = self.count_mismatches(search_seq, template, i)
            if mismatches <= self.max_mismatches and end_mismatches <= self.allow_3prime_mismatches:
                matches.append((i, i + primer_len, template[i : i + primer_len], mismatches))

        return matches

    @staticmethod
    def find_potential_products(
        forward_matches: List[Tuple[int, int, str, int]],
        reverse_matches: List[Tuple[int, int, str, int]],
    ) -> List[Tuple[int, int, int]]:
        """Find potential amplification products from primer matches."""
        products = []
        for fw_start, _, _, _ in forward_matches:
            for rv_start, _, _, _ in reverse_matches:
                if fw_start < rv_start:
                    products.append((fw_start, rv_start, rv_start - fw_start + 1))
        return products


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------

def read_fasta_file(filepath: str) -> str:
    """Read FASTA file and return concatenated sequence."""
    with open(filepath, "r") as f:
        lines = f.readlines()
    seq = ""
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            continue
        seq += line
    return seq.replace(" ", "").replace("\n", "").replace("\r", "").upper()


def read_genbank_file(filepath: str) -> str:
    """Read GenBank file and return sequence from ORIGIN section."""
    with open(filepath, "r") as f:
        content = f.read()

    origin_start = content.find("ORIGIN")
    if origin_start == -1:
        raise ValueError("Invalid GenBank file: ORIGIN not found")

    sequence_start = origin_start + len("ORIGIN")
    end_marker = content.find("//", sequence_start)
    if end_marker == -1:
        raise ValueError("Invalid GenBank file: // end marker not found")

    sequence_part = content[sequence_start:end_marker]
    seq = ""
    for line in sequence_part.split("\n"):
        seq_line = re.sub(r"^\s*\d+\s+", "", line)
        seq_line = re.sub(r"\s+", "", seq_line)
        seq += seq_line
    return seq.upper()


def read_template_file(filepath: str) -> str:
    """Read template from FASTA or GenBank file."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Template file not found: {filepath}")

    ext = os.path.splitext(filepath)[1].lower()
    if ext in (".fasta", ".fa"):
        return read_fasta_file(filepath)
    elif ext in (".gb", ".genbank"):
        return read_genbank_file(filepath)
    else:
        try:
            return read_fasta_file(filepath)
        except Exception:
            try:
                return read_genbank_file(filepath)
            except Exception as e:
                raise ValueError(
                    f"Unsupported file format: {ext}. Supported: .fasta, .fa, .gb, .genbank"
                ) from e


# ---------------------------------------------------------------------------
# Primer3Validator
# ---------------------------------------------------------------------------

class Primer3Validator:
    """Primer evaluation tool based on primer3-py."""

    def __init__(self, **kwargs):
        self.params = {
            "mv_conc": 50.0,
            "dv_conc": 1.5,
            "dntp_conc": 0.6,
            "dna_conc": 50.0,
            "temp_c": 37.0,
            "max_loop": 30,
            "output_structure": True,
            "max_nn_length": 60,
            "tm_method": "santalucia",
            "salt_corrections_method": "santalucia",
        }

        self.thresholds = {
            "min_gc": 40.0,
            "max_gc": 60.0,
            "max_tm_diff": 5.0,
            "max_hairpin_dg": -3.0,
            "max_homodimer_dg": -5.0,
            "max_heterodimer_dg": -5.0,
            "max_three_prime_dg": -3.0,
            # specificity thresholds
            "max_binding_sites": 5,
            "max_potential_products": 3,
        }

        for key, value in kwargs.items():
            if key in self.params:
                self.params[key] = value
            elif key in self.thresholds:
                self.thresholds[key] = value

    def _validate_sequence(self, seq: str, name: str = "Primer") -> str:
        if not seq:
            raise ValueError(f"{name} sequence cannot be empty")
        seq = seq.upper().strip()
        valid_bases = set("ATGCNRYKMSWBDHV")
        if not all(base in valid_bases for base in seq):
            invalid = set(seq) - valid_bases
            raise ValueError(f"{name} contains invalid bases: {invalid}")
        if len(seq) >= 60:
            raise ValueError(
                f"{name} must be < 60 bp for thermodynamic calculations (got {len(seq)} bp)"
            )
        if len(seq) < 15:
            raise ValueError(
                f"{name} should be >= 15 bp for PCR (got {len(seq)} bp)"
            )
        return seq

    def _validate_template(self, seq: str) -> str:
        if not seq:
            raise ValueError("Template sequence cannot be empty")
        seq = seq.upper().strip()
        valid_bases = set("ATGCNRYKMSWBDHV")
        if not all(base in valid_bases for base in seq):
            invalid = set(seq) - valid_bases
            raise ValueError(f"Template contains invalid bases: {invalid}")
        return seq

    def _calculate_gc_content(self, seq: str) -> float:
        gc = sum(1 for b in seq if b in "GC")
        gc += sum(0.5 for b in seq if b == "S")  # S = G or C
        gc += sum(0.5 for b in seq if b == "N")  # N = any, 50% GC probability
        return (gc / len(seq)) * 100 if seq else 0

    @staticmethod
    def _convert_cal_to_kcal(value: float) -> float:
        return value / CAL_TO_KCAL

    def _safe_thermo(self, func, *args, **kwargs) -> Tuple[float, bool]:
        try:
            result = func(*args, **kwargs)
            if hasattr(result, "structure_found") and result.structure_found:
                if hasattr(result, "dg"):
                    return self._convert_cal_to_kcal(result.dg), True
            return 0.0, False
        except RuntimeError as e:
            if "sequence too long" in str(e).lower():
                raise ValueError(f"Sequence too long for thermodynamic analysis: {args[0]}") from e
            return 0.0, False
        except Exception:
            return 0.0, False

    def _analyze_single_primer(self, seq: str) -> PrimerAnalysis:
        tm_params = {
            k: self.params[k] for k in [
                "mv_conc", "dv_conc", "dntp_conc", "dna_conc",
                "max_nn_length", "tm_method", "salt_corrections_method",
            ]
        }
        tm = calc_tm(seq, **tm_params)

        thermo_params = {
            k: self.params[k] for k in [
                "mv_conc", "dv_conc", "dntp_conc", "dna_conc",
                "temp_c", "max_loop", "output_structure",
            ]
        }
        hairpin_dg, _ = self._safe_thermo(calc_hairpin, seq, **thermo_params)
        homodimer_dg, _ = self._safe_thermo(calc_homodimer, seq, **thermo_params)

        return PrimerAnalysis(
            sequence=seq,
            length=len(seq),
            gc_content=round(self._calculate_gc_content(seq), 2),
            tm=round(tm, 2),
            hairpin_dg=round(hairpin_dg, 2),
            homodimer_dg=round(homodimer_dg, 2),
        )

    def _check_three_prime_stability(self, seq1: str, seq2: str) -> Tuple[bool, float, str]:
        thermo_params = {
            k: self.params[k] for k in [
                "mv_conc", "dv_conc", "dntp_conc", "dna_conc", "temp_c", "max_loop",
            ]
        }
        try:
            result = calc_end_stability(seq1, seq2, **thermo_params)
            if hasattr(result, "structure_found") and result.structure_found and hasattr(result, "dg"):
                dg_kcal = self._convert_cal_to_kcal(result.dg)
                is_risky = dg_kcal < self.thresholds["max_three_prime_dg"]
                detail = f"3'-end stability \u0394G = {dg_kcal:.2f} kcal/mol"
                return is_risky, dg_kcal, detail
        except Exception as e:
            logger.warning("3'-end stability check failed: %s", e)

        return False, 0.0, ""

    def analyze_primer_pair(self, forward_seq: str, reverse_seq: str) -> PrimerPairAnalysis:
        fwd = self._validate_sequence(forward_seq, "Forward primer")
        rev = self._validate_sequence(reverse_seq, "Reverse primer")
        warnings = []

        fwd_analysis = self._analyze_single_primer(fwd)
        rev_analysis = self._analyze_single_primer(rev)

        thermo_params = {
            k: self.params[k] for k in [
                "mv_conc", "dv_conc", "dntp_conc", "dna_conc",
                "temp_c", "max_loop", "output_structure",
            ]
        }
        heterodimer_dg, _ = self._safe_thermo(calc_heterodimer, fwd, rev, **thermo_params)

        fwd_risk, fwd_dg, fwd_detail = self._check_three_prime_stability(fwd, rev)
        rev_risk, rev_dg, rev_detail = self._check_three_prime_stability(rev, fwd)

        three_prime_risk = fwd_risk or rev_risk
        three_prime_details = []
        if fwd_risk:
            three_prime_details.append(f"Forward 3'-end: {fwd_detail}")
        if rev_risk:
            three_prime_details.append(f"Reverse 3'-end: {rev_detail}")

        tm_diff = abs(fwd_analysis.tm - rev_analysis.tm)

        if not (self.thresholds["min_gc"] <= fwd_analysis.gc_content <= self.thresholds["max_gc"]):
            warnings.append(
                f"[GC%] Forward: {fwd_analysis.gc_content:.1f}% "
                f"(expected {self.thresholds['min_gc']}-{self.thresholds['max_gc']}%)"
            )
        if not (self.thresholds["min_gc"] <= rev_analysis.gc_content <= self.thresholds["max_gc"]):
            warnings.append(
                f"[GC%] Reverse: {rev_analysis.gc_content:.1f}% "
                f"(expected {self.thresholds['min_gc']}-{self.thresholds['max_gc']}%)"
            )
        if tm_diff > self.thresholds["max_tm_diff"]:
            warnings.append(
                f"[Tm] Tm difference: {tm_diff:.1f}\u00b0C "
                f"(max {self.thresholds['max_tm_diff']}\u00b0C)"
            )
        if fwd_analysis.hairpin_dg < self.thresholds["max_hairpin_dg"]:
            warnings.append(
                f"[Hairpin] Forward: \u0394G = {fwd_analysis.hairpin_dg:.2f} kcal/mol "
                f"(threshold {self.thresholds['max_hairpin_dg']} kcal/mol)"
            )
        if rev_analysis.hairpin_dg < self.thresholds["max_hairpin_dg"]:
            warnings.append(
                f"[Hairpin] Reverse: \u0394G = {rev_analysis.hairpin_dg:.2f} kcal/mol "
                f"(threshold {self.thresholds['max_hairpin_dg']} kcal/mol)"
            )
        if fwd_analysis.homodimer_dg < self.thresholds["max_homodimer_dg"]:
            warnings.append(
                f"[Homodimer] Forward: \u0394G = {fwd_analysis.homodimer_dg:.2f} kcal/mol "
                f"(threshold {self.thresholds['max_homodimer_dg']} kcal/mol)"
            )
        if rev_analysis.homodimer_dg < self.thresholds["max_homodimer_dg"]:
            warnings.append(
                f"[Homodimer] Reverse: \u0394G = {rev_analysis.homodimer_dg:.2f} kcal/mol "
                f"(threshold {self.thresholds['max_homodimer_dg']} kcal/mol)"
            )
        if heterodimer_dg < self.thresholds["max_heterodimer_dg"]:
            warnings.append(
                f"[Heterodimer] Pair: \u0394G = {heterodimer_dg:.2f} kcal/mol "
                f"(threshold {self.thresholds['max_heterodimer_dg']} kcal/mol)"
            )
        if three_prime_risk:
            warnings.extend(f"[3'-end risk] {d}" for d in three_prime_details)

        return PrimerPairAnalysis(
            forward=fwd_analysis,
            reverse=rev_analysis,
            tm_difference=round(tm_diff, 2),
            heterodimer_dg=round(heterodimer_dg, 2),
            three_prime_risk=three_prime_risk,
            three_prime_details=three_prime_details,
            warnings=warnings,
            passed=len(warnings) == 0,
        )

    def analyze_specificity(
        self,
        forward_seq: str,
        reverse_seq: str,
        template_seq: str,
        max_mismatches: int = 3,
        allow_3prime_mismatches: int = 1,
    ) -> SpecificityResult:
        fwd = self._validate_sequence(forward_seq, "Forward primer")
        rev = self._validate_sequence(reverse_seq, "Reverse primer")
        template = self._validate_template(template_seq)

        matcher = SequenceMatcher(
            max_mismatches=max_mismatches,
            allow_3prime_mismatches=allow_3prime_mismatches,
        )
        fwd_matches = matcher.find_matches(fwd, template, is_reverse=False)
        rev_matches = matcher.find_matches(rev, template, is_reverse=True)

        max_sites = self.thresholds["max_binding_sites"]
        forward_specific = len(fwd_matches) <= max_sites
        reverse_specific = len(rev_matches) <= max_sites

        potential_products = matcher.find_potential_products(fwd_matches, rev_matches)

        specificity_warnings = []
        if len(fwd_matches) > max_sites:
            specificity_warnings.append(
                f"[Specificity] Forward has {len(fwd_matches)} binding sites (expected <={max_sites})"
            )
        if len(rev_matches) > max_sites:
            specificity_warnings.append(
                f"[Specificity] Reverse has {len(rev_matches)} binding sites (expected <={max_sites})"
            )
        if len(potential_products) > self.thresholds["max_potential_products"]:
            specificity_warnings.append(
                f"[Specificity] {len(potential_products)} potential products "
                f"(expected <={self.thresholds['max_potential_products']}, non-specific amplification possible)"
            )
        if potential_products:
            lengths = [p[2] for p in potential_products]
            if any(l < 100 or l > 10000 for l in lengths):
                specificity_warnings.append(
                    f"[Specificity] Unusual product lengths: {lengths}"
                )

        return SpecificityResult(
            forward_matches=fwd_matches,
            reverse_matches=rev_matches,
            forward_specific=forward_specific,
            reverse_specific=reverse_specific,
            potential_products=potential_products,
            specificity_warnings=specificity_warnings,
        )

    def analyze_primer_pair_with_template(
        self,
        forward_seq: str,
        reverse_seq: str,
        template_seq: str,
        max_mismatches: int = 3,
        allow_3prime_mismatches: int = 1,
    ) -> PrimerPairAnalysis:
        pair_analysis = self.analyze_primer_pair(forward_seq, reverse_seq)
        spec = self.analyze_specificity(
            forward_seq, reverse_seq, template_seq,
            max_mismatches, allow_3prime_mismatches,
        )
        pair_analysis.specificity = spec
        pair_analysis.warnings.extend(spec.specificity_warnings)
        pair_analysis.passed = len(pair_analysis.warnings) == 0
        return pair_analysis

    def generate_report(self, analysis: PrimerPairAnalysis) -> str:
        from datetime import datetime

        lines = []
        lines.append("=" * 60)
        lines.append("PRIMER PAIR QUALITY ASSESSMENT REPORT")
        lines.append("=" * 60)
        lines.append(f"Analysis time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("-" * 60)

        for name, primer in [("FORWARD PRIMER", analysis.forward), ("REVERSE PRIMER", analysis.reverse)]:
            lines.append(f"\n{name}")
            lines.append("-" * 30)
            lines.append(f"Sequence:  5'-{primer.sequence}-3'")
            lines.append(f"Length:    {primer.length} bp")
            lines.append(f"GC%:       {primer.gc_content:.2f}%")
            lines.append(f"Tm:        {primer.tm:.2f} \u00b0C")
            lines.append(f"Hairpin \u0394G:  {primer.hairpin_dg:.2f} kcal/mol")
            lines.append(f"Homodimer \u0394G: {primer.homodimer_dg:.2f} kcal/mol")

        lines.append("\nPRIMER PAIR")
        lines.append("-" * 30)
        lines.append(f"Tm difference:    {analysis.tm_difference:.2f} \u00b0C")
        lines.append(f"Heterodimer \u0394G: {analysis.heterodimer_dg:.2f} kcal/mol")

        if analysis.three_prime_risk:
            lines.append("3'-end stability risks:")
            lines.extend(f"  \u2022 {d}" for d in analysis.three_prime_details)
        else:
            lines.append("3'-end stability: OK")

        if analysis.specificity:
            lines.append("\nSPECIFICITY ANALYSIS")
            lines.append("-" * 30)
            lines.append(f"Forward binding sites: {len(analysis.specificity.forward_matches)}")
            lines.append(f"Reverse binding sites: {len(analysis.specificity.reverse_matches)}")
            lines.append(f"Forward specific: {'Yes' if analysis.specificity.forward_specific else 'No'}")
            lines.append(f"Reverse specific: {'Yes' if analysis.specificity.reverse_specific else 'No'}")
            lines.append(f"Potential products: {len(analysis.specificity.potential_products)}")
            if analysis.specificity.potential_products:
                lines.append("Product lengths:")
                for i, (_, _, length) in enumerate(analysis.specificity.potential_products[:5]):
                    lines.append(f"  Product {i+1}: {length} bp")
                if len(analysis.specificity.potential_products) > 5:
                    lines.append(f"  ... and {len(analysis.specificity.potential_products) - 5} more")

        lines.append("\nASSESSMENT RESULT")
        lines.append("-" * 30)
        if analysis.passed:
            lines.append("\u2705 PASSED - No significant issues detected")
        else:
            lines.append("\u274c FAILED - Issues detected:")
            for i, w in enumerate(analysis.warnings, 1):
                lines.append(f"  {i}. {w}")

        lines.append("=" * 60)
        return "\n".join(lines)
