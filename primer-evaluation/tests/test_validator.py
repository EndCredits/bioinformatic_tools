"""Tests for primer evaluation core library."""

import pytest
from primer_eval.validator import (
    Primer3Validator,
    PrimerAnalysis,
    PrimerPairAnalysis,
    SpecificityResult,
    SequenceMatcher,
)


class TestPrimer3Validator:
    def test_basic_analysis(self):
        validator = Primer3Validator()
        result = validator.analyze_primer_pair(
            "ATGCCCTGAGCTAAAGCTG",
            "GTGAGCTTTGTCTCGGTGA",
        )
        assert isinstance(result, PrimerPairAnalysis)
        assert result.forward.length == 19
        assert result.reverse.length == 19
        assert 40 <= result.forward.gc_content <= 60
        assert 40 <= result.reverse.gc_content <= 60

    def test_warnings_on_gc_out_of_range(self):
        validator = Primer3Validator()
        # Very GC-rich primers
        result = validator.analyze_primer_pair(
            "GCGCGCGCGCGCGCGCGC",
            "GCGCGCGCGCGCGCGCGC",
        )
        assert len(result.warnings) > 0
        assert any("GC%" in w for w in result.warnings)

    def test_validation_too_long_sequence(self):
        validator = Primer3Validator()
        with pytest.raises(ValueError, match="must be < 60 bp"):
            validator.analyze_primer_pair(
                "A" * 60,
                "A" * 60,
            )

    def test_validation_too_short_sequence(self):
        validator = Primer3Validator()
        with pytest.raises(ValueError, match="should be >= 15 bp"):
            validator.analyze_primer_pair(
                "ATGC",
                "GCAT",
            )

    def test_invalid_bases(self):
        validator = Primer3Validator()
        with pytest.raises(ValueError, match="invalid bases"):
            validator.analyze_primer_pair(
                "ATGCTGCATGCX",
                "GCATGCATGCA",
            )

    def test_case_insensitive(self):
        validator = Primer3Validator()
        result1 = validator.analyze_primer_pair("ATGCATGC", "GCATGCAT")
        result2 = validator.analyze_primer_pair("atgcatgc", "gcatgcat")
        assert result1.forward.tm == result2.forward.tm
        assert result1.forward.gc_content == result2.forward.gc_content

    def test_specificity_analysis(self):
        validator = Primer3Validator()
        template = "ATGCCCTGAGCTAAAGCTGTGATGGCCCTGGCTGTCCTCCTGCTACTCTGCCTGCTGCTCAAGCTCTGGGGCACCGGCTTCAGCTGGATGTGCTGCGAGGCCTATGAGCAGGCCACCCAGCAGCTGAGCGAGAAGCTGCAGAGGGCCGAGGACGCCGAGCTGCCCGAGGACGAGCTGGACGAGCTGGACGAGGAGCTGGACGAGGAGCTGGACGAGGAGCTGGACGAGGAGCTGGACGAGGAGCTGGACGAGGAGCTGGACGAGGAGCTGGACGAGGAG"
        result = validator.analyze_primer_pair_with_template(
            "ATGCCCTGAGCTAAAGCTG",
            "GTGAGCTTTGTCTCGGTGA",
            template,
        )
        assert result.specificity is not None
        assert isinstance(result.specificity, SpecificityResult)


class TestSequenceMatcher:
    def test_exact_match(self):
        matcher = SequenceMatcher(max_mismatches=0)
        template = "ATCGATCG"
        matches = matcher.find_matches("ATCG", template, is_reverse=False)
        assert len(matches) == 1
        assert matches[0] == (0, 4, "ATCG", 0)

    def test_mismatch_count(self):
        matcher = SequenceMatcher(max_mismatches=1)
        template = "ATTGATCG"
        matches = matcher.find_matches("ATCG", template, is_reverse=False)
        assert len(matches) == 1
        assert matches[0][3] == 1  # 1 mismatch

    def test_reverse_complement(self):
        seq = "ATCG"
        rc = SequenceMatcher.reverse_complement(seq)
        assert rc == "CGAT"

    def test_no_match(self):
        matcher = SequenceMatcher(max_mismatches=2)
        template = "AAAAAAAA"
        matches = matcher.find_matches("GCGCGCGC", template, is_reverse=False)
        assert len(matches) == 0

    def test_potential_products(self):
        fwd = [(0, 10, "AAAAAAAAAA", 0)]
        rev = [(100, 110, "TTTTTTTTTT", 0)]
        products = SequenceMatcher.find_potential_products(fwd, rev)
        assert len(products) == 1
        assert products[0] == (0, 100, 101)
