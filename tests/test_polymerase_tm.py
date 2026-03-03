"""Unit tests for polymerase-tm package."""

import os
import tempfile

import pytest

from polymerase_tm import (
    tm,
    ta,
    reverse_complement,
    batch_tm,
    optimal_binding_length,
    check_pair,
    pcr_protocol,
    gc_content,
    primer_dimer,
    gibson_overlaps,
    restriction_scan,
    primer_quality,
    list_polymerases,
    from_csv,
    to_csv,
    RESTRICTION_ENZYMES,
)


# =====================================================================
# Core Tm/Ta — verified against NEB Tm Calculator
# =====================================================================

class TestTm:
    """Tm values verified against https://tmcalculator.neb.com/"""

    def test_q5_basic(self):
        assert tm("ATGTCCCTGCTCTTCTCTCGATGCAA") == 72

    def test_q5_short(self):
        assert tm("ATCGATCGATCG") == 47

    def test_taq_basic(self):
        assert tm("ATGTCCCTGCTCTTCTCTCGATGCAA", polymerase="taq") == 63

    def test_phusion(self):
        assert tm("ATGTCCCTGCTCTTCTCTCGATGCAA", polymerase="phusion_hf") == 74

    def test_gc_rich(self):
        t = tm("GCGCGCGCGCGCGCGCGCGC")
        assert t > 80  # very high GC → high Tm

    def test_at_rich(self):
        t = tm("AATTAATTAATTAATTAATT")
        assert t < 50  # very AT-rich → low Tm


class TestTa:
    def test_basic_pair(self):
        result_ta, t1, t2 = ta(
            "ATGTCCCTGCTCTTCTCTCGATGCAA",
            "GTGCCTCCGAGCCAGCACC",
        )
        assert result_ta == 72
        assert t1 == 72
        assert t2 == 75

    def test_dmso_correction(self):
        ta_no_dmso, _, _ = ta("ATGTCCCTGCTCTTCTCTCGATGCAA", "GTGCCTCCGAGCCAGCACC")
        ta_with_dmso, _, _ = ta(
            "ATGTCCCTGCTCTTCTCTCGATGCAA", "GTGCCTCCGAGCCAGCACC", dmso_pct=3
        )
        assert ta_with_dmso <= ta_no_dmso

    def test_taq_pair(self):
        result_ta, _, _ = ta(
            "ATGTCCCTGCTCTTCTCTCGATGCAA",
            "GTGCCTCCGAGCCAGCACC",
            polymerase="taq",
        )
        assert result_ta <= 68  # Taq cap is 68


# =====================================================================
# Automation functions
# =====================================================================

class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self):
        assert reverse_complement("GAATTC") == "GAATTC"  # EcoRI is palindromic

    def test_poly_a(self):
        assert reverse_complement("AAAA") == "TTTT"


class TestBatchTm:
    def test_multiple(self):
        seqs = ["ATCGATCGATCG", "GCGCGCGCGCGC", "ATGTCCCTGCTCTTCTCTCGATGCAA"]
        results = batch_tm(seqs)
        assert len(results) == 3
        for r in results:
            assert "tm" in r
            assert "gc_pct" in r
            assert "sequence" in r


class TestOptimalBindingLength:
    def test_finds_length(self):
        result = optimal_binding_length(
            "ATGTCCCTGCTCTTCTCTCGATGCAA", target_tm=60
        )
        assert result["tm"] >= 55  # should be near target
        assert result["length"] <= 26

    def test_returns_full_if_needed(self):
        result = optimal_binding_length("ATCGATCG", target_tm=80)
        # Can't reach 80 with 8 nt, should return full sequence
        assert result["length"] == 8


class TestCheckPair:
    def test_compatible_pair(self):
        result = check_pair(
            "ATGTCCCTGCTCTTCTCTCGATGCAA",
            "GTGCCTCCGAGCCAGCACC",
        )
        assert "fwd_tm" in result
        assert "rev_tm" in result
        assert "ta" in result
        assert "compatible" in result
        assert "additive" in result

    def test_additive_recommendation_gc_rich(self):
        result = check_pair("GCGCGCGCGCGCGCGCGCGC", "GCGCGCGCGCGCGCGC")
        assert result["additive"]["recommended"] is True
        assert "GC Enhancer" in result["additive"]["additive"]

    def test_additive_taq_dmso(self):
        result = check_pair(
            "GCGCGCGCGCGCGCGCGCGC",
            "GCGCGCGCGCGCGCGC",
            polymerase="taq",
        )
        assert result["additive"]["additive"] == "DMSO"


class TestPcrProtocol:
    def test_basic(self):
        result = pcr_protocol(
            "ATGTCCCTGCTCTTCTCTCGATGCAA",
            "GTGCCTCCGAGCCAGCACC",
            amplicon_length=1000,
        )
        assert "cycling" in result
        assert len(result["cycling"]) == 6  # init, denat, anneal, ext, final, hold
        assert result["total_time_min"] > 0
        assert "additive" in result

    def test_step_names(self):
        result = pcr_protocol(
            "ATGTCCCTGCTCTTCTCTCGATGCAA",
            "GTGCCTCCGAGCCAGCACC",
        )
        steps = [s["step"] for s in result["cycling"]]
        assert "Initial Denaturation" in steps
        assert "Hold" in steps


# =====================================================================
# New Analysis Functions
# =====================================================================

class TestPrimerDimer:
    def test_complementary_pair(self):
        # These primers have 3' complementarity
        result = primer_dimer("ATCGATCGATCG", "CGATCGATCGAT")
        assert "max_score" in result
        assert "risk_level" in result
        assert result["risk_level"] in ("none", "low", "moderate", "high")

    def test_no_dimer(self):
        result = primer_dimer("AAAAAAAAAA", "CCCCCCCCCC")
        assert result["risk_level"] == "none"

    def test_self_dimer(self):
        result = primer_dimer("ATCGATCGATCG", "TTTTTTTTTTTT")
        assert "fwd_self_dimer" in result


class TestGibsonOverlaps:
    def test_basic_design(self):
        result = gibson_overlaps(
            fwd_bind="ATGTCCCTGCTCTTCTCTCGATGCAA",
            rev_bind="GTGCCTCCGAGCCAGCACC",
            left_seq="AATTCCGGTTAACCGGATCC" * 3,
            right_seq="GCGATCGATCGATCGATCGAT" * 3,
            overlap_len=20,
        )
        assert len(result["fwd_full"]) == 20 + 26  # overlap + binding
        assert len(result["rev_full"]) == 20 + 19
        assert result["fwd_bind_tm"] == 72
        assert result["overlap_len"] == 20

    def test_custom_overlap(self):
        result = gibson_overlaps(
            fwd_bind="ATGTCCCTGCTCTTCTCTCGATGCAA",
            rev_bind="GTGCCTCCGAGCCAGCACC",
            left_seq="AATTCCGGTTAACCGGATCC" * 3,
            right_seq="GCGATCGATCGATCGATCGAT" * 3,
            overlap_len=25,
        )
        assert result["overlap_len"] == 25


class TestRestrictionScan:
    def test_ecori(self):
        hits = restriction_scan("ATGAATTCGATCG")
        ecori_hits = [h for h in hits if h["enzyme"] == "EcoRI"]
        assert len(ecori_hits) >= 1
        assert ecori_hits[0]["position"] == 2

    def test_no_sites(self):
        hits = restriction_scan("AAAAAAAAAAA")
        assert len(hits) == 0

    def test_multiple_sites(self):
        # Contains both EcoRI (GAATTC) and BamHI (GGATCC)
        hits = restriction_scan("GAATTCATCGGGATCC")
        enzymes = set(h["enzyme"] for h in hits)
        assert "EcoRI" in enzymes
        assert "BamHI" in enzymes

    def test_custom_enzymes(self):
        custom = {"TestEnzyme": "AAAA"}
        hits = restriction_scan("GGGAAAAGGGG", enzymes=custom)
        assert len(hits) >= 1
        assert hits[0]["enzyme"] == "TestEnzyme"


class TestPrimerQuality:
    def test_good_primer(self):
        result = primer_quality("ATGTCCCTGCTCTTCTCTCGATGCAA")
        assert result["score"] >= 60
        assert result["grade"] in ("A", "B", "C")
        assert result["length_ok"] is True

    def test_too_short(self):
        result = primer_quality("ATCGATC")
        assert result["score"] < 80
        assert any("short" in i.lower() for i in result["issues"])

    def test_gc_rich(self):
        result = primer_quality("GCGCGCGCGCGCGCGCGCGCGCGC")
        assert any("gc" in i.lower() for i in result["issues"])

    def test_homopolymer_run(self):
        result = primer_quality("ATCGAAAAATCGATCGATCG")
        assert result["max_run"] >= 5
        assert any("homopolymer" in i.lower() for i in result["issues"])

    def test_perfect_primer(self):
        # A well-designed 20mer with 50% GC and GC clamp
        result = primer_quality("ATCGATCGGCTAGCATCGAC")
        assert result["score"] >= 70  # Near ideal, minor deductions possible


# =====================================================================
# Utility
# =====================================================================

class TestGcContent:
    def test_all_gc(self):
        assert gc_content("GCGCGCGC") == 1.0

    def test_all_at(self):
        assert gc_content("ATATATAT") == 0.0

    def test_half(self):
        assert gc_content("ATGC") == 0.5


class TestListPolymerases:
    def test_returns_list(self):
        polys = list_polymerases()
        assert len(polys) >= 22
        assert all("key" in p for p in polys)
        assert all("description" in p for p in polys)


class TestCsvRoundTrip:
    def test_write_and_read(self):
        # Create a temp CSV with primer data
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False, newline=""
        ) as f:
            f.write("name,fwd,rev\n")
            f.write("pair1,ATGTCCCTGCTCTTCTCTCGATGCAA,GTGCCTCCGAGCCAGCACC\n")
            f.write("pair2,ATCGATCGATCGATCGATCG,GCTAGCTAGCTAGCTAGCT\n")
            tmp_input = f.name

        try:
            results = from_csv(tmp_input)
            assert len(results) == 2
            assert "fwd_tm" in results[0]
            assert "ta" in results[0]

            # Write results
            tmp_output = tmp_input + "_out.csv"
            to_csv(results, tmp_output)
            assert os.path.exists(tmp_output)
            os.unlink(tmp_output)
        finally:
            os.unlink(tmp_input)
