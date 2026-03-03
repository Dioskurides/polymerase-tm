# -*- coding: utf-8 -*-
"""Tests for the mutagenesis module (NEB Base Changer integration)."""

import pytest
from polymerase_tm import (
    BaseChanger, SDMPrimer, MutagenesisResult,
    select_codon, parse_aa_mutation,
    calc_sdm_tm, owczarzy_bivariate, calc_nn_raw,
    GENETIC_CODES, get_codon_table, get_aa_to_codons,
    CODON_TABLE, ECOLI_CODON_USAGE,
)


# =====================================================================
# Test Template — a realistic ORF for testing
# pUC19 lacZ-alpha fragment (first 240 nt, simplified)
# =====================================================================
TEMPLATE = (
    "ATGACCATGATTACGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTG"
    "GCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGA"
    "AGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCCTG"
    "ATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATGGTGCACTCTC"
)  # 240 nt


class TestCodonSelection:
    """Tests for codon selection logic."""

    def test_usage_mode_ecoli(self):
        """E. coli usage mode should pick highest-frequency codon."""
        # Leucine: CTG has highest usage (52.0) in E. coli
        codon = select_codon("L", mode="usage")
        assert codon == "CTG"

    def test_parsimony_mode(self):
        """Parsimony mode should pick codon with fewest base changes."""
        # AAA (Lys) -> AGA (Arg) = 1 change; CGT (Arg) = 3 changes
        codon = select_codon("R", original_codon="AAA", mode="parsimony")
        assert codon == "AGA"

    def test_usage_mode_alanine(self):
        """Most used Ala codon in E. coli should be GCG."""
        codon = select_codon("A", mode="usage")
        assert codon == "GCG"

    def test_unknown_aa_raises(self):
        with pytest.raises(ValueError, match="No codons found"):
            select_codon("X")

    def test_stop_codon(self):
        codon = select_codon("*", mode="usage")
        assert codon in ["TAA", "TAG", "TGA"]


class TestMutationParsing:
    """Tests for mutation format parsing."""

    def test_single_mutation(self):
        result = parse_aa_mutation("M1A")
        assert len(result) == 1
        assert result[0]["old_aa"] == "M"
        assert result[0]["position"] == 1
        assert result[0]["new_aa"] == "A"

    def test_multiple_mutations(self):
        result = parse_aa_mutation("M1A K2R F32H")
        assert len(result) == 3

    def test_comma_separated(self):
        result = parse_aa_mutation("M1A,K2R,F32H")
        assert len(result) == 3

    def test_explicit_codon(self):
        result = parse_aa_mutation("K2R:CGC")
        assert result[0]["explicit_codon"] == "CGC"

    def test_invalid_format_raises(self):
        with pytest.raises(ValueError, match="Invalid mutation format"):
            parse_aa_mutation("INVALID")

    def test_stop_codon_mutation(self):
        result = parse_aa_mutation("Q10*")
        assert result[0]["new_aa"] == "*"


class TestGeneticCodes:
    """Tests for genetic code tables."""

    def test_standard_code_is_default(self):
        table = get_codon_table(1)
        assert table["ATG"] == "M"
        assert table["TAA"] == "*"

    def test_vertebrate_mitochondrial(self):
        table = get_codon_table(2)
        assert table["AGA"] == "*"  # Stop in vertebrate mito
        assert table["TGA"] == "W"  # Trp instead of Stop

    def test_yeast_mitochondrial(self):
        table = get_codon_table(3)
        assert table["CTG"] == "T"  # Thr instead of Leu

    def test_all_12_codes_exist(self):
        expected_ids = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14]
        for code_id in expected_ids:
            table = get_codon_table(code_id)
            assert len(table) == 64

    def test_invalid_code_raises(self):
        with pytest.raises(ValueError, match="Unknown genetic code"):
            get_codon_table(99)

    def test_aa_to_codons_reverse(self):
        aa_codons = get_aa_to_codons(1)
        # Every codon for Met should be ATG
        assert aa_codons["M"] == ["ATG"]
        # Leucine has 6 codons
        assert len(aa_codons["L"]) == 6


class TestCalcSdmTm:
    """Tests for the bivariate Tm calculation."""

    def test_basic_20mer(self):
        """20-mer with mixed GC content should give reasonable Tm."""
        seq = "ATGACCATGATTACGAATTC"  # 20 nt, 40% GC
        tm = calc_sdm_tm(seq)
        assert 45 < tm < 70, f"Tm={tm} out of expected range"

    def test_gc_rich_primer(self):
        """GC-rich primer should have higher Tm."""
        gc_rich = "GCGCCCGCGCCGGCGCGCCG"
        at_rich = "AATATAAATAATTTAAATAA"
        tm_gc = calc_sdm_tm(gc_rich)
        tm_at = calc_sdm_tm(at_rich)
        assert tm_gc > tm_at

    def test_bivariate_vs_mono(self):
        """Bivariate correction with Mg2+ should differ from no Mg2+."""
        seq = "ATGACCATGATTACGAATTC"
        raw = calc_nn_raw(seq, 500)
        tm_biv = owczarzy_bivariate(raw, seq, 50, 2)
        tm_mono = owczarzy_bivariate(raw, seq, 50, 0)
        # They should differ because one has Mg2+
        assert tm_biv != tm_mono

    def test_dmso_correction(self):
        """DMSO should decrease Tm by 0.6 per percent."""
        seq = "ATGACCATGATTACGAATTC"
        tm_no_dmso = calc_sdm_tm(seq, dmso_pct=0)
        tm_5_dmso = calc_sdm_tm(seq, dmso_pct=5)
        assert abs((tm_no_dmso - tm_5_dmso) - 3.0) < 0.1


class TestBaseChanger:
    """Integration tests for the BaseChanger class."""

    def test_init(self):
        bc = BaseChanger(TEMPLATE)
        assert bc.template == TEMPLATE
        assert bc.min_tm == 55.0

    def test_empty_template_raises(self):
        with pytest.raises(ValueError, match="must not be empty"):
            BaseChanger("")

    def test_invalid_bases_raises(self):
        with pytest.raises(ValueError, match="Invalid characters"):
            BaseChanger("ATGXYZ")

    def test_point_mutation_basic(self):
        """K2R on the template: ATG AAA -> ATG CGT (or best codon)."""
        bc = BaseChanger(TEMPLATE)
        # Template starts ATG ACC (M T ...), codon 2 is ACC (Thr)
        result = bc.point_mutation("T2A")
        assert isinstance(result, MutagenesisResult)
        assert result.forward.length > 0
        assert result.reverse.length > 0
        assert result.forward.tm > 0
        assert result.reverse.tm > 0

    def test_point_mutation_preserves_template(self):
        """Original template should not be modified."""
        bc = BaseChanger(TEMPLATE)
        original = TEMPLATE
        bc.point_mutation("T2A")
        assert bc.template == original

    def test_point_mutation_codon_change(self):
        """Verify the codon is correctly changed in the mutated sequence."""
        bc = BaseChanger(TEMPLATE)
        result = bc.point_mutation("T2A")
        # Codon 2 starts at nt index 3
        new_codon = result.mutated_sequence[3:6]
        # Should encode Alanine
        assert CODON_TABLE.get(new_codon) == "A"

    def test_explicit_codon(self):
        bc = BaseChanger(TEMPLATE)
        result = bc.point_mutation("T2A:GCT")
        assert result.new_codon == "GCT"

    def test_ta_calculation(self):
        """Ta should be min(Tm_fwd, Tm_rev) + 1."""
        bc = BaseChanger(TEMPLATE)
        result = bc.point_mutation("T2A")
        assert result.ta == min(result.forward.tm, result.reverse.tm) + 1

    def test_deletion(self):
        bc = BaseChanger(TEMPLATE)
        result = bc.deletion(start=10, length=3)
        assert isinstance(result, MutagenesisResult)
        assert len(result.mutated_sequence) == len(TEMPLATE) - 3

    def test_insertion(self):
        bc = BaseChanger(TEMPLATE)
        result = bc.insertion(position=10, insert_seq="AAAAAA")
        assert isinstance(result, MutagenesisResult)
        assert len(result.mutated_sequence) == len(TEMPLATE) + 6

    def test_substitution(self):
        bc = BaseChanger(TEMPLATE)
        result = bc.substitution(start=10, replacement="GGG", length=3)
        assert isinstance(result, MutagenesisResult)
        assert len(result.mutated_sequence) == len(TEMPLATE)
        assert result.mutated_sequence[9:12] == "GGG"

    def test_batch(self):
        bc = BaseChanger(TEMPLATE)
        results = bc.batch("T2A A3G")
        assert isinstance(results, list)
        assert len(results) == 2

    def test_confine_to_tails(self):
        bc = BaseChanger(TEMPLATE, confine_to_tails=True)
        result = bc.point_mutation("T2A")
        assert isinstance(result, MutagenesisResult)

    def test_codon_mode_parsimony(self):
        bc = BaseChanger(TEMPLATE, codon_mode="parsimony")
        result = bc.point_mutation("T2A")
        assert isinstance(result, MutagenesisResult)

    def test_genetic_code_2(self):
        """Vertebrate mitochondrial code should work."""
        bc = BaseChanger(TEMPLATE, genetic_code=2)
        result = bc.point_mutation("T2A")
        assert isinstance(result, MutagenesisResult)

    def test_min_length_parameter(self):
        bc = BaseChanger(TEMPLATE, min_length=20)
        result = bc.point_mutation("T2A")
        # All primers should be >= 20 nt
        assert result.forward.length >= 20

    def test_sdm_primer_dataclass(self):
        bc = BaseChanger(TEMPLATE)
        result = bc.point_mutation("T2A")
        fwd = result.forward
        assert isinstance(fwd, SDMPrimer)
        assert fwd.direction == "FWD"
        assert 0 <= fwd.gc_pct <= 100
        assert len(fwd.sequence) > 0

    def test_summary_table(self):
        bc = BaseChanger(TEMPLATE)
        result = bc.point_mutation("T2A")
        table = result.summary_table()
        assert "FWD" in table
        assert "REV" in table
        assert "Mutation" in table
