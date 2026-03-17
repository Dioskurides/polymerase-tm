# -*- coding: utf-8 -*-
"""NEBaseChanger — Site-directed mutagenesis primer design.

Exact reproduction of the NEB Base Changer tool (v2.7.2):
    https://nebasechanger.neb.com/

Designs primers for site-directed mutagenesis using the
Q5 Site-Directed Mutagenesis Kit (NEB #E0554).

Algorithms
----------
- Tm: SantaLucia (1998) nearest-neighbor + Owczarzy (2008) bivariate salt.
- Primer extension: iterative 3' extension to minimum Tm, GC clamp.
- Mutation placement: ≤6 nt centred in FWD, >6 nt split to 5' tails.

References
----------
1. SantaLucia J Jr (1998) PNAS 95:1460-5.
2. Allawi HT, SantaLucia J Jr (1997-1999) Mismatch parameters.
3. Owczarzy R et al. (2008) Biochemistry 47:5336-53.
"""

from __future__ import annotations

import re
import math
from dataclasses import dataclass, field
from typing import Optional

import primer3
from .core import calc_nn_raw, calc_sdm_tm
from .batch import reverse_complement
from .constants import (
    CODON_TABLE, AA_TO_CODONS, ECOLI_CODON_USAGE,
    AA3_TO_AA1, AA1_TO_AA3,
    Q5_SDM_MONO_MM, Q5_SDM_DIVALENT_MM, Q5_SDM_PRIMER_NM,
    Q5_SDM_MIN_TM, Q5_SDM_MAX_TM, Q5_SDM_MIN_LEN, Q5_SDM_MIN_ANNEAL,
    get_codon_table, get_aa_to_codons,
)


# =====================================================================
# Data Classes
# =====================================================================

@dataclass
class SDMPrimer:
    """A single site-directed mutagenesis primer."""
    sequence: str
    length: int
    tm: float
    gc_pct: float
    direction: str      # "FWD" or "REV"
    bind_region: str     # template-binding portion only
    mutation_bases: str  # mutated bases (lowercase in full sequence)
    risk: Optional[str] = None  # secondary structure risk warning

    def __repr__(self) -> str:
        risk_str = f", RISK={self.risk}" if self.risk else ""
        return (
            f"SDMPrimer({self.direction}, {self.length} nt, "
            f"Tm={self.tm:.0f}°C, GC={self.gc_pct:.0f}%{risk_str})"
        )


@dataclass
class MutagenesisResult:
    """Result of a single mutagenesis primer design."""
    description: str
    forward: SDMPrimer
    reverse: SDMPrimer
    ta: float
    mutated_sequence: str
    original_codon: Optional[str] = None
    new_codon: Optional[str] = None
    warnings: list[str] = field(default_factory=list)

    def __repr__(self) -> str:
        return (
            f"MutagenesisResult({self.description}, "
            f"Ta={self.ta:.0f}°C)"
        )

    def summary_table(self) -> str:
        """Return a formatted summary string."""
        lines = [
            f"Mutation: {self.description}",
            f"{'Primer':<10} {'ID':<12} {'Oligo':<45} {'Len':>4} {'Tm':>4} {'Ta':>4}",
            "-" * 80,
        ]
        for p in [self.forward, self.reverse]:
            # Show mutation bases in lowercase
            oligo = p.sequence
            lines.append(
                f"{'':10} {p.direction + ':1':<12} {oligo:<45} "
                f"{p.length:>4d} {p.tm:>4.0f} {self.ta:>4.0f}"
            )
        if self.original_codon and self.new_codon:
            lines.append(
                f"Codon change: {self.original_codon} → {self.new_codon}"
            )
        if self.forward.risk:
            lines.append(f"FWD Risk: {self.forward.risk}")
        if self.reverse.risk:
            lines.append(f"REV Risk: {self.reverse.risk}")
        if self.warnings:
            lines.append("Warnings: " + "; ".join(self.warnings))
        return "\n".join(lines)


# =====================================================================
# Codon Selection Helpers
# =====================================================================

def _hamming(a: str, b: str) -> int:
    """Count base differences between two equal-length strings."""
    return sum(1 for x, y in zip(a.upper(), b.upper()) if x != y)


def select_codon(
    target_aa: str,
    original_codon: str = "",
    mode: str = "usage",
    genetic_code: int = 1,
) -> str:
    """Select the best codon for a target amino acid.

    Parameters
    ----------
    target_aa : str
        One-letter amino acid code (e.g., "A" for alanine).
    original_codon : str
        Original codon in the template (for parsimony comparison).
    mode : str
        "usage" — highest E. coli codon usage frequency.
        "parsimony" — fewest base changes vs. original_codon.
    genetic_code : int
        NCBI genetic code ID. Default 1 (Standard).

    Returns
    -------
    str
        Selected codon (uppercase).
    """
    target_aa = target_aa.upper()
    codons = get_aa_to_codons(genetic_code).get(target_aa, [])
    if not codons:
        raise ValueError(
            f"No codons found for amino acid '{target_aa}' "
            f"in genetic code {genetic_code}."
        )

    if mode == "parsimony" and original_codon:
        original_codon = original_codon.upper()
        # Find codon(s) with minimum Hamming distance
        min_dist = min(_hamming(c, original_codon) for c in codons)
        best = [c for c in codons if _hamming(c, original_codon) == min_dist]
        if len(best) == 1:
            return best[0]
        # Tie-break by usage frequency
        return max(best, key=lambda c: ECOLI_CODON_USAGE.get(c, 0))

    # Default: codon usage mode
    return max(codons, key=lambda c: ECOLI_CODON_USAGE.get(c, 0))


# =====================================================================
# Mutation Parsing
# =====================================================================

# Pattern: [old AA][position][new AA] or [old AA][position][new AA]:[codon]
_MUTATION_RE = re.compile(
    r"([A-Z*])(\d+)([A-Z*])(?::([ATCG]{3}))?",
    re.IGNORECASE,
)


def parse_aa_mutation(mutation_str: str) -> list[dict]:
    """Parse amino acid mutation notation.

    Formats: "M2A", "F32H", "K2R:CGC", "M1A F32H C45Y"
    Multiple mutations separated by spaces, commas, or semicolons.

    Returns
    -------
    list[dict]
        Each dict has keys: old_aa, position, new_aa, explicit_codon.
    """
    mutations = []
    # Split on whitespace, commas, semicolons
    parts = re.split(r"[,;\s]+", mutation_str.strip())
    for part in parts:
        if not part:
            continue
        m = _MUTATION_RE.match(part.upper())
        if not m:
            raise ValueError(
                f"Invalid mutation format: '{part}'. "
                "Expected format: [old AA][position][new AA] "
                "e.g., M2A, F32H, K2R:CGC"
            )
        mutations.append({
            "old_aa": m.group(1),
            "position": int(m.group(2)),
            "new_aa": m.group(3),
            "explicit_codon": m.group(4) if m.group(4) else None,
        })
    return mutations


# =====================================================================
# Primer Design Engine
# =====================================================================

def _gc_frac(seq: str) -> float:
    """GC fraction of a sequence."""
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def _build_primer(
    template: str,
    start: int,
    direction: str,
    min_tm: float = Q5_SDM_MIN_TM,
    max_tm: float = Q5_SDM_MAX_TM,
    min_len: int = Q5_SDM_MIN_LEN,
    mono_mM: float = Q5_SDM_MONO_MM,
    divalent_mM: float = Q5_SDM_DIVALENT_MM,
    primer_conc_nM: float = Q5_SDM_PRIMER_NM,
    method: Optional[str] = None,
) -> tuple[str, float]:
    """Extend a primer from `start` in the given direction until Tm target.

    Parameters
    ----------
    template : str
        Full template sequence (uppercase).
    start : int
        0-based start position of the annealing region.
    direction : str
        "forward" (extend 3') or "reverse" (extend 5', sequence is revcomp).
    min_tm : float
        Minimum Tm target.
    max_tm : float
        Maximum Tm (for GC clamp extension).
    min_len : int
        Minimum primer length.

    Returns
    -------
    (sequence, tm) : tuple[str, float]
    """
    template = template.upper()
    tlen = len(template)

    if direction == "forward":
        # Extend rightward from start
        best_seq = ""
        best_tm = 0.0
        for end in range(start + min_len, min(start + 60, tlen) + 1):
            seq = template[start:end]
            try:
                if method:
                    from .core import tm
                    t = tm(seq, salt_mM=int(mono_mM), method=method)
                else:
                    t = calc_sdm_tm(seq, mono_mM, divalent_mM, primer_conc_nM)
            except ValueError:
                continue
            if t >= min_tm:
                best_seq = seq
                best_tm = t
                # Try GC clamp
                if seq[-1] in "GC":
                    break
                # Keep extending for GC clamp, up to max_tm
                if t > max_tm:
                    break
        if not best_seq:
            # Could not reach min_tm — return longest possible
            seq = template[start:min(start + 60, tlen)]
            try:
                if method:
                    from .core import tm
                    t = tm(seq, salt_mM=int(mono_mM), method=method)
                else:
                    t = calc_sdm_tm(seq, mono_mM, divalent_mM, primer_conc_nM)
            except ValueError:
                t = 0.0
            return seq, t
        return best_seq, best_tm

    else:
        # Reverse primer — extend leftward from start
        best_seq = ""
        best_tm = 0.0
        for begin in range(start - min_len, max(start - 60, -1), -1):
            if begin < 0:
                break
            seq = template[begin:start]
            rc = reverse_complement(seq)
            try:
                if method:
                    from .core import tm
                    t = tm(rc, salt_mM=int(mono_mM), method=method)
                else:
                    t = calc_sdm_tm(rc, mono_mM, divalent_mM, primer_conc_nM)
            except ValueError:
                continue
            if t >= min_tm:
                best_seq = rc
                best_tm = t
                # GC clamp check (3' of revcomp = 5' of original)
                if rc[-1] in "GC":
                    break
                if t > max_tm:
                    break
        if not best_seq:
            seq = template[max(0, start - 60):start]
            rc = reverse_complement(seq)
            try:
                if method:
                    from .core import tm
                    t = tm(rc, salt_mM=int(mono_mM), method=method)
                else:
                    t = calc_sdm_tm(rc, mono_mM, divalent_mM, primer_conc_nM)
            except ValueError:
                t = 0.0
            return rc, t
        return best_seq, best_tm


def _balance_primers(
    fwd_seq: str, rev_seq: str,
    fwd_tm: float, rev_tm: float,
    template: str,
    fwd_start: int, rev_end: int,
    mono_mM: float = Q5_SDM_MONO_MM,
    divalent_mM: float = Q5_SDM_DIVALENT_MM,
    primer_conc_nM: float = Q5_SDM_PRIMER_NM,
    max_diff: float = 5.0,
    method: Optional[str] = None,
) -> tuple[str, float, str, float]:
    """Balance Tm between FWD and REV primers by extending the shorter one.

    Returns (fwd_seq, fwd_tm, rev_seq, rev_tm).
    """
    if abs(fwd_tm - rev_tm) <= max_diff:
        return fwd_seq, fwd_tm, rev_seq, rev_tm

    template = template.upper()

    # Extend the lower-Tm primer
    if fwd_tm < rev_tm:
        # Extend FWD 3'
        for ext in range(1, 15):
            end = fwd_start + len(fwd_seq) + ext
            if end > len(template):
                break
            seq = template[fwd_start:end]
            try:
                if method:
                    from .core import tm
                    t = tm(seq, salt_mM=int(mono_mM), method=method)
                else:
                    t = calc_sdm_tm(seq, mono_mM, divalent_mM, primer_conc_nM)
            except ValueError:
                continue
            if abs(t - rev_tm) <= max_diff or t >= rev_tm:
                return seq, t, rev_seq, rev_tm
    else:
        # Extend REV 5' (which means extending leftward in template)
        rc_orig = reverse_complement(rev_seq)
        # Find where rev binds in template
        rev_start = template.find(rc_orig)
        if rev_start >= 0:
            for ext in range(1, 15):
                begin = rev_start - ext
                if begin < 0:
                    break
                seq = template[begin:rev_start + len(rc_orig)]
                rc = reverse_complement(seq)
                try:
                    if method:
                        from .core import tm
                        t = tm(rc, salt_mM=int(mono_mM), method=method)
                    else:
                        t = calc_sdm_tm(rc, mono_mM, divalent_mM, primer_conc_nM)
                except ValueError:
                    continue
                if abs(t - fwd_tm) <= max_diff or t >= fwd_tm:
                    return fwd_seq, fwd_tm, rc, t

    return fwd_seq, fwd_tm, rev_seq, rev_tm


# =====================================================================
# BaseChanger — Main API
# =====================================================================

class BaseChanger:
    """NEB Base Changer — SDM primer designer for Q5 SDM Kit.

    Parameters
    ----------
    template : str
        Full plasmid/template DNA sequence (5'→3').
    orf_start : int
        1-based position of the ORF start (first base of ATG). Default 1.
    genetic_code : int
        NCBI genetic code ID. Default 1 (Standard).
    circular : bool
        If True (default), treat template as circular plasmid. Primers
        can wrap around the sequence junction. Matches NEB Base Changer
        default behavior.
    min_length : int
        Minimum primer length in nt. Default 15.
    min_tm : float
        Minimum primer Tm in °C. Default 55.
    confine_to_tails : bool
        If True, mutations are always placed as 5' non-annealing tails.
    preserve_degenerate : bool
        If True, IUPAC ambiguity codes are kept as-is (not expanded).
    merge_mutations : bool
        If True, nearby mutations (≤6 bases apart) are merged into one
        primer pair.
    codon_mode : str
        "usage" or "parsimony". Default "usage".
    dmso_pct : float
        DMSO percentage for Tm correction. Default 0.

    Examples
    --------
    >>> from polymerase_tm import BaseChanger
    >>> bc = BaseChanger("ATGAAAGCAATTTTCGTACTG...")
    >>> result = bc.point_mutation("K2R")
    >>> print(result.forward.sequence)
    >>> print(result.ta)
    """

    def __init__(
        self,
        template: str,
        orf_start: int = 1,
        genetic_code: int = 1,
        circular: bool = True,
        min_length: int = Q5_SDM_MIN_LEN,
        min_tm: float = Q5_SDM_MIN_TM,
        confine_to_tails: bool = False,
        preserve_degenerate: bool = False,
        merge_mutations: bool = False,
        codon_mode: str = "usage",
        dmso_pct: float = 0.0,
        method: Optional[str] = None,
    ):
        # Clean and validate template
        self.template = re.sub(r"\s+", "", template).upper()
        if not self.template:
            raise ValueError("Template sequence must not be empty.")
        invalid = set(self.template) - set("ATCG")
        if invalid and not preserve_degenerate:
            raise ValueError(
                f"Invalid characters in template: {sorted(invalid)}. "
                "Only A, T, G, C allowed (or enable preserve_degenerate)."
            )

        self.orf_start = orf_start  # 1-based
        self.genetic_code = genetic_code
        self.circular = circular
        self.min_length = min_length
        self.min_tm = min_tm
        self.max_tm = Q5_SDM_MAX_TM
        self.confine_to_tails = confine_to_tails
        self.preserve_degenerate = preserve_degenerate
        self.merge_mutations = merge_mutations
        self.codon_mode = codon_mode
        self.dmso_pct = dmso_pct
        self.method = method

        # Build codon table for this genetic code
        self._codon_table = get_codon_table(genetic_code)
        self._aa_to_codons = get_aa_to_codons(genetic_code)

    def _check_primer_risk(self, seq: str) -> Optional[str]:
        """Check for hairpins or homodimers using Primer3."""
        try:
            hp = primer3.calc_hairpin(seq)
            if hp.tm > 45:
                # Use .tm for result objects
                return f"Hairpin (Tm={hp.tm:.1f}°C)"
            
            dimer = primer3.calc_homodimer(seq)
            if dimer.tm > 40:
                return f"Homodimer (Tm={dimer.tm:.1f}°C)"
        except (ImportError, Exception):
            pass
        return None

    def _codon_at(self, aa_pos: int) -> tuple[str, int]:
        """Get the codon and its 0-based template index for a 1-based AA position.

        Returns (codon, template_idx).
        """
        # orf_start is 1-based
        nt_idx = (self.orf_start - 1) + (aa_pos - 1) * 3
        if nt_idx < 0 or nt_idx + 3 > len(self.template):
            raise ValueError(
                f"Codon position {aa_pos} (nt {nt_idx + 1}-{nt_idx + 3}) "
                f"is out of template range (1-{len(self.template)})."
            )
        codon = self.template[nt_idx:nt_idx + 3]
        return codon, nt_idx

    def point_mutation(self, mutations: str) -> MutagenesisResult | list[MutagenesisResult]:
        """Design primers for amino acid point mutation(s).

        Parameters
        ----------
        mutations : str
            Mutation(s) in format: "K2R", "M1A F32H", "K2R:CGC" (with
            explicit codon). Multiple mutations separated by spaces,
            commas, or semicolons.

        Returns
        -------
        MutagenesisResult or list[MutagenesisResult]
            Single result if one mutation, list if multiple.
        """
        parsed = parse_aa_mutation(mutations)
        if not parsed:
            raise ValueError("No valid mutations found in input.")

        results = []
        for mut in parsed:
            result = self._design_point_mutation(mut)
            results.append(result)

        if len(results) == 1:
            return results[0]
        return results

    def _design_point_mutation(self, mut: dict) -> MutagenesisResult:
        """Design primers for a single amino acid point mutation."""
        old_aa = mut["old_aa"]
        new_aa = mut["new_aa"]
        pos = mut["position"]
        explicit_codon = mut.get("explicit_codon")

        # Get original codon
        original_codon, nt_idx = self._codon_at(pos)

        # Verify old AA (warning, not error — matches NEB behavior)
        warnings = []
        expected_aa = self._codon_table.get(original_codon, "?")
        if expected_aa != old_aa:
            warnings.append(
                f"Expected {old_aa} at position {pos}, "
                f"but template has {expected_aa} ({original_codon})."
            )

        # Select new codon
        if explicit_codon:
            new_codon = explicit_codon.upper()
        else:
            new_codon = select_codon(
                new_aa, original_codon, self.codon_mode, self.genetic_code
            )

        # Build mutated template
        mutated = (
            self.template[:nt_idx]
            + new_codon
            + self.template[nt_idx + 3:]
        )

        # Determine the mutation region (nt changes)
        n_changes = _hamming(original_codon, new_codon)
        description = (
            f"Replace {original_codon.lower()} with "
            f"{new_codon.upper()} at codon {pos} ({old_aa}{pos}{new_aa})"
        )

        # Design primers
        fwd, rev = self._design_substitution_primers(
            mutated, nt_idx, 3, new_codon, method=self.method
        )

        # Ta is min(Tm_fwd, Tm_rev) + 1, capped at 72°C (matching NEB)
        ta = min(fwd.tm, rev.tm) + 1
        ta = min(ta, 72)

        # Check for secondary structure risks
        fwd.risk = self._check_primer_risk(fwd.sequence)
        rev.risk = self._check_primer_risk(rev.sequence)
        
        warnings = []
        if fwd.risk:
            warnings.append(f"Forward primer risk: {fwd.risk}")
        if rev.risk:
            warnings.append(f"Reverse primer risk: {rev.risk}")

        return MutagenesisResult(
            description=description,
            forward=fwd,
            reverse=rev,
            ta=ta,
            mutated_sequence=mutated,
            original_codon=original_codon,
            new_codon=new_codon,
            warnings=warnings,
        )



    def substitution(
        self,
        start: int,
        replacement: str,
        length: int = 0,
    ) -> MutagenesisResult:
        """Design primers for a nucleotide substitution.

        Parameters
        ----------
        start : int
            1-based start position in the template.
        replacement : str
            Replacement nucleotide sequence.
        length : int
            Number of template bases to replace. 0 = same as replacement length.

        Returns
        -------
        MutagenesisResult
        """
        replacement = replacement.upper()
        idx = start - 1
        if length == 0:
            length = len(replacement)

        if idx < 0 or idx + length > len(self.template):
            raise ValueError(
                f"Substitution region ({start}..{start + length - 1}) "
                f"exceeds template length ({len(self.template)})."
            )

        original = self.template[idx:idx + length]
        mutated = self.template[:idx] + replacement + self.template[idx + length:]

        description = f"Replace {original.lower()} with {replacement} at position {start}"

        fwd, rev = self._design_substitution_primers(
            mutated, idx, len(replacement), replacement, method=self.method
        )

        ta = min(fwd.tm, rev.tm) + 1
        ta = min(ta, 72)

        # Check for secondary structure risks
        warnings = []
        fwd_risk = self._check_primer_risk(fwd.sequence)
        rev_risk = self._check_primer_risk(rev.sequence)
        if fwd_risk:
            warnings.append(f"Forward primer risk: {fwd_risk}")
        if rev_risk:
            warnings.append(f"Reverse primer risk: {rev_risk}")

        return MutagenesisResult(
            description=description,
            forward=fwd,
            reverse=rev,
            ta=ta,
            mutated_sequence=mutated,
            warnings=warnings,
        )

    def deletion(self, start: int, length: int) -> MutagenesisResult:
        """Design primers for a nucleotide deletion.

        Parameters
        ----------
        start : int
            1-based start position of the deletion.
        length : int
            Number of bases to delete.
        """
        idx = start - 1
        if idx < 0 or idx + length > len(self.template):
            raise ValueError(
                f"Deletion region ({start}..{start + length - 1}) "
                f"exceeds template length."
            )

        deleted = self.template[idx:idx + length]
        mutated = self.template[:idx] + self.template[idx + length:]

        description = f"Delete {deleted.lower()} ({length} nt) at position {start}"

        # FWD primer: starts right after deletion, extends 3'
        fwd_seq, fwd_tm = _build_primer(
            mutated, idx, "forward",
            min_tm=self.min_tm, max_tm=self.max_tm, min_len=self.min_length,
            method=self.method,
        )

        # REV primer: ends right before deletion, extends 5'
        rev_seq, rev_tm = _build_primer(
            mutated, idx, "reverse",
            min_tm=self.min_tm, max_tm=self.max_tm, min_len=self.min_length,
            method=self.method,
        )

        # Balance
        fwd_seq, fwd_tm, rev_seq, rev_tm = _balance_primers(
            fwd_seq, rev_seq, fwd_tm, rev_tm, mutated, idx, idx,
            method=self.method,
        )

        ta = min(fwd_tm, rev_tm) + 1
        ta = min(ta, 72)

        # Check for secondary structure risks
        fwd_risk = self._check_primer_risk(fwd_seq)
        rev_risk = self._check_primer_risk(rev_seq)
        
        warnings = []
        if fwd_risk:
            warnings.append(f"Forward primer risk: {fwd_risk}")
        if rev_risk:
            warnings.append(f"Reverse primer risk: {rev_risk}")

        fwd_primer = SDMPrimer(
            sequence=fwd_seq,
            length=len(fwd_seq),
            tm=math.ceil(fwd_tm),
            gc_pct=round(_gc_frac(fwd_seq) * 100),
            direction="FWD",
            bind_region=fwd_seq,
            mutation_bases="",
            risk=fwd_risk,
        )

        rev_primer = SDMPrimer(
            sequence=rev_seq,
            length=len(rev_seq),
            tm=math.ceil(rev_tm),
            gc_pct=round(_gc_frac(rev_seq) * 100),
            direction="REV",
            bind_region=rev_seq,
            mutation_bases="",
            risk=rev_risk,
        )

        return MutagenesisResult(
            description=description,
            forward=fwd_primer,
            reverse=rev_primer,
            ta=ta,
            mutated_sequence=mutated,
            warnings=warnings,
        )

    def insertion(self, position: int, insert_seq: str) -> MutagenesisResult:
        """Design primers for a nucleotide insertion.

        Parameters
        ----------
        position : int
            1-based position after which to insert.
        insert_seq : str
            Sequence to insert.
        """
        insert_seq = insert_seq.upper()
        idx = position  # insert after this 0-based position

        mutated = (
            self.template[:idx]
            + insert_seq
            + self.template[idx:]
        )

        description = f"Insert {insert_seq} ({len(insert_seq)} nt) after position {position}"

        # The insertion is a 5' flap on the FWD primer
        # Annealing starts from the insertion point going 3'
        anneal_start = idx + len(insert_seq)

        fwd_anneal, fwd_anneal_tm = _build_primer(
            mutated, anneal_start, "forward",
            min_tm=self.min_tm, max_tm=self.max_tm, min_len=self.min_length,
            method=self.method,
        )

        # Full FWD = insertion flap + annealing
        fwd_full = insert_seq.lower() + fwd_anneal

        # REV primer: right before insertion point
        rev_seq, rev_tm = _build_primer(
            mutated, idx, "reverse",
            min_tm=self.min_tm, max_tm=self.max_tm, min_len=self.min_length,
            method=self.method,
        )

        ta = min(fwd_anneal_tm, rev_tm) + 1
        ta = min(ta, 72)

        # Check for secondary structure risks
        fwd_risk = self._check_primer_risk(fwd_full)
        rev_risk = self._check_primer_risk(rev_seq)
        
        warnings = []
        if fwd_risk:
            warnings.append(f"Forward primer risk: {fwd_risk}")
        if rev_risk:
            warnings.append(f"Reverse primer risk: {rev_risk}")

        fwd_primer = SDMPrimer(
            sequence=fwd_full,
            length=len(fwd_full),
            tm=math.ceil(fwd_anneal_tm),
            gc_pct=round(_gc_frac(fwd_full) * 100),
            direction="FWD",
            bind_region=fwd_anneal,
            mutation_bases=insert_seq.lower(),
            risk=fwd_risk,
        )

        rev_primer = SDMPrimer(
            sequence=rev_seq,
            length=len(rev_seq),
            tm=math.ceil(rev_tm),
            gc_pct=round(_gc_frac(rev_seq) * 100),
            direction="REV",
            bind_region=rev_seq,
            mutation_bases="",
            risk=rev_risk,
        )

        return MutagenesisResult(
            description=description,
            forward=fwd_primer,
            reverse=rev_primer,
            ta=ta,
            mutated_sequence=mutated,
            warnings=warnings,
        )

    def batch(self, mutations: str) -> list[MutagenesisResult]:
        """Process multiple mutations.

        Parameters
        ----------
        mutations : str
            Space/comma/semicolon-separated mutation list.
            Supports AA point mutations (M1A, K2R) and mixed types.
        """
        parsed = parse_aa_mutation(mutations)
        results = []
        for mut in parsed:
            result = self._design_point_mutation(mut)
            results.append(result)
        return results

    def _design_substitution_primers(
        self,
        mutated_template: str,
        mut_start: int,
        mut_len: int,
        mut_bases: str,
        method: Optional[str] = None,
    ) -> tuple[SDMPrimer, SDMPrimer]:
        """Design primers for a substitution (codon change or nt substitution).

        In the NEB Base Changer model, primers are BACK-TO-BACK:
        - FWD primer: contains the mutation, flanked by annealing bases
        - REV primer: reverse complement, binding immediately adjacent
          to the FWD primer's 5' end (upstream on template)

        For circular templates, we create a doubled-template so primers
        can wrap around the junction (matching NEB Base Changer behavior).
        """
        orig_len = len(mutated_template)

        if self.confine_to_tails or mut_len > 6:
            return self._design_tail_primers(
                mutated_template, mut_start, mut_len, mut_bases, method=method
            )

        # For circular templates, double the sequence so primers can
        # wrap around the junction. Place mutation in the second copy.
        if self.circular:
            circ = mutated_template + mutated_template
            circ_mut_start = orig_len + mut_start  # mutation in 2nd copy
        else:
            circ = mutated_template
            circ_mut_start = mut_start

        tlen = len(circ)
        min_flank = Q5_SDM_MIN_ANNEAL

        # FWD primer: mutation flanked by annealing bases
        left_start = circ_mut_start - min_flank
        if not self.circular:
            left_start = max(0, left_start)
        right_end = circ_mut_start + mut_len

        # NEB computes FWD Tm on the right-side annealing region ONLY:
        # from the mutation start to the 3' end of the primer.
        # The left flank (before mutation) is non-annealing extension.
        fwd_right_anneal = ""
        fwd_tm_val = 0.0
        for end in range(right_end + min_flank, min(right_end + 60, tlen) + 1):
            full_candidate = circ[left_start:end]
            # Tm is computed on right-side only (mutation + right flank)
            right_region = circ[circ_mut_start:end]
            try:
                if method:
                    from .core import tm
                    t = tm(right_region, salt_mM=int(Q5_SDM_MONO_MM), method=method)
                else:
                    t = calc_sdm_tm(right_region)
            except ValueError:
                continue
            if t >= self.min_tm:
                fwd_right_anneal = full_candidate
                fwd_tm_val = t
                if right_region[-1] in "GC":
                    break
                if t > self.max_tm:
                    break

        if not fwd_right_anneal:
            fwd_right_anneal = circ[left_start:min(right_end + 40, tlen)]
            right_region = circ[circ_mut_start:left_start + len(fwd_right_anneal)]
            try:
                if method:
                    from .core import tm
                    fwd_tm_val = tm(right_region, salt_mM=int(Q5_SDM_MONO_MM), method=method) if right_region else 0.0
                else:
                    fwd_tm_val = calc_sdm_tm(right_region) if right_region else 0.0
            except ValueError:
                fwd_tm_val = 0.0

        # Mark mutation bases as lowercase in FWD display
        mut_offset = circ_mut_start - left_start
        fwd_display = (
            fwd_right_anneal[:mut_offset]
            + fwd_right_anneal[mut_offset:mut_offset + mut_len].lower()
            + fwd_right_anneal[mut_offset + mut_len:]
        )

        # REV primer: binds immediately upstream of FWD's 5' end
        if self.circular or left_start >= self.min_length:
            rev_seq, rev_tm = _build_primer(
                circ, left_start, "reverse",
                min_tm=self.min_tm, max_tm=self.max_tm,
                min_len=self.min_length,
                method=self.method,
            )
        else:
            # Linear mode: not enough room upstream — REV binds downstream
            fwd_end_pos = left_start + len(fwd_right_anneal)
            rev_seq, rev_tm = _build_primer(
                circ, fwd_end_pos, "forward",
                min_tm=self.min_tm, max_tm=self.max_tm,
                min_len=self.min_length,
                method=self.method,
            )
            rev_seq = reverse_complement(rev_seq)
            try:
                rev_tm = calc_sdm_tm(rev_seq)
            except ValueError:
                rev_tm = 0.0

        fwd_primer = SDMPrimer(
            sequence=fwd_display,
            length=len(fwd_right_anneal),
            tm=math.ceil(fwd_tm_val),
            gc_pct=round(_gc_frac(fwd_right_anneal) * 100),
            direction="FWD",
            bind_region=fwd_right_anneal,
            mutation_bases=mut_bases.lower(),
        )

        rev_primer = SDMPrimer(
            sequence=rev_seq,
            length=len(rev_seq),
            tm=math.ceil(rev_tm),
            gc_pct=round(_gc_frac(rev_seq) * 100),
            direction="REV",
            bind_region=rev_seq,
            mutation_bases="",
        )

        return fwd_primer, rev_primer

    def _design_tail_primers(
        self,
        mutated_template: str,
        mut_start: int,
        mut_len: int,
        mut_bases: str,
        method: Optional[str] = None,
    ) -> tuple[SDMPrimer, SDMPrimer]:
        """Design primers with mutation as 5' non-annealing tails.

        Used for substitutions >6 nt or when confine_to_tails is enabled.
        """
        tlen = len(mutated_template)

        # Split mutation between FWD and REV 5' tails
        half = mut_len // 2
        fwd_tail = mut_bases[:half].lower() if half > 0 else ""
        rev_tail_bases = mut_bases[half:].lower() if half < mut_len else ""

        # FWD annealing: starts after the mutation region, extends 3'
        fwd_anneal_start = mut_start + mut_len
        fwd_anneal, fwd_tm = _build_primer(
            mutated_template, fwd_anneal_start, "forward",
            min_tm=self.min_tm, max_tm=self.max_tm, min_len=self.min_length,
            method=self.method,
        )

        # REV annealing: ends before the mutation region, extends 5'
        rev_anneal_end = mut_start
        rev_anneal, rev_tm = _build_primer(
            mutated_template, rev_anneal_end, "reverse",
            min_tm=self.min_tm, max_tm=self.max_tm, min_len=self.min_length,
            method=self.method,
        )

        # Full primer sequences = tail + annealing
        fwd_full = fwd_tail + fwd_anneal
        # REV tail needs to be revcomp of the remaining mutation
        rev_tail_rc = reverse_complement(rev_tail_bases.upper()).lower() if rev_tail_bases else ""
        rev_full = rev_tail_rc + rev_anneal

        fwd_primer = SDMPrimer(
            sequence=fwd_full,
            length=len(fwd_full),
            tm=math.ceil(fwd_tm),
            gc_pct=round(_gc_frac(fwd_full) * 100),
            direction="FWD",
            bind_region=fwd_anneal,
            mutation_bases=fwd_tail,
        )

        rev_primer = SDMPrimer(
            sequence=rev_full,
            length=len(rev_full),
            tm=math.ceil(rev_tm),
            gc_pct=round(_gc_frac(rev_full) * 100),
            direction="REV",
            bind_region=rev_anneal,
            mutation_bases=rev_tail_rc,
        )

        return fwd_primer, rev_primer
