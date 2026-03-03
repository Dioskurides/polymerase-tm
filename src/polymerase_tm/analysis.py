"""Primer analysis: dimer check, Gibson overlaps, restriction scan, quality score.

Contains RESTRICTION_ENZYMES database (~120 NEB enzymes),
primer_dimer, gibson_overlaps, restriction_scan, and primer_quality.
"""

from __future__ import annotations

import re
from typing import Optional

from .core import tm
from .dmso import gc_content, primer_hairpin, _reverse_complement


# IUPAC ambiguity code to regex character class mapping
_IUPAC_MAP: dict[str, str] = {
    "A": "A", "T": "T", "G": "G", "C": "C",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "M": "[AC]", "K": "[GT]", "B": "[CGT]", "V": "[ACG]",
    "D": "[AGT]", "H": "[ACT]", "N": "[ATGC]",
}


def _iupac_to_regex(site: str) -> str:
    """Convert an IUPAC DNA recognition site to a regex pattern."""
    return "".join(_IUPAC_MAP.get(c, re.escape(c)) for c in site.upper())


def _has_degenerate(site: str) -> bool:
    """Check if a recognition site contains degenerate IUPAC bases."""
    return bool(set(site.upper()) - set("ATGC"))


# Comprehensive restriction enzyme recognition sites (5'->3')
# All common NEB enzymes, including those with degenerate IUPAC bases
# (N, R, Y, S, W, M, K, etc.).  Degenerate sites are matched via
# regex expansion at scan time.
RESTRICTION_ENZYMES: dict[str, str] = {
    # --- 4-cutters ---
    "AluI": "AGCT",
    "CviAII": "CATG",
    "DpnI": "GATC",
    "DpnII": "GATC",
    "FatI": "CATG",
    "HaeIII": "GGCC",
    "HhaI": "GCGC",
    "HpyCH4V": "TGCA",
    "MboI": "GATC",
    "MluCI": "AATT",
    "MseI": "TTAA",
    "MspI": "CCGG",
    "NlaIII": "CATG",
    "RsaI": "GTAC",
    "Sau3AI": "GATC",
    "TaqI": "TCGA",
    # --- 5-cutters ---
    "AvaII": "GGACC",
    "EcoRII": "CCAGG",
    # --- 6-cutters (most common) ---
    "AatII": "GACGTC",
    "AccI": "GTCGAC",
    "AclI": "AACGTT",
    "AfeI": "AGCGCT",
    "AflII": "CTTAAG",
    "AgeI": "ACCGGT",
    "ApaI": "GGGCCC",
    "ApaLI": "GTGCAC",
    "AscI": "GGCGCGCC",
    "AseI": "ATTAAT",
    "AvaI": "CTCGAG",
    "AvrII": "CCTAGG",
    "BamHI": "GGATCC",
    "BclI": "TGATCA",
    "BglI": "GCCGC",
    "BglII": "AGATCT",
    "BlpI": "GCTNAGC",
    "BmtI": "GCTAGC",
    "BsaBI": "GATC",
    "BsiWI": "CGTACG",
    "BspEI": "TCCGGA",
    "BspHI": "TCATGA",
    "BsrGI": "TGTACA",
    "BssHII": "GCGCGC",
    "BstBI": "TTCGAA",
    "BstEII": "GGTNACC",
    "BstXI": "CCANNNNTGG",
    "BstZ17I": "GTATAC",
    "ClaI": "ATCGAT",
    "DraI": "TTTAAA",
    "EagI": "CGGCCG",
    "EcoNI": "CCTNNNNNAGG",
    "EcoRI": "GAATTC",
    "EcoRV": "GATATC",
    "FseI": "GGCCGGCC",
    "FspI": "TGCGCA",
    "HincII": "GTYRAC",
    "HindIII": "AAGCTT",
    "HpaI": "GTTAAC",
    "KasI": "GGCGCC",
    "KpnI": "GGTACC",
    "MfeI": "CAATTG",
    "MluI": "ACGCGT",
    "MscI": "TGGCCA",
    "NaeI": "GCCGGC",
    "NarI": "GGCGCC",
    "NciI": "CCSGG",
    "NcoI": "CCATGG",
    "NdeI": "CATATG",
    "NgoMIV": "GCCGGC",
    "NheI": "GCTAGC",
    "NotI": "GCGGCCGC",
    "NruI": "TCGCGA",
    "NsiI": "ATGCAT",
    "NspI": "CATG",
    "PacI": "TTAATTAA",
    "PciI": "ACATGT",
    "PmeI": "GTTTAAAC",
    "PmlI": "CACGTG",
    "PpuMI": "AGGCCT",
    "PsiI": "TTATAA",
    "PspOMI": "GGGCCC",
    "PstI": "CTGCAG",
    "PvuI": "CGATCG",
    "PvuII": "CAGCTG",
    "SacI": "GAGCTC",
    "SacII": "CCGCGG",
    "SalI": "GTCGAC",
    "SbfI": "CCTGCAGG",
    "ScaI": "AGTACT",
    "SfoI": "GGCGCC",
    "SgrAI": "CRCCGGYG",
    "SmaI": "CCCGGG",
    "SmlI": "CTYRAG",
    "SnaBI": "TACGTA",
    "SpeI": "ACTAGT",
    "SphI": "GCATGC",
    "SspI": "AATATT",
    "StuI": "AGGCCT",
    "SwaI": "ATTTAAAT",
    "TliI": "CTCGAG",
    "XbaI": "TCTAGA",
    "XhoI": "CTCGAG",
    "XmaI": "CCCGGG",
    # --- 8-cutters (rare cutters) ---
    "AsiSI": "GCGATCGC",
    # --- Common cloning enzymes ---
    "SfiI": "GGCCNNNNNGGCC",
    "BsaI": "GGTCTC",
    "BbsI": "GAAGAC",
    "BsmBI": "CGTCTC",
    "SapI": "GAAGAGC",
    "BtgZI": "GCGATG",
    "EarI": "CTCTTC",
    "BspMI": "ACCTGC",
    "PaqCI": "CACCTGC",
}


def primer_dimer(
    fwd: str,
    rev: str,
    min_complementarity: int = 4,
) -> dict:
    """Check 3'-end complementarity between two primers (primer dimer risk).

    Slides the reverse complement of each primer along the other and
    scores consecutive base-pair matches at the 3' ends.

    Parameters
    ----------
    fwd, rev : str
        Primer sequences.
    min_complementarity : int
        Minimum consecutive 3'-end matches to report (default 4).

    Returns
    -------
    dict
        Keys: ``max_score`` (int), ``risk_level`` (str),
        ``fwd_3prime_matches`` (int), ``rev_3prime_matches`` (int),
        ``details`` (list of dicts with match info).

    Examples
    --------
    >>> result = primer_dimer("ATCGATCGATCG", "CGATCGATCGAT")
    >>> print(result["risk_level"])
    """
    fwd = fwd.strip().upper()
    rev = rev.strip().upper()
    rev_rc = _reverse_complement(rev)
    fwd_rc = _reverse_complement(fwd)

    def _count_3prime_matches(a: str, b_rc: str) -> list[dict]:
        """Find consecutive matches at 3' end of `a` with `b_rc`."""
        matches = []
        len_a, len_b = len(a), len(b_rc)
        # Slide b_rc along a: offset = position where b_rc[0] aligns with a[offset]
        for offset in range(-(len_b - 1), len_a):
            consec = 0
            best_consec = 0
            last_match_pos = -1
            for i in range(max(0, offset), min(len_a, offset + len_b)):
                j = i - offset
                if a[i] == b_rc[j]:
                    consec += 1
                    if consec > best_consec:
                        best_consec = consec
                        last_match_pos = i
                else:
                    consec = 0
            # Check if 3' end of a is involved
            if best_consec >= min_complementarity and last_match_pos >= len_a - 3:
                matches.append({
                    "complementary_bases": best_consec,
                    "position_in_primer": last_match_pos - best_consec + 1,
                    "at_3prime": last_match_pos >= len_a - best_consec,
                })
        return matches

    fwd_matches = _count_3prime_matches(fwd, rev_rc)
    rev_matches = _count_3prime_matches(rev, fwd_rc)

    # Also check self-dimers (3' of fwd with itself, 3' of rev with itself)
    fwd_self = _count_3prime_matches(fwd, fwd_rc)
    rev_self = _count_3prime_matches(rev, rev_rc)

    all_matches = fwd_matches + rev_matches
    max_score = max((m["complementary_bases"] for m in all_matches), default=0)

    # Determine risk level
    if max_score >= 8:
        risk = "high"
    elif max_score >= 6:
        risk = "moderate"
    elif max_score >= min_complementarity:
        risk = "low"
    else:
        risk = "none"

    return {
        "max_score": max_score,
        "risk_level": risk,
        "fwd_3prime_matches": max(
            (m["complementary_bases"] for m in fwd_matches), default=0
        ),
        "rev_3prime_matches": max(
            (m["complementary_bases"] for m in rev_matches), default=0
        ),
        "fwd_self_dimer": max(
            (m["complementary_bases"] for m in fwd_self), default=0
        ),
        "rev_self_dimer": max(
            (m["complementary_bases"] for m in rev_self), default=0
        ),
        "details": all_matches,
    }


def gibson_overlaps(
    fwd_bind: str,
    rev_bind: str,
    left_seq: str,
    right_seq: str,
    overlap_len: int = 20,
    polymerase: str = "q5",
) -> dict:
    """Design Gibson Assembly / HiFi Assembly primers with overlaps.

    Parameters
    ----------
    fwd_bind : str
        Forward primer binding region (anneals to template).
    rev_bind : str
        Reverse primer binding region (anneals to template).
    left_seq : str
        Sequence of the left neighboring fragment (for reverse primer overlap).
    right_seq : str
        Sequence of the right neighboring fragment (for forward primer overlap).
    overlap_len : int
        Overlap length in bp (default 20, NEB recommends 15-25).
    polymerase : str
        Polymerase for Tm calculation of binding regions.

    Returns
    -------
    dict
        Keys: ``fwd_full``, ``rev_full`` (complete primers with overhangs),
        ``fwd_overlap``, ``rev_overlap``, ``fwd_bind_tm``, ``rev_bind_tm``,
        ``fwd_overlap_tm``, ``rev_overlap_tm``.

    Examples
    --------
    >>> result = gibson_overlaps(
    ...     fwd_bind="ATGTCCCTGCTCTTCTCTCGATGCAA",
    ...     rev_bind="GTGCCTCCGAGCCAGCACC",
    ...     left_seq="AATTCCGGTTAACCDDEE" * 3,
    ...     right_seq="GCGATCGATCGATCGATCG" * 3,
    ... )
    >>> print(f"Forward: {result['fwd_full']}")
    """
    fwd_bind = fwd_bind.strip().upper()
    rev_bind = rev_bind.strip().upper()
    left_seq = left_seq.strip().upper()
    right_seq = right_seq.strip().upper()

    # Forward primer: overlap from RIGHT fragment's 5' end + fwd binding
    fwd_overlap = right_seq[:overlap_len]
    fwd_full = fwd_overlap + fwd_bind

    # Reverse primer: overlap from LEFT fragment's 3' end (reverse complement) + rev binding
    left_end = left_seq[-overlap_len:]
    rev_overlap = _reverse_complement(left_end)
    rev_full = rev_overlap + rev_bind

    # Calculate Tm for all parts
    fwd_bind_tm = tm(fwd_bind, polymerase=polymerase)
    rev_bind_tm = tm(rev_bind, polymerase=polymerase)

    # Overlap Tm at 50 nM (typical Gibson concentration)
    fwd_overlap_tm = tm(fwd_overlap, polymerase=polymerase)
    rev_overlap_tm = tm(rev_overlap, polymerase=polymerase)

    return {
        "fwd_full": fwd_full,
        "rev_full": rev_full,
        "fwd_bind": fwd_bind,
        "rev_bind": rev_bind,
        "fwd_overlap": fwd_overlap,
        "rev_overlap": rev_overlap,
        "fwd_length": len(fwd_full),
        "rev_length": len(rev_full),
        "overlap_len": overlap_len,
        "fwd_bind_tm": fwd_bind_tm,
        "rev_bind_tm": rev_bind_tm,
        "fwd_overlap_tm": fwd_overlap_tm,
        "rev_overlap_tm": rev_overlap_tm,
    }


def restriction_scan(
    seq: str,
    enzymes: "Optional[dict[str, str] | list[str]]" = None,
) -> list[dict]:
    """Scan a DNA sequence for restriction enzyme recognition sites.

    Supports IUPAC ambiguity codes (N, R, Y, S, W, M, K, B, V, D, H)
    in recognition sites — degenerate sites are matched via regex.

    Parameters
    ----------
    seq : str
        DNA sequence to scan (primer or amplicon).
    enzymes : list[str] | dict | None
        - ``None``: scan all ~120 built-in NEB enzymes (default)
        - ``list[str]``: enzyme names to check, e.g. ``["EcoRI", "BamHI"]``
        - ``dict``: custom ``{name: recognition_site}`` mapping

    Returns
    -------
    list[dict]
        Each dict: ``enzyme``, ``site``, ``position`` (0-based),
        ``strand`` ("+" or "-").

    Examples
    --------
    >>> hits = restriction_scan("ATGAATTCGATCG")
    >>> restriction_scan("ATGAATTCGATCG", enzymes=["EcoRI", "BamHI"])
    >>> restriction_scan("ATGAATTCGATCG", enzymes={"Custom": "AATTC"})
    """
    seq = seq.strip().upper()
    rc = _reverse_complement(seq)

    # Resolve enzymes argument
    if enzymes is None:
        enzyme_dict = RESTRICTION_ENZYMES
    elif isinstance(enzymes, list):
        # Look up enzyme names (case-insensitive)
        name_map = {k.upper(): k for k in RESTRICTION_ENZYMES}
        enzyme_dict = {}
        for name in enzymes:
            key = name.upper()
            if key in name_map:
                original = name_map[key]
                enzyme_dict[original] = RESTRICTION_ENZYMES[original]
            else:
                raise ValueError(
                    f"Unknown enzyme '{name}'. Use RESTRICTION_ENZYMES.keys() "
                    f"to see available enzymes, or pass a dict for custom sites."
                )
    else:
        enzyme_dict = enzymes

    hits = []
    for name, site in enzyme_dict.items():
        site = site.upper()
        degenerate = _has_degenerate(site)

        if degenerate:
            # Use regex matching for IUPAC-degenerate sites
            pattern = _iupac_to_regex(site)
            for m in re.finditer(f"(?=({pattern}))", seq):
                hits.append({
                    "enzyme": name,
                    "site": m.group(1),
                    "position": m.start(),
                    "strand": "+",
                })
            # Reverse strand
            site_rc = _reverse_complement(site)
            if site_rc != site:
                pattern_rc = _iupac_to_regex(site_rc)
                for m in re.finditer(f"(?=({pattern_rc}))", seq):
                    hits.append({
                        "enzyme": name,
                        "site": m.group(1),
                        "position": m.start(),
                        "strand": "-",
                    })
        else:
            # Fast exact matching for unambiguous sites
            pos = 0
            while True:
                idx = seq.find(site, pos)
                if idx == -1:
                    break
                hits.append({
                    "enzyme": name,
                    "site": site,
                    "position": idx,
                    "strand": "+",
                })
                pos = idx + 1

            # Reverse strand (check if site is palindromic first)
            site_rc = _reverse_complement(site)
            if site_rc != site:
                pos = 0
                while True:
                    idx = seq.find(site_rc, pos)
                    if idx == -1:
                        break
                    hits.append({
                        "enzyme": name,
                        "site": site_rc,
                        "position": idx,
                        "strand": "-",
                    })
                    pos = idx + 1

    return sorted(hits, key=lambda x: x["position"])


def primer_quality(seq: str) -> dict:
    """Comprehensive primer quality assessment with a 0-100 score.

    Evaluates: length, GC content, GC clamp, homopolymer runs,
    dinucleotide repeats, hairpins, and 3' stability.

    Parameters
    ----------
    seq : str
        Primer sequence.

    Returns
    -------
    dict
        Keys: ``score`` (0-100), ``grade`` (A/B/C/D/F),
        ``gc_pct``, ``gc_clamp``, ``max_run``, ``max_repeat``,
        ``hairpin_count``, ``length_ok``, ``issues`` (list).

    Examples
    --------
    >>> result = primer_quality("ATGTCCCTGCTCTTCTCTCGATGCAA")
    >>> print(f"Score: {result['score']}/100 ({result['grade']})")
    """
    seq = seq.strip().upper()
    issues: list[str] = []
    score = 100

    # --- Length ---
    length = len(seq)
    length_ok = 18 <= length <= 30
    if length < 15:
        score -= 25
        issues.append(f"Too short ({length} nt, ideal 18-30)")
    elif length < 18:
        score -= 10
        issues.append(f"Short ({length} nt, ideal 18-30)")
    elif length > 35:
        score -= 15
        issues.append(f"Too long ({length} nt, ideal 18-30)")
    elif length > 30:
        score -= 5
        issues.append(f"Long ({length} nt, ideal 18-30)")

    # --- GC content ---
    gc = gc_content(seq) * 100
    if gc < 30:
        score -= 20
        issues.append(f"Low GC ({gc:.0f}%, ideal 40-60%)")
    elif gc < 40:
        score -= 10
        issues.append(f"GC slightly low ({gc:.0f}%, ideal 40-60%)")
    elif gc > 70:
        score -= 20
        issues.append(f"High GC ({gc:.0f}%, ideal 40-60%)")
    elif gc > 60:
        score -= 10
        issues.append(f"GC slightly high ({gc:.0f}%, ideal 40-60%)")

    # --- GC clamp (last 2 bases) ---
    last2 = seq[-2:]
    gc_clamp_count = last2.count("G") + last2.count("C")
    gc_clamp = gc_clamp_count >= 1
    if gc_clamp_count == 0:
        score -= 10
        issues.append("No GC clamp at 3' end (last 2 bases are AT/T)")
    elif gc_clamp_count == 2:
        # Strong clamp, slight bonus (no penalty)
        pass

    # --- 3' stability (last 5 bases GC) ---
    last5 = seq[-5:] if len(seq) >= 5 else seq
    gc_3prime = (last5.count("G") + last5.count("C")) / len(last5) * 100
    if gc_3prime > 80:
        score -= 5
        issues.append(f"3' end too GC-rich ({gc_3prime:.0f}%, may cause mispriming)")

    # --- Homopolymer runs ---
    max_run = 1
    current_run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1

    if max_run >= 5:
        score -= 15
        issues.append(f"Long homopolymer run ({max_run}x same base)")
    elif max_run >= 4:
        score -= 5
        issues.append(f"Homopolymer run ({max_run}x same base)")

    # --- Dinucleotide repeats ---
    max_repeat = 0
    for i in range(len(seq) - 3):
        dinuc = seq[i:i + 2]
        repeat_count = 1
        j = i + 2
        while j + 1 < len(seq) and seq[j:j + 2] == dinuc:
            repeat_count += 1
            j += 2
        max_repeat = max(max_repeat, repeat_count)

    if max_repeat >= 4:
        score -= 15
        issues.append(f"Dinucleotide repeat ({max_repeat}x)")
    elif max_repeat >= 3:
        score -= 5
        issues.append(f"Dinucleotide repeat ({max_repeat}x)")

    # --- Hairpins ---
    hairpins = primer_hairpin(seq)
    strong_hairpins = [h for h in hairpins if h["stem_length"] >= 4]
    hairpin_count = len(strong_hairpins)
    if hairpin_count >= 2:
        score -= 15
        issues.append(f"{hairpin_count} strong hairpins detected")
    elif hairpin_count == 1:
        score -= 5
        issues.append("1 hairpin detected")

    # --- Restriction sites ---
    re_sites = restriction_scan(seq)
    if re_sites:
        # Informational, not penalized
        enzymes_found = list(set(h["enzyme"] for h in re_sites))
        issues.append(f"Contains restriction site(s): {', '.join(enzymes_found)}")

    # Clamp score
    score = max(0, min(100, score))

    # Grade
    if score >= 90:
        grade = "A"
    elif score >= 75:
        grade = "B"
    elif score >= 60:
        grade = "C"
    elif score >= 40:
        grade = "D"
    else:
        grade = "F"

    return {
        "sequence": seq,
        "length": length,
        "score": score,
        "grade": grade,
        "gc_pct": round(gc, 1),
        "gc_clamp": gc_clamp,
        "max_run": max_run,
        "max_repeat": max_repeat,
        "hairpin_count": hairpin_count,
        "length_ok": length_ok,
        "issues": issues if issues else ["No issues detected"],
    }
