"""
neb_tm_calculator.py — NEB Tm Calculator, Python Implementation

Reproduces the melting temperature (Tm) and annealing temperature (Ta)
calculations of the NEB Tm Calculator (https://tmcalculator.neb.com/).

Algorithm
---------
- Tm:  SantaLucia (1998) nearest-neighbor model.
- Salt correction:  Owczarzy et al. (2004, 2008) monovalent-ion method.
- Ta:  Polymerase-specific offset + cap (see POLYMERASES table).
- DMSO correction:  -0.6 degC per 1 % DMSO (v/v).

Buffer salt concentrations and Ta rules were extracted from the NEB
Tm Calculator front-end source (main-90540d67e2.js, version 1.16.10).

References
----------
1. SantaLucia J Jr. (1998) PNAS 95:1460-5.
2. Owczarzy R et al. (2004) Biochemistry 43:3537-54.
3. Owczarzy R et al. (2008) Biochemistry 47:5336-53.
"""

from __future__ import annotations

import math
from typing import Optional

# Try to import Biopython (optional, only needed for DMSO analysis with
# GenBank files and reverse-complement operations)
try:
    from Bio.Seq import Seq
    from Bio import SeqIO
    _HAS_BIO = True
except ImportError:
    _HAS_BIO = False

# =====================================================================
# Nearest-Neighbor Parameters — SantaLucia 1998, Table 2
# (dH in kcal/mol, dS in cal/(mol*K))
# =====================================================================

NN_PARAMS: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

# Terminal-base initiation corrections (SantaLucia 1998)
TERMINAL: dict[str, tuple[float, float]] = {
    "A": (2.3, 4.1), "T": (2.3, 4.1),
    "G": (0.1, -2.8), "C": (0.1, -2.8),
}

R = 1.987  # Gas constant, cal/(mol*K)

# =====================================================================
# NEB Buffer Salt Concentrations (effective monovalent cation, mM)
# Extracted from the ``buffs`` table in the NEB Tm Calculator JS source.
# =====================================================================

BUFFERS: dict[str, int] = {
    "standard_taq":   55,
    "thermopol":      40,
    "hemo_klentaq":   70,
    "crimson_taq":    55,
    "longamp":       100,
    "multiplex":      90,
    "phusion_hf":    222,
    "phusion_gc":    222,
    "phusionflex_hf": 222,
    "phusionflex_gc": 222,
    "onetaq_std":     54,
    "onetaq_gc":      80,
    "q5":            150,
    "q5mm":          150,
    "q5u":           170,
    "q5bd":          150,
    "longamphs":     100,
}

# =====================================================================
# Polymerase Product Definitions
#
# Each entry maps a human-readable product key to:
#   buffer       — key into BUFFERS (determines salt correction)
#   conc         — default primer concentration (nM)
#   ta_rule      — one of the Ta calculation strategies (see _calc_ta)
#   ta_cap       — maximum Ta (degC)
#   min_len      — minimum primer length to apply the Ta offset
#   description  — full product name
# =====================================================================

POLYMERASES: dict[str, dict] = {
    # --- Q5 family ---
    "q5": {
        "buffer": "q5", "conc": 500, "ta_rule": "offset",
        "ta_offset": 1, "ta_cap": 72, "min_len": 8,
        "description": "Q5 High-Fidelity DNA Polymerase",
    },
    "q5_hot_start": {
        "buffer": "q5", "conc": 500, "ta_rule": "offset",
        "ta_offset": 1, "ta_cap": 72, "min_len": 8,
        "description": "Q5 Hot Start High-Fidelity DNA Polymerase",
    },
    "q5_master_mix": {
        "buffer": "q5mm", "conc": 500, "ta_rule": "offset",
        "ta_offset": 1, "ta_cap": 72, "min_len": 8,
        "description": "Q5 High-Fidelity 2X Master Mix",
    },
    "q5u_hot_start": {
        "buffer": "q5u", "conc": 500, "ta_rule": "offset",
        "ta_offset": 2, "ta_cap": 72, "min_len": 8,
        "description": "Q5U Hot Start High-Fidelity DNA Polymerase",
    },
    "q5_blood_direct": {
        "buffer": "q5bd", "conc": 500, "ta_rule": "offset",
        "ta_offset": 1, "ta_cap": 72, "min_len": 8,
        "description": "Q5 Blood Direct 2X Master Mix",
    },
    # --- Phusion family ---
    "phusion_hf": {
        "buffer": "phusion_hf", "conc": 500, "ta_rule": "phusion",
        "ta_cap": 72, "min_len": 0,
        "description": "Phusion High-Fidelity DNA Polymerase (HF Buffer)",
    },
    "phusion_gc": {
        "buffer": "phusion_gc", "conc": 500, "ta_rule": "phusion",
        "ta_cap": 72, "min_len": 0,
        "description": "Phusion High-Fidelity DNA Polymerase (GC Buffer)",
    },
    "phusion_flex_hf": {
        "buffer": "phusionflex_hf", "conc": 500, "ta_rule": "phusion",
        "ta_cap": 72, "min_len": 0,
        "description": "Phusion Hot Start Flex DNA Polymerase (HF Buffer)",
    },
    "phusion_flex_gc": {
        "buffer": "phusionflex_gc", "conc": 500, "ta_rule": "phusion",
        "ta_cap": 72, "min_len": 0,
        "description": "Phusion Hot Start Flex DNA Polymerase (GC Buffer)",
    },
    # --- Taq family ---
    "taq": {
        "buffer": "standard_taq", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "Taq DNA Polymerase (Standard Buffer)",
    },
    "taq_thermopol": {
        "buffer": "thermopol", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "Taq DNA Polymerase (ThermoPol Buffer)",
    },
    "hot_start_taq": {
        "buffer": "standard_taq", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "Hot Start Taq DNA Polymerase",
    },
    "crimson_taq": {
        "buffer": "crimson_taq", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "Crimson Taq DNA Polymerase",
    },
    "hemo_klentaq": {
        "buffer": "hemo_klentaq", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "Hemo KlenTaq",
    },
    "epimark": {
        "buffer": "standard_taq", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "EpiMark Hot Start Taq DNA Polymerase",
    },
    # --- OneTaq family ---
    "onetaq_std": {
        "buffer": "onetaq_std", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "OneTaq DNA Polymerase (Standard Buffer)",
    },
    "onetaq_gc": {
        "buffer": "onetaq_gc", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "OneTaq DNA Polymerase (GC Buffer)",
    },
    # --- LongAmp family ---
    "longamp": {
        "buffer": "longamp", "conc": 400, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 65, "min_len": 8,
        "description": "LongAmp Taq DNA Polymerase",
    },
    "longamp_hot_start": {
        "buffer": "longamphs", "conc": 400, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 65, "min_len": 8,
        "description": "LongAmp Hot Start Taq DNA Polymerase",
    },
    # --- Vent / Deep Vent ---
    "vent": {
        "buffer": "thermopol", "conc": 200, "ta_rule": "offset",
        "ta_offset": -2, "ta_cap": 72, "min_len": 21,
        "description": "Vent DNA Polymerase",
    },
    "deep_vent": {
        "buffer": "thermopol", "conc": 200, "ta_rule": "offset",
        "ta_offset": -2, "ta_cap": 72, "min_len": 21,
        "description": "Deep Vent DNA Polymerase",
    },
    # --- Multiplex ---
    "multiplex": {
        "buffer": "multiplex", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "Multiplex PCR 5X Master Mix",
    },
}


# =====================================================================
# Core Tm Calculations
# =====================================================================

def calc_nn_raw(seq: str, primer_conc_nM: float = 500) -> float:
    """SantaLucia 1998 nearest-neighbor Tm without salt correction.

    Parameters
    ----------
    seq : str
        Primer sequence (A/T/G/C only).
    primer_conc_nM : float
        Total primer concentration in nM. Default 500 nM.

    Returns
    -------
    float
        Raw Tm in degrees Celsius (1 M NaCl reference state).
    """
    seq = seq.upper().strip()
    conc_M = primer_conc_nM * 1e-9
    dH, dS = 0.0, 0.0
    for i in range(len(seq) - 1):
        h, s = NN_PARAMS[seq[i:i + 2]]
        dH += h
        dS += s
    for base in (seq[0], seq[-1]):
        h, s = TERMINAL[base]
        dH += h
        dS += s
    return (dH * 1000.0) / (dS + R * math.log(conc_M)) - 273.15


def owczarzy_correction(tm_raw: float, seq: str, salt_mM: float) -> float:
    """Owczarzy et al. (2004) monovalent-ion Tm correction.

    Parameters
    ----------
    tm_raw : float
        Raw Tm from ``calc_nn_raw`` (1 M reference).
    seq : str
        Primer sequence.
    salt_mM : float
        Effective monovalent salt concentration in mM.

    Returns
    -------
    float
        Salt-corrected Tm in degrees Celsius.
    """
    fgc = (seq.upper().count("G") + seq.upper().count("C")) / len(seq)
    salt_M = salt_mM / 1000.0
    ln_s = math.log(salt_M)
    corr = (4.29e-5 * fgc - 3.95e-5) * ln_s + 9.40e-6 * ln_s ** 2
    return 1.0 / (1.0 / (tm_raw + 273.15) + corr) - 273.15


def tm(seq: str, polymerase: str = "q5") -> int:
    """Calculate Tm for a primer sequence with a given NEB polymerase.

    Parameters
    ----------
    seq : str
        Primer binding sequence (A/T/G/C only, no overhangs).
    polymerase : str
        Key into ``POLYMERASES``. Default ``"q5"``.

    Returns
    -------
    int
        Melting temperature rounded to the nearest integer (degC).
    """
    poly = POLYMERASES[polymerase]
    salt = BUFFERS[poly["buffer"]]
    conc = poly["conc"]
    raw = calc_nn_raw(seq, conc)
    return round(owczarzy_correction(raw, seq, salt))


def ta(
    seq1: str,
    seq2: str,
    polymerase: str = "q5",
    dmso_pct: float = 0,
) -> tuple[int, int, int]:
    """Calculate the recommended annealing temperature for a primer pair.

    Parameters
    ----------
    seq1, seq2 : str
        Primer binding sequences (no overhangs).
    polymerase : str
        Key into ``POLYMERASES``.
    dmso_pct : float
        DMSO percentage (v/v). Correction: -0.6 degC per 1 %.

    Returns
    -------
    (Ta, Tm1, Tm2) : tuple[int, int, int]
        Recommended annealing temperature, Tm of seq1, Tm of seq2.
    """
    t1 = tm(seq1, polymerase)
    t2 = tm(seq2, polymerase)
    poly = POLYMERASES[polymerase]
    min_tm = min(t1, t2)
    min_len = min(len(seq1.strip()), len(seq2.strip()))

    if poly["ta_rule"] == "phusion":
        result = 0.93 * min_tm + 7.5
    elif min_len >= poly.get("min_len", 8):
        result = min_tm + poly.get("ta_offset", 0)
    else:
        result = min_tm

    result = min(result, poly["ta_cap"])
    result -= 0.6 * dmso_pct
    return round(result), t1, t2


def list_polymerases() -> list[dict]:
    """Return a summary of all available polymerases.

    Returns
    -------
    list[dict]
        Each dict has keys: key, description, buffer_salt_mM, conc_nM.
    """
    out = []
    for key, poly in POLYMERASES.items():
        out.append({
            "key": key,
            "description": poly["description"],
            "buffer_salt_mM": BUFFERS[poly["buffer"]],
            "primer_conc_nM": poly["conc"],
        })
    return out


# =====================================================================
# DMSO / Secondary-Structure Analysis
# =====================================================================

def gc_content(seq: str) -> float:
    """Return GC content as a fraction (0.0 – 1.0)."""
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)


def gc_windows(
    seq: str, window: int = 50, step: int = 10
) -> list[tuple[int, float]]:
    """Sliding-window GC analysis.

    Returns
    -------
    list of (position, gc_fraction)
    """
    seq = seq.upper()
    return [
        (i, gc_content(seq[i:i + window]))
        for i in range(0, len(seq) - window + 1, step)
    ]


def find_gc_hotspots(
    seq: str, window: int = 50, step: int = 10, threshold: float = 0.70
) -> list[tuple[int, int, float]]:
    """Find contiguous regions where GC content exceeds *threshold*.

    Returns
    -------
    list of (start, end, max_gc_fraction)
    """
    windows = gc_windows(seq, window, step)
    hotspots: list[tuple[int, int, float]] = []
    in_hs = False
    start = 0
    for pos, gc in windows:
        if gc >= threshold and not in_hs:
            start, in_hs = pos, True
        elif gc < threshold and in_hs:
            peak = max(g for p, g in windows if start <= p < pos + window)
            hotspots.append((start, pos + window, peak))
            in_hs = False
    if in_hs:
        peak = max(g for p, g in windows if p >= start)
        hotspots.append((start, len(seq), peak))
    return hotspots


def find_hairpins(
    seq: str,
    stem_min: int = 6,
    loop_min: int = 3,
    loop_max: int = 8,
) -> list[dict]:
    """Detect potential hairpin (stem-loop) structures in *seq*.

    Only perfectly self-complementary stems are reported.

    Returns
    -------
    list of dict
        Keys: position, stem_length, loop_length, stem_seq, stem_gc,
        tm_estimate (Wallace rule).
    """
    seq = seq.upper()
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    hits: list[dict] = []
    n = len(seq)

    for i in range(n):
        for loop_len in range(loop_min, loop_max + 1):
            for stem_len in range(stem_min, min(20, i + 1, n - i - loop_len) + 1):
                j = i + loop_len
                if j + stem_len > n:
                    break
                left_start = i - stem_len + 1
                if left_start < 0:
                    break
                left = seq[left_start:i + 1]
                right = seq[j:j + stem_len]
                if all(comp.get(a) == b for a, b in zip(left, reversed(right))):
                    gc = gc_content(left)
                    tm_est = 2 * (left.count("A") + left.count("T")) + \
                             4 * (left.count("G") + left.count("C"))
                    hits.append({
                        "position": left_start,
                        "stem_length": stem_len,
                        "loop_length": loop_len,
                        "stem_seq": left,
                        "stem_gc": gc,
                        "tm_estimate": tm_est,
                    })

    # Deduplicate: keep strongest hairpin per 5-bp window
    seen: set[int] = set()
    unique: list[dict] = []
    for hp in sorted(hits, key=lambda x: -x["stem_length"]):
        bucket = (hp["position"] // 5) * 5
        if bucket not in seen:
            seen.add(bucket)
            unique.append(hp)
    return unique


def primer_hairpin(seq: str, stem_min: int = 4) -> list[dict]:
    """Convenience wrapper: detect hairpins within a single primer."""
    return find_hairpins(seq, stem_min=stem_min, loop_min=2, loop_max=6)


def _reverse_complement(seq: str) -> str:
    """Simple reverse complement (no Biopython dependency)."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(comp[b] for b in reversed(seq.upper()))


def analyze_amplicon(
    template_seq: str,
    fwd_bind: str,
    rev_bind: str,
) -> tuple[Optional[str], int, int]:
    """Locate an amplicon on *template_seq* given primer binding regions.

    Parameters
    ----------
    template_seq : str
        Full template sequence (linear or circular).
    fwd_bind : str
        Forward primer binding sequence.
    rev_bind : str
        Reverse primer binding sequence.

    Returns
    -------
    (amplicon, fwd_pos, rev_pos)
        *amplicon* is ``None`` when one or both primers are not found.
    """
    template = template_seq.upper()
    fwd = fwd_bind.upper()
    rev_rc = _reverse_complement(rev_bind)

    fwd_pos = template.find(fwd)
    rev_pos = template.find(rev_rc)

    if fwd_pos == -1 or rev_pos == -1:
        return None, fwd_pos, rev_pos

    if fwd_pos < rev_pos:
        amplicon = template[fwd_pos:rev_pos + len(rev_rc)]
    else:
        amplicon = template[fwd_pos:] + template[:rev_pos + len(rev_rc)]

    return amplicon, fwd_pos, rev_pos


def dmso_recommendation(
    fwd_bind: str,
    rev_bind: str,
    template_seq: Optional[str] = None,
    template_file: Optional[str] = None,
) -> dict:
    """Analyse whether DMSO should be added to the PCR reaction.

    Checks performed
    ----------------
    1. Primer self-hairpins (stem >= 5 bp, Wallace Tm >= 40 degC).
    2. Amplicon overall GC content (> 65 % triggers recommendation).
    3. GC-rich hotspots in amplicon (sliding window > 70 %).
    4. Strong hairpins in amplicon (stem >= 8 bp, Tm >= 50 degC).

    Parameters
    ----------
    fwd_bind, rev_bind : str
        Primer binding sequences (no Gibson overhangs).
    template_seq : str, optional
        Raw template sequence.
    template_file : str, optional
        Path to a GenBank file (requires Biopython).

    Returns
    -------
    dict
        Keys: dmso_recommended (bool), dmso_pct (int),
        reasons (list[str]), primer_analysis (dict),
        amplicon_analysis (dict).
    """
    result: dict = {
        "dmso_recommended": False,
        "dmso_pct": 0,
        "reasons": [],
        "primer_analysis": {},
        "amplicon_analysis": {},
    }

    # --- Primer hairpins ---
    for label, seq in [("fwd", fwd_bind), ("rev", rev_bind)]:
        hps = primer_hairpin(seq)
        strong = [h for h in hps if h["stem_length"] >= 5 and h["tm_estimate"] >= 40]
        result["primer_analysis"][label] = {
            "sequence": seq,
            "length": len(seq),
            "gc_pct": round(gc_content(seq) * 100, 1),
            "hairpins_total": len(hps),
            "strong_hairpins": len(strong),
            "details": strong,
        }
        if strong:
            result["dmso_recommended"] = True
            result["reasons"].append(
                f"Primer {label} has {len(strong)} strong hairpin(s) "
                f"(stem >= 5 bp, Tm >= 40 degC): {strong[0]['stem_seq']}"
            )

    # --- Amplicon analysis ---
    if template_file and not template_seq:
        if not _HAS_BIO:
            raise ImportError(
                "Biopython is required for GenBank file reading. "
                "Install with: pip install biopython"
            )
        record = SeqIO.read(template_file, "genbank")
        template_seq = str(record.seq)

    if template_seq:
        amplicon, fwd_pos, rev_pos = analyze_amplicon(
            template_seq, fwd_bind, rev_bind
        )
        if amplicon:
            amp_gc = gc_content(amplicon)
            hotspots = find_gc_hotspots(amplicon, window=30, step=5, threshold=0.70)
            amp_hps = find_hairpins(amplicon, stem_min=8)
            strong_hp = [h for h in amp_hps if h["tm_estimate"] >= 50]

            result["amplicon_analysis"] = {
                "length": len(amplicon),
                "gc_pct": round(amp_gc * 100, 1),
                "gc_hotspots": len(hotspots),
                "gc_hotspot_details": [
                    (s, e, round(g * 100, 1)) for s, e, g in hotspots
                ],
                "hairpins_total": len(amp_hps),
                "strong_hairpins": len(strong_hp),
            }

            if amp_gc > 0.65:
                result["dmso_recommended"] = True
                result["dmso_pct"] = 3
                result["reasons"].append(
                    f"Amplicon GC content is {amp_gc * 100:.1f} % (> 65 %)"
                )
            if hotspots:
                result["dmso_recommended"] = True
                result["dmso_pct"] = max(result["dmso_pct"], 3)
                result["reasons"].append(
                    f"{len(hotspots)} GC-rich region(s) (> 70 %) in amplicon"
                )
            if strong_hp:
                result["dmso_recommended"] = True
                result["dmso_pct"] = max(result["dmso_pct"], 3)
                result["reasons"].append(
                    f"{len(strong_hp)} strong hairpin(s) in amplicon "
                    f"(stem >= 8 bp, Tm >= 50 degC)"
                )
        else:
            result["amplicon_analysis"] = {
                "error": "Primers not found on template",
                "fwd_pos": fwd_pos,
                "rev_pos": rev_pos,
            }

    if not result["dmso_recommended"]:
        result["reasons"].append("No structural issues detected -- DMSO not required")

    return result


def print_dmso_report(report: dict) -> None:
    """Print a human-readable DMSO recommendation report."""
    print("=" * 60)
    print("DMSO RECOMMENDATION")
    print("=" * 60)

    if report["dmso_recommended"]:
        print(f"\n  [!] DMSO RECOMMENDED: {report['dmso_pct']} %")
    else:
        print(f"\n  [OK] DMSO NOT REQUIRED")

    print(f"\n  Reasons:")
    for r in report["reasons"]:
        print(f"    - {r}")

    pa = report["primer_analysis"]
    print(f"\n  Primer analysis:")
    for label in ("fwd", "rev"):
        p = pa[label]
        hp = f", {p['strong_hairpins']} strong" if p["strong_hairpins"] else ""
        print(
            f"    {label.upper()}: {p['length']} bp, GC = {p['gc_pct']} %, "
            f"{p['hairpins_total']} hairpins{hp}"
        )

    aa = report.get("amplicon_analysis", {})
    if aa and "length" in aa:
        print(f"\n  Amplicon analysis:")
        print(f"    Length: {aa['length']} bp, GC: {aa['gc_pct']} %")
        print(f"    GC hotspots (> 70 %): {aa['gc_hotspots']}")
        for s, e, g in aa.get("gc_hotspot_details", []):
            print(f"      Position {s}-{e}: {g} % GC")
        print(
            f"    Hairpins (stem >= 8 bp): {aa['hairpins_total']}"
            f" ({aa['strong_hairpins']} strong)"
        )
    print()



# =====================================================================
# Automation & Batch Processing
# =====================================================================

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence (A/T/G/C, case-insensitive).

    Returns
    -------
    str
        Reverse complement in uppercase.

    Examples
    --------
    >>> reverse_complement("ATCGATCG")
    'CGATCGAT'
    >>> reverse_complement("AAGCTTT")
    'AAAGCTT'
    """
    return _reverse_complement(seq)


def batch_tm(
    sequences: list[str],
    polymerase: str = "q5",
) -> list[dict]:
    """Calculate Tm for multiple primer sequences at once.

    Parameters
    ----------
    sequences : list[str]
        List of primer binding sequences.
    polymerase : str
        Polymerase key (default ``"q5"``).

    Returns
    -------
    list[dict]
        Each dict has keys: ``sequence``, ``length``, ``gc_pct``, ``tm``.

    Examples
    --------
    >>> results = batch_tm(["ATCGATCGATCG", "GCGCGCGCGCGC"])
    >>> for r in results:
    ...     print(f"{r['sequence']}: Tm={r['tm']} degC")
    """
    results = []
    for seq in sequences:
        s = seq.strip().upper()
        t = tm(s, polymerase=polymerase)
        gc = (s.count("G") + s.count("C")) / len(s) * 100
        results.append({
            "sequence": s,
            "length": len(s),
            "gc_pct": round(gc, 1),
            "tm": t,
        })
    return results


def optimal_binding_length(
    seq: str,
    target_tm: int = 72,
    polymerase: str = "q5",
    direction: str = "3prime",
) -> dict:
    """Find the shortest primer binding region that reaches a target Tm.

    Starting from the 3' or 5' end, progressively extends the binding
    region until the target Tm is reached or the full sequence is used.

    Parameters
    ----------
    seq : str
        Full primer or template sequence to derive the binding region from.
    target_tm : int
        Desired Tm in degC (default 72 for Q5).
    polymerase : str
        Polymerase key (default ``"q5"``).
    direction : str
        ``"3prime"`` (default) extends from the 3' end (standard primer
        design). ``"5prime"`` extends from the 5' end.

    Returns
    -------
    dict
        Keys: ``binding_seq``, ``length``, ``tm``, ``gc_pct``,
        ``target_reached`` (bool), ``full_sequence``.

    Examples
    --------
    >>> result = optimal_binding_length("ATGTCCCTGCTCTTCTCTCGATGCAA", target_tm=65)
    >>> print(f"{result['binding_seq']} ({result['length']} nt, Tm={result['tm']})")
    """
    seq = seq.strip().upper()

    # Minimum 10 nt to get a meaningful Tm
    for length in range(10, len(seq) + 1):
        if direction == "3prime":
            subseq = seq[-length:]
        else:
            subseq = seq[:length]

        t = tm(subseq, polymerase=polymerase)
        if t >= target_tm:
            gc = (subseq.count("G") + subseq.count("C")) / len(subseq) * 100
            return {
                "binding_seq": subseq,
                "length": length,
                "tm": t,
                "gc_pct": round(gc, 1),
                "target_reached": True,
                "full_sequence": seq,
            }

    # Target not reached with full sequence
    subseq = seq
    gc = (subseq.count("G") + subseq.count("C")) / len(subseq) * 100
    return {
        "binding_seq": subseq,
        "length": len(subseq),
        "tm": tm(subseq, polymerase=polymerase),
        "gc_pct": round(gc, 1),
        "target_reached": False,
        "full_sequence": seq,
    }


def check_pair(
    fwd: str,
    rev: str,
    polymerase: str = "q5",
    dmso_pct: float = 0,
) -> dict:
    """Comprehensive primer pair compatibility check.

    Parameters
    ----------
    fwd, rev : str
        Primer binding sequences.
    polymerase : str
        Polymerase key.
    dmso_pct : float
        DMSO percentage.

    Returns
    -------
    dict
        Keys: ``fwd_tm``, ``rev_tm``, ``ta``, ``tm_diff``,
        ``gc_fwd``, ``gc_rev``, ``compatible`` (bool), ``warnings`` (list).

    Examples
    --------
    >>> result = check_pair("ATGTCCCTGCTCTTCTCTCGATGCAA", "GTGCCTCCGAGCCAGCACC")
    >>> if not result["compatible"]:
    ...     for w in result["warnings"]: print(w)
    """
    fwd = fwd.strip().upper()
    rev = rev.strip().upper()

    result_ta, t1, t2 = ta(fwd, rev, polymerase=polymerase, dmso_pct=dmso_pct)
    diff = abs(t1 - t2)
    gc_f = (fwd.count("G") + fwd.count("C")) / len(fwd) * 100
    gc_r = (rev.count("G") + rev.count("C")) / len(rev) * 100

    warnings = []
    if diff > 5:
        warnings.append(
            f"Tm difference is {diff} degC (> 5 degC). Consider redesigning "
            f"primers to reduce the Tm gap."
        )
    if diff > 2:
        warnings.append(
            f"Tm difference is {diff} degC (> 2 degC). Minor mismatch in "
            f"annealing efficiency possible."
        )
    if len(fwd) < 15:
        warnings.append(f"Forward primer is short ({len(fwd)} nt). Risk of non-specific binding.")
    if len(rev) < 15:
        warnings.append(f"Reverse primer is short ({len(rev)} nt). Risk of non-specific binding.")
    if gc_f < 30 or gc_f > 70:
        warnings.append(f"Forward primer GC content ({gc_f:.0f}%) is outside ideal range (30-70%).")
    if gc_r < 30 or gc_r > 70:
        warnings.append(f"Reverse primer GC content ({gc_r:.0f}%) is outside ideal range (30-70%).")

    return {
        "fwd_seq": fwd,
        "rev_seq": rev,
        "fwd_length": len(fwd),
        "rev_length": len(rev),
        "fwd_tm": t1,
        "rev_tm": t2,
        "tm_diff": diff,
        "ta": result_ta,
        "gc_fwd": round(gc_f, 1),
        "gc_rev": round(gc_r, 1),
        "polymerase": polymerase,
        "dmso_pct": dmso_pct,
        "compatible": len(warnings) == 0 or all("minor" in w.lower() for w in warnings),
        "warnings": warnings,
    }


def pcr_protocol(
    fwd: str,
    rev: str,
    polymerase: str = "q5",
    amplicon_length: Optional[int] = None,
    dmso_pct: float = 0,
    num_cycles: int = 30,
) -> dict:
    """Generate a complete PCR cycling protocol.

    Extension times are calculated from the amplicon length and the
    polymerase's extension rate. Denaturation and annealing temperatures
    are polymerase-specific.

    Parameters
    ----------
    fwd, rev : str
        Primer binding sequences.
    polymerase : str
        Polymerase key (default ``"q5"``).
    amplicon_length : int, optional
        Expected amplicon length in bp. If ``None``, extension time
        defaults to 30 s.
    dmso_pct : float
        DMSO percentage.
    num_cycles : int
        Number of PCR cycles (default 30).

    Returns
    -------
    dict
        Full cycling protocol with temperatures and times.

    Examples
    --------
    >>> protocol = pcr_protocol(
    ...     "ATGTCCCTGCTCTTCTCTCGATGCAA",
    ...     "GTGCCTCCGAGCCAGCACC",
    ...     amplicon_length=2500,
    ... )
    >>> for step in protocol["cycling"]:
    ...     print(f"{step['step']}: {step['temp']} degC, {step['time']}")
    """
    # Extension rates (seconds per kb) and denaturation temps
    EXTENSION_RATES = {
        "q5": {"rate_s_per_kb": 30, "denat_temp": 98, "denat_init": 30,
               "denat_cycle": 10, "ext_temp": 72, "final_ext": 120},
        "phusion": {"rate_s_per_kb": 30, "denat_temp": 98, "denat_init": 30,
                    "denat_cycle": 10, "ext_temp": 72, "final_ext": 120},
        "taq": {"rate_s_per_kb": 60, "denat_temp": 95, "denat_init": 30,
                "denat_cycle": 30, "ext_temp": 68, "final_ext": 300},
        "onetaq": {"rate_s_per_kb": 60, "denat_temp": 94, "denat_init": 30,
                   "denat_cycle": 30, "ext_temp": 68, "final_ext": 300},
        "longamp": {"rate_s_per_kb": 50, "denat_temp": 94, "denat_init": 30,
                    "denat_cycle": 30, "ext_temp": 65, "final_ext": 600},
        "vent": {"rate_s_per_kb": 120, "denat_temp": 95, "denat_init": 30,
                 "denat_cycle": 30, "ext_temp": 72, "final_ext": 300},
    }

    # Map polymerase key to family
    poly_key = polymerase.lower()
    if "q5" in poly_key:
        family = "q5"
    elif "phusion" in poly_key:
        family = "phusion"
    elif "longamp" in poly_key:
        family = "longamp"
    elif "onetaq" in poly_key:
        family = "onetaq"
    elif "vent" in poly_key or "deep_vent" in poly_key:
        family = "vent"
    else:
        family = "taq"

    params = EXTENSION_RATES[family]
    result_ta, t1, t2 = ta(fwd, rev, polymerase=polymerase, dmso_pct=dmso_pct)

    # Calculate extension time
    if amplicon_length:
        ext_time = max(10, round(amplicon_length / 1000 * params["rate_s_per_kb"]))
    else:
        ext_time = 30  # default

    def fmt_time(seconds: int) -> str:
        if seconds >= 60:
            m, s = divmod(seconds, 60)
            return f"{m} min {s:02d} s" if s else f"{m} min"
        return f"{seconds} s"

    cycling = [
        {
            "step": "Initial Denaturation",
            "temp": params["denat_temp"],
            "time": fmt_time(params["denat_init"]),
            "seconds": params["denat_init"],
            "cycles": 1,
        },
        {
            "step": "Denaturation",
            "temp": params["denat_temp"],
            "time": fmt_time(params["denat_cycle"]),
            "seconds": params["denat_cycle"],
            "cycles": num_cycles,
        },
        {
            "step": "Annealing",
            "temp": result_ta,
            "time": fmt_time(30),
            "seconds": 30,
            "cycles": num_cycles,
        },
        {
            "step": "Extension",
            "temp": params["ext_temp"],
            "time": fmt_time(ext_time),
            "seconds": ext_time,
            "cycles": num_cycles,
        },
        {
            "step": "Final Extension",
            "temp": params["ext_temp"],
            "time": fmt_time(params["final_ext"]),
            "seconds": params["final_ext"],
            "cycles": 1,
        },
        {
            "step": "Hold",
            "temp": 4,
            "time": "indefinite",
            "seconds": 0,
            "cycles": 1,
        },
    ]

    total_seconds = (
        params["denat_init"]
        + num_cycles * (params["denat_cycle"] + 30 + ext_time)
        + params["final_ext"]
    )

    return {
        "polymerase": polymerase,
        "polymerase_family": family,
        "fwd_tm": t1,
        "rev_tm": t2,
        "ta": result_ta,
        "dmso_pct": dmso_pct,
        "amplicon_length": amplicon_length,
        "extension_rate_s_per_kb": params["rate_s_per_kb"],
        "num_cycles": num_cycles,
        "cycling": cycling,
        "total_time_min": round(total_seconds / 60, 1),
    }


def from_csv(
    path: str,
    polymerase: str = "q5",
    fwd_col: str = "fwd",
    rev_col: str = "rev",
    name_col: Optional[str] = "name",
    sep: str = ",",
) -> list[dict]:
    """Read primer pairs from a CSV file and compute Tm/Ta for each.

    Parameters
    ----------
    path : str
        Path to the CSV file.
    polymerase : str
        Polymerase key.
    fwd_col, rev_col : str
        Column names for forward and reverse primer sequences.
    name_col : str or None
        Optional column name for pair labels.
    sep : str
        CSV delimiter (default ``","``).

    Returns
    -------
    list[dict]
        One dict per row with keys: ``name``, ``fwd``, ``rev``,
        ``fwd_tm``, ``rev_tm``, ``ta``, ``tm_diff``, ``gc_fwd``,
        ``gc_rev``, ``compatible``, ``warnings``.

    Examples
    --------
    >>> results = from_csv("primers.csv", polymerase="q5")
    >>> for r in results:
    ...     print(f"{r['name']}: Ta={r['ta']} degC")
    """
    import csv

    results = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter=sep)
        for i, row in enumerate(reader):
            fwd_seq = row.get(fwd_col, "").strip()
            rev_seq = row.get(rev_col, "").strip()
            name = row.get(name_col, f"pair_{i+1}") if name_col else f"pair_{i+1}"

            if not fwd_seq or not rev_seq:
                results.append({"name": name, "error": "Missing sequence"})
                continue

            pair_result = check_pair(fwd_seq, rev_seq, polymerase=polymerase)
            pair_result["name"] = name
            results.append(pair_result)

    return results


def to_csv(
    results: list[dict],
    path: str,
    sep: str = ",",
) -> None:
    """Write batch results to a CSV file.

    Parameters
    ----------
    results : list[dict]
        Output from ``from_csv``, ``batch_tm``, or any list of dicts.
    path : str
        Output CSV file path.
    sep : str
        Delimiter (default ``","``).

    Examples
    --------
    >>> results = from_csv("input.csv")
    >>> to_csv(results, "output_with_tm.csv")
    """
    import csv

    if not results:
        return

    fieldnames = list(results[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=sep)
        writer.writeheader()
        for row in results:
            # Flatten lists to strings for CSV
            flat = {}
            for k, v in row.items():
                if isinstance(v, list):
                    flat[k] = "; ".join(str(x) for x in v)
                else:
                    flat[k] = v
            writer.writerow(flat)


__version__ = "0.3.0"

__all__ = [
    # Core
    "tm",
    "ta",
    "list_polymerases",
    "calc_nn_raw",
    "owczarzy_correction",
    # Automation
    "reverse_complement",
    "batch_tm",
    "optimal_binding_length",
    "check_pair",
    "pcr_protocol",
    "from_csv",
    "to_csv",
    # DMSO
    "dmso_recommendation",
    "print_dmso_report",
    "gc_content",
    "gc_windows",
    "find_gc_hotspots",
    "find_hairpins",
    "primer_hairpin",
    "analyze_amplicon",
    # Data
    "NN_PARAMS",
    "TERMINAL",
    "BUFFERS",
    "POLYMERASES",
]

