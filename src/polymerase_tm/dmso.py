"""DMSO / secondary-structure analysis.

GC analysis, hairpin detection, amplicon analysis,
and DMSO recommendation logic.
"""

from __future__ import annotations

from typing import Optional

# Biopython is optional (only for GenBank template files)
try:
    from Bio import SeqIO
    _HAS_BIO = True
except ImportError:
    _HAS_BIO = False


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


def _hairpin_nn_tm(stem_seq: str, loop_length: int,
                   salt_mM: float = 50.0) -> float:
    """Compute hairpin Tm using nearest-neighbor thermodynamics.

    Uses SantaLucia (1998) NN parameters with a loop penalty
    of dS = -1.75 * R * ln(loop_length) (Jacobson-Stockmayer).
    Salt correction via Owczarzy (2004) simplified form.

    Parameters
    ----------
    stem_seq : str
        5'->3' sequence of one arm of the hairpin stem.
    loop_length : int
        Number of unpaired bases in the loop.
    salt_mM : float
        Monovalent cation concentration in mM (default 50).

    Returns
    -------
    float
        Estimated Tm in degC.  Returns 0.0 if stem is too short.
    """
    import math
    from .constants import NN_PARAMS, R

    seq = stem_seq.upper()
    if len(seq) < 2:
        return 0.0

    # Accumulate NN dH, dS for the stem duplex
    dH = 0.0  # kcal/mol
    dS = 0.0  # cal/(mol*K)

    for i in range(len(seq) - 1):
        dinuc = seq[i:i + 2]
        if dinuc in NN_PARAMS:
            h, s = NN_PARAMS[dinuc]
            dH += h
            dS += s

    # Initiation for first/last base of stem (hairpin = unimolecular)
    # Use a simplified +0.2 kcal/mol initiation for hairpin stems
    dH += 0.2

    # Loop entropy penalty  (Jacobson-Stockmayer approximation)
    # dS_loop ≈ -1.75 * R * ln(loop_length)  in cal/(mol*K)
    if loop_length >= 1:
        dS += -1.75 * R * math.log(loop_length)

    if dS == 0:
        return 0.0

    # Unimolecular Tm: Tm = dH / dS  (no concentration term)
    # Convert dH from kcal to cal
    tm_kelvin = (dH * 1000.0) / dS

    # Salt correction (simplified Owczarzy 2004 for monovalent ions)
    if salt_mM > 0:
        log_salt = math.log(salt_mM / 1000.0)
        # Owczarzy Eq. 22 simplified: 1/Tm_corrected = 1/Tm + (4.29*f_GC - 3.95)*1e-5*ln[Na+] + 9.40e-6*(ln[Na+])^2
        f_gc = gc_content(seq)
        inv_tm = (1.0 / tm_kelvin
                  + (4.29 * f_gc - 3.95) * 1e-5 * log_salt
                  + 9.40e-6 * log_salt ** 2)
        if inv_tm > 0:
            tm_kelvin = 1.0 / inv_tm

    return round(tm_kelvin - 273.15, 1)


def find_hairpins(
    seq: str,
    stem_min: int = 6,
    loop_min: int = 3,
    loop_max: int = 8,
) -> list[dict]:
    """Detect potential hairpin (stem-loop) structures in *seq*.

    Uses nearest-neighbor thermodynamics (SantaLucia 1998) for Tm
    estimation instead of the Wallace rule.  Perfect and single-mismatch
    stems are reported; G-T wobble pairs are tolerated.

    Returns
    -------
    list of dict
        Keys: position, stem_length, loop_length, stem_seq, stem_gc,
        tm_estimate (NN-based), mismatches.
    """
    seq = seq.upper()
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    # G-T wobble is a common non-Watson-Crick pair in stems
    wobble = {("G", "T"), ("T", "G")}
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

                # Count mismatches, allowing G-T wobble
                mismatches = 0
                for a, b in zip(left, reversed(right)):
                    if comp.get(a) == b:
                        continue  # Watson-Crick
                    elif (a, b) in wobble:
                        mismatches += 1  # wobble counts as half-mismatch
                    else:
                        mismatches = 999  # real mismatch
                        break

                # Allow 0 mismatches for short stems, up to 1 wobble for stems >= 6
                max_mm = 1 if stem_len >= 6 else 0
                if mismatches > max_mm:
                    continue

                gc = gc_content(left)
                tm_est = _hairpin_nn_tm(left, loop_len)
                hits.append({
                    "position": left_start,
                    "stem_length": stem_len,
                    "loop_length": loop_len,
                    "stem_seq": left,
                    "stem_gc": gc,
                    "tm_estimate": tm_est,
                    "mismatches": mismatches,
                })

    # Deduplicate: keep strongest hairpin per 5-bp window (by Tm, then stem length)
    seen: set[int] = set()
    unique: list[dict] = []
    for hp in sorted(hits, key=lambda x: (-x["tm_estimate"], -x["stem_length"])):
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
    comp = {
        "A": "T", "T": "A", "G": "C", "C": "G",
        "N": "N", "R": "Y", "Y": "R", "M": "K", "K": "M",
        "S": "S", "W": "W", "B": "V", "V": "B", "D": "H", "H": "D",
    }
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


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
    rev = rev_bind.upper()

    fwd_match = fwd
    fwd_pos = template.find(fwd_match)
    if fwd_pos == -1 and len(fwd) > 15:
        fwd_match = fwd[-15:]
        fwd_pos = template.find(fwd_match)

    rev_rc = _reverse_complement(rev)
    rev_match = rev_rc
    rev_pos = template.find(rev_match)
    if rev_pos == -1 and len(rev) > 15:
        rev_match = rev_rc[:15]
        rev_pos = template.find(rev_match)

    if fwd_pos == -1 or rev_pos == -1:
        return None, fwd_pos, rev_pos

    if fwd_pos <= rev_pos:
        a_len = int(len(fwd_match))
        template_part = template[fwd_pos + a_len:rev_pos]
        amplicon = fwd + template_part + rev_rc
    else:
        a_len = int(len(fwd_match))
        template_part = template[fwd_pos + a_len:] + template[:rev_pos]
        amplicon = fwd + template_part + rev_rc

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
    1. Primer self-hairpins (stem >= 5 bp, NN Tm >= 40 degC).
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
                "Install with: pip install polymerase-tm[bio]"
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
