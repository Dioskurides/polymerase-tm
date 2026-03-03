"""Batch processing, automation, and PCR protocol generation.

Contains reverse_complement, batch_tm, optimal_binding_length,
check_pair, pcr_protocol, from_csv, to_csv, and _additive_recommendation.
"""

from __future__ import annotations

from typing import Optional

from .core import tm, ta
from .dmso import gc_content, primer_hairpin, _reverse_complement, analyze_amplicon


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


def additive_recommendation(
    fwd: str,
    rev: str,
    polymerase: str = "q5",
    amplicon_gc: Optional[float] = None,
) -> dict:
    """Analyse whether DMSO or GC enhancer should be added.

    Logic
    -----
    - Q5 family: recommend **Q5 High GC Enhancer** (1x) for GC-rich targets.
      DMSO is NOT recommended with Q5 (NEB advises against it).
    - Other polymerases: recommend **DMSO** (2-5 % v/v).

    Triggers (any one is sufficient):
    - Primer GC > 70 %
    - Average primer GC > 65 %
    - Strong self-hairpin in primer (stem >= 5 bp, Tm >= 40 degC)
    - Amplicon GC > 65 % (if provided)

    Parameters
    ----------
    fwd, rev : str
        Primer binding sequences.
    polymerase : str
        Polymerase key (default ``"q5"``).
    amplicon_gc : float, optional
        Amplicon GC content as percentage (0-100).

    Returns
    -------
    dict
        Keys: ``recommended`` (bool), ``additive`` (str or None),
        ``concentration`` (str or None), ``reasons`` (list[str]).
    """
    fwd = fwd.strip().upper()
    rev = rev.strip().upper()
    gc_f = gc_content(fwd) * 100
    gc_r = gc_content(rev) * 100
    avg_gc = (gc_f + gc_r) / 2

    is_q5 = "q5" in polymerase.lower()
    reasons = []

    # Check primer GC
    if gc_f > 70:
        reasons.append(f"Forward primer GC is {gc_f:.0f}% (> 70%)")
    if gc_r > 70:
        reasons.append(f"Reverse primer GC is {gc_r:.0f}% (> 70%)")
    if avg_gc > 65 and gc_f <= 70 and gc_r <= 70:
        reasons.append(f"Average primer GC is {avg_gc:.0f}% (> 65%)")

    # Check primer hairpins
    for label, seq in [("fwd", fwd), ("rev", rev)]:
        hps = primer_hairpin(seq)
        strong = [h for h in hps if h["stem_length"] >= 5 and h["tm_estimate"] >= 40]
        if strong:
            reasons.append(
                f"Primer {label} has {len(strong)} strong hairpin(s) "
                f"(stem >= 5 bp, Tm >= 40 degC)"
            )

    # Check amplicon GC if provided
    if amplicon_gc is not None and amplicon_gc > 65:
        reasons.append(f"Amplicon GC is {amplicon_gc:.0f}% (> 65%)")

    if not reasons:
        return {
            "recommended": False,
            "additive": None,
            "concentration": None,
            "reasons": ["No GC-related issues detected"],
        }

    if is_q5:
        return {
            "recommended": True,
            "additive": "Q5 High GC Enhancer",
            "concentration": "1x",
            "reasons": reasons,
            "note": "NEB does not recommend DMSO with Q5. Use GC Enhancer instead.",
        }
    else:
        # Recommend DMSO percentage based on severity
        pct = 3 if len(reasons) <= 2 else 5
        return {
            "recommended": True,
            "additive": "DMSO",
            "concentration": f"{pct}%",
            "reasons": reasons,
        }


# Backward-compatible alias (was private in v0.9.x)
_additive_recommendation = additive_recommendation


def check_pair(
    fwd: str,
    rev: str,
    polymerase: str = "q5",
    dmso_pct: float = 0,
) -> dict:
    """Comprehensive primer pair compatibility check.

    Automatically includes DMSO / GC Enhancer recommendation.

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
        ``gc_fwd``, ``gc_rev``, ``compatible`` (bool), ``warnings`` (list),
        ``additive`` (dict with DMSO/GC enhancer recommendation).

    Examples
    --------
    >>> result = check_pair("ATGTCCCTGCTCTTCTCTCGATGCAA", "GTGCCTCCGAGCCAGCACC")
    >>> if result["additive"]["recommended"]:
    ...     print(f"Use {result['additive']['additive']} ({result['additive']['concentration']})")
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

    additive = additive_recommendation(fwd, rev, polymerase=polymerase)

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
        "additive": additive,
    }


def pcr_protocol(
    fwd: str,
    rev: str,
    polymerase: str = "q5",
    amplicon_length: Optional[int] = None,
    dmso_pct: float = 0,
    num_cycles: int = 30,
    template: Optional[str] = None,
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
    template : str, optional
        Template sequence. If provided, amplicon_length is derived
        automatically.

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

    # Automatically determine amplicon length if template is provided
    expected_amplicon_seq = None
    if template:
        amp_seq, _, _ = analyze_amplicon(template, fwd, rev)
        if amp_seq:
            expected_amplicon_seq = amp_seq
            if amplicon_length is None:
                amplicon_length = len(amp_seq)

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

    additive = additive_recommendation(fwd, rev, polymerase=polymerase)

    return {
        "polymerase": polymerase,
        "polymerase_family": family,
        "fwd_tm": t1,
        "rev_tm": t2,
        "ta": result_ta,
        "dmso_pct": dmso_pct,
        "amplicon_length": amplicon_length,
        "amplicon_sequence": expected_amplicon_seq,
        "extension_rate_s_per_kb": params["rate_s_per_kb"],
        "num_cycles": num_cycles,
        "cycling": cycling,
        "total_time_min": round(total_seconds / 60, 1),
        "additive": additive,
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
