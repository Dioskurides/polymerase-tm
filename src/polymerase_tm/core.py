"""Core Tm/Ta calculation functions.

Implements the SantaLucia (1998) nearest-neighbor model with
Owczarzy (2004) salt correction.
"""

from __future__ import annotations

import math
from typing import Optional

from .constants import NN_PARAMS, TERMINAL, R, BUFFERS, POLYMERASES


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
    if not seq:
        raise ValueError("Primer sequence must not be empty.")
    invalid = set(seq) - {"A", "T", "G", "C"}
    if invalid:
        raise ValueError(
            f"Invalid character(s) in primer sequence: {', '.join(sorted(invalid))}. "
            "Only A, T, G, C are allowed (no overhangs, no IUPAC ambiguity codes)."
        )
    dH = 0.0
    dS = 0.0
    for i in range(len(seq) - 1):
        dinuc = seq[i:i + 2]
        h, s = NN_PARAMS[dinuc]
        dH += h
        dS += s

    # Terminal corrections
    dH += TERMINAL[seq[0]][0] + TERMINAL[seq[-1]][0]
    dS += TERMINAL[seq[0]][1] + TERMINAL[seq[-1]][1]

    # Convert primer concentration: nM -> M
    conc_M = primer_conc_nM * 1e-9
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
    if salt_mM <= 0:
        raise ValueError(f"Salt concentration must be positive, got {salt_mM} mM.")
    salt_M = salt_mM / 1000.0
    ln_s = math.log(salt_M)
    corr = (4.29e-5 * fgc - 3.95e-5) * ln_s + 9.40e-6 * ln_s ** 2
    return 1.0 / (1.0 / (tm_raw + 273.15) + corr) - 273.15


def _resolve_salt(
    polymerase: str,
    buffer: Optional[str] = None,
    salt_mM: Optional[int] = None,
) -> int:
    """Resolve the effective salt concentration from the given parameters."""
    if salt_mM is not None:
        return salt_mM
    if buffer is not None:
        key = buffer.lower().replace(" ", "_").replace("-", "_")
        if key not in BUFFERS:
            avail = ", ".join(sorted(BUFFERS.keys()))
            raise ValueError(
                f"Unknown buffer '{buffer}'. Available buffers: {avail}"
            )
        return BUFFERS[key]
    poly = POLYMERASES[polymerase]
    return BUFFERS[poly["buffer"]]


def list_buffers() -> list[dict]:
    """Return all available NEB buffer names and their salt concentrations.

    Returns
    -------
    list[dict]
        Each dict has keys: ``name``, ``salt_mM``.

    Examples
    --------
    >>> from polymerase_tm import list_buffers
    >>> for b in list_buffers():
    ...     print(f"{b['name']:20s} {b['salt_mM']:>4d} mM")
    """
    return [{"name": k, "salt_mM": v} for k, v in sorted(BUFFERS.items())]


def tm(
    seq: str,
    polymerase: str = "q5",
    buffer: Optional[str] = None,
    salt_mM: Optional[int] = None,
) -> int:
    """Calculate Tm for a primer sequence.

    Parameters
    ----------
    seq : str
        Primer binding sequence (A/T/G/C only, no overhangs).
    polymerase : str
        Key into ``POLYMERASES``. Default ``"q5"``.
    buffer : str, optional
        NEB buffer name to override the polymerase default.
        Use ``list_buffers()`` to see all 17 available buffers.
    salt_mM : int, optional
        Direct salt concentration (mM) override. Takes priority
        over both ``polymerase`` and ``buffer``.

    Returns
    -------
    int
        Melting temperature rounded to the nearest integer (degC).

    Examples
    --------
    >>> tm("ATCGATCGATCG")                    # Q5 default
    >>> tm("ATCGATCGATCG", buffer="thermopol") # Thermopol buffer
    >>> tm("ATCGATCGATCG", salt_mM=50)         # custom 50 mM
    """
    polymerase = polymerase.lower().replace(" ", "_").replace("-", "_")
    if polymerase not in POLYMERASES:
        avail = ", ".join(sorted(POLYMERASES.keys()))
        raise ValueError(
            f"Unknown polymerase '{polymerase}'. Available: {avail}"
        )
    poly = POLYMERASES[polymerase]
    salt = _resolve_salt(polymerase, buffer, salt_mM)
    conc = poly["conc"]
    raw = calc_nn_raw(seq, conc)
    return round(owczarzy_correction(raw, seq, salt))


def ta(
    seq1: str,
    seq2: str,
    polymerase: str = "q5",
    dmso_pct: float = 0,
    buffer: Optional[str] = None,
    salt_mM: Optional[int] = None,
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
    buffer : str, optional
        NEB buffer name to override the polymerase default.
    salt_mM : int, optional
        Direct salt concentration (mM) override.

    Returns
    -------
    (Ta, Tm1, Tm2) : tuple[int, int, int]
        Recommended annealing temperature, Tm of seq1, Tm of seq2.
    """
    t1 = tm(seq1, polymerase, buffer=buffer, salt_mM=salt_mM)
    t2 = tm(seq2, polymerase, buffer=buffer, salt_mM=salt_mM)
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
        Each dict has keys: key, description, buffer_salt_mM, primer_conc_nM, ta_rule.
    """
    out = []
    for key, poly in POLYMERASES.items():
        out.append({
            "key": key,
            "description": poly["description"],
            "buffer_salt_mM": BUFFERS[poly["buffer"]],
            "primer_conc_nM": poly["conc"],
            "ta_rule": poly.get("ta_rule", ""),
        })
    return out
