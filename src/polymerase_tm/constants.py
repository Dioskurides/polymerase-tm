"""Thermodynamic and polymerase constants for polymerase-tm.

Contains nearest-neighbor parameters (SantaLucia 1998),
NEB buffer salt concentrations, and polymerase product definitions.
"""

from __future__ import annotations

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
        "description": "Taq DNA Polymerase",
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
    "onetaq": {
        "buffer": "onetaq_std", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "OneTaq DNA Polymerase (Standard Buffer)",
    },
    "onetaq_gc": {
        "buffer": "onetaq_gc", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "OneTaq DNA Polymerase (GC Buffer)",
    },
    "onetaq_hot_start": {
        "buffer": "onetaq_std", "conc": 200, "ta_rule": "offset",
        "ta_offset": -5, "ta_cap": 68, "min_len": 8,
        "description": "OneTaq Hot Start DNA Polymerase",
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
