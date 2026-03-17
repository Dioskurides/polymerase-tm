# -*- coding: utf-8 -*-
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
    "standard_taq":   51,
    "thermopol":      40,
    "hemo_klentaq":   70,
    "crimson_taq":    51,
    "longamp":       100,
    "multiplex":      90,
    "phusion_hf":    222,
    "phusion_gc":    222,
    "phusionflex_hf": 222,
    "phusionflex_gc": 222,
    "onetaq_std":     51,
    "onetaq_gc":      51,
    "q5":            150,
    "q5mm":          150,
    "q5u":           180,
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
        "ta_cap": 72, "min_len": 0, "tm_method": "schildkraut",
        "description": "Phusion High-Fidelity DNA Polymerase (HF Buffer)",
    },
    "phusion_gc": {
        "buffer": "phusion_gc", "conc": 500, "ta_rule": "phusion",
        "ta_cap": 72, "min_len": 0, "tm_method": "schildkraut",
        "description": "Phusion High-Fidelity DNA Polymerase (GC Buffer)",
    },
    "phusion_flex_hf": {
        "buffer": "phusionflex_hf", "conc": 500, "ta_rule": "phusion",
        "ta_cap": 72, "min_len": 0, "tm_method": "schildkraut",
        "description": "Phusion Hot Start Flex DNA Polymerase (HF Buffer)",
    },
    "phusion_flex_gc": {
        "buffer": "phusionflex_gc", "conc": 500, "ta_rule": "phusion",
        "ta_cap": 72, "min_len": 0, "tm_method": "schildkraut",
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


# =====================================================================
# Owczarzy et al. (2008) Bivariate Salt Correction Coefficients
# Used by NEBaseChanger for Mg²⁺-aware Tm correction.
# Biochemistry 47:5336-53.
# =====================================================================

OWCZARZY_BIVARIATE: dict[str, float] = {
    "a": -3.22e-5,   # coefficient a
    "b": 6.39e-5,    # coefficient b
    "c": 1.75e-5,    # coefficient c
    "d": -2.63e-4,   # coefficient d
    "e": 4.29e-5,    # coefficient e (same as mono)
    "f": 3.95e-5,    # coefficient f
    "g": 9.40e-6,    # coefficient g
    # Divalent-only coefficients for sc_ow
    "l": -0.911e-5,
    "r": 6.26e-5,
    "ff": -48.2e-5,
    "u": 52.5e-5,
    "n": 3.92e-5,
    "aa": 1.42e-5,
    "h": 8.31e-5,
}

# NEBaseChanger Q5 SDM default buffer chemistry
Q5_SDM_MONO_MM = 50     # monovalent salt (mM)
Q5_SDM_DIVALENT_MM = 2  # Mg²⁺ (mM)
Q5_SDM_PRIMER_NM = 500  # primer concentration (nM)
Q5_SDM_MIN_TM = 55.0    # minimum primer Tm (°C)
Q5_SDM_MAX_TM = 72.0    # maximum primer Tm (°C)
Q5_SDM_MIN_LEN = 15     # minimum primer length (nt)
Q5_SDM_MIN_ANNEAL = 10  # minimum annealing bases flanking mutation


# =====================================================================
# Mismatch Nearest-Neighbor Parameters
# Allawi & SantaLucia (1997-1999) — internal single mismatches
# Key format: "XY/WZ" → X-W mismatch flanked by Y and Z context
# Values: (dH kcal/mol, dS cal/(mol*K))
#
# References:
#  - Allawi & SantaLucia 1997 Biochemistry 36:10581 (G-T)
#  - Allawi & SantaLucia 1998 Biochemistry 37:9435 (A-C)
#  - Allawi & SantaLucia 1998 Biochemistry 37:2170 (G-A)
#  - Allawi & SantaLucia 1998 Nucleic Acids Res 26:2694 (C-T)
#  - Peyret et al. 1999 Biochemistry 38:3468 (A-A, C-C, G-G, T-T)
# =====================================================================

MISMATCH_NN: dict[str, tuple[float, float]] = {
    # A-C mismatches
    "AA/CA": (2.3, 4.6),   "AC/AG": (-0.7, -3.8),
    "AG/AC": (-7.2, -19.6),"AT/AA": (-2.3, -8.0),
    "CA/GA": (-0.9, -4.2), "CC/GG": (-0.6, -3.8),
    "CG/GA": (-4.0, -13.2),"CT/GC": (-1.5, -6.1),
    "GA/CA": (5.3, 14.6),  "GC/CG": (-0.6, -1.0),
    "GG/CC": (0.5, 3.2),   "GT/CA": (-2.2, -8.4),
    "TA/AA": (3.4, 8.0),   "TC/AG": (-1.2, -3.7),
    "TG/AC": (1.0, 0.7),   "TT/AC": (-0.2, -1.5),
    # G-T mismatches
    "AA/TG": (-3.1, -9.5), "AC/TT": (-0.6, -0.6),
    "AG/TC": (-4.3, -12.3),"AT/TA": (0.7, 0.7),
    "CA/GG": (-0.7, -2.3), "CC/GT": (-3.1, -10.6),
    "CG/GC": (-4.4, -12.5),"CT/GA": (2.9, 8.3),
    "GA/TG": (-1.6, -5.6), "GC/TT": (-4.1, -13.0),
    "GG/TC": (-6.0, -15.8),"GT/TA": (-1.4, -4.2),
    "TA/AG": (0.7, 0.7),   "TC/AT": (-0.1, -1.7),
    "TG/AC": (1.0, 0.7),   "TT/AA": (1.0, 0.7),
    # G-A mismatches
    "AA/AG": (-2.9, -9.8), "AC/AT": (4.7, 12.9),
    "AG/AA": (-0.6, -2.3), "AT/AC": (-0.7, -2.3),
    "CA/GG": (-0.7, -2.3), "CC/GA": (5.2, 14.2),
    "CG/GT": (-0.6, -1.0), "CT/GG": (3.6, 8.9),
    "GA/AG": (-2.9, -9.8), "GC/AA": (-6.0, -15.8),
    "GG/CA": (-3.1, -10.6),"GT/GA": (-1.4, -4.2),
    "TA/TG": (-1.0, -2.3), "TC/TA": (-4.1, -13.0),
    "TG/TA": (-1.6, -5.6), "TT/TG": (-0.6, -0.6),
    # C-T mismatches
    "AA/TC": (4.7, 12.9),  "AC/TA": (0.7, 0.7),
    "AG/TG": (-3.1, -9.5), "AT/TT": (-2.7, -10.8),
    "CA/GC": (6.1, 16.4),  "CC/GC": (-0.8, -4.5),
    "CG/GG": (-4.4, -12.5),"CT/GT": (-0.1, -1.7),
    "GA/TC": (-2.9, -9.8), "GC/TA": (-6.0, -15.8),
    "GG/TG": (-6.0, -15.8),"GT/TC": (-4.1, -13.0),
    "TA/AC": (3.4, 8.0),   "TC/AA": (-0.1, -1.7),
    "TG/AG": (-4.3, -12.3),"TT/AT": (0.2, -1.5),
    # A-A mismatches
    "AA/AA": (1.2, 1.7),   "AC/AT": (4.7, 12.9),
    "AG/AT": (-0.6, -2.3), "AT/AA": (-2.3, -8.0),
    "CA/GA": (-0.9, -4.2), "CC/GT": (-3.1, -10.6),
    "CG/GA": (-4.0, -13.2),"CT/GA": (2.9, 8.3),
    "GA/AA": (-2.9, -9.8), "GC/AT": (-6.0, -15.8),
    "GG/AT": (-3.1, -10.6),"GT/AA": (-2.2, -8.4),
    "TA/TA": (-1.3, -5.3), "TC/TT": (-0.1, -1.7),
    "TG/TA": (-1.6, -5.6), "TT/TA": (-0.6, -0.6),
    # G-G mismatches
    "AG/GC": (-4.0, -13.2),"CG/GC": (-4.4, -12.5),
    "GG/GC": (3.3, 10.4),  "TG/GA": (-1.6, -5.6),
    "AG/GA": (-0.6, -2.3), "CG/GA": (-4.0, -13.2),
    "GG/GA": (-6.0, -15.8),"TG/GC": (-4.4, -12.5),
    # C-C mismatches
    "AC/CG": (0.6, -0.6),  "CC/CG": (-0.8, -4.5),
    "GC/CA": (-0.6, -1.0), "TC/CG": (-1.5, -6.1),
    "AC/CA": (-0.7, -3.8), "CC/CA": (-0.6, -3.8),
    "GC/CG": (-0.6, -1.0), "TC/CA": (-1.2, -3.7),
    # T-T mismatches
    "AT/TT": (-2.7, -10.8),"CT/TG": (-0.1, -1.7),
    "GT/TA": (-1.4, -4.2), "TT/TA": (-0.6, -0.6),
    "AT/TA": (0.7, 0.7),   "CT/TA": (2.9, 8.3),
    "GT/TG": (-1.4, -4.2), "TT/TT": (-0.2, -1.5),
}


# =====================================================================
# Internal Loop and Bulge Entropies
# SantaLucia & Hicks 2004, Table 4
# (dS in cal/(mol*K) — dH assumed 0 for loops/bulges)
# =====================================================================

LOOP_PENALTIES: dict[int, float] = {
    3: -10.3, 4: -11.0, 5: -12.1, 6: -13.2, 7: -13.7,
    8: -14.1, 9: -15.0, 10: -15.8, 12: -17.1, 14: -18.1,
    16: -18.9, 18: -19.5, 20: -20.1, 25: -21.3, 30: -22.2,
}

BULGE_PENALTIES: dict[int, float] = {
    1: -12.9, 2: -7.0, 3: -10.3, 4: -11.0, 5: -12.1,
    6: -13.2, 7: -13.7, 8: -14.1, 9: -15.0, 10: -15.8,
    12: -17.1, 14: -18.1, 16: -18.9, 18: -19.5, 20: -20.1,
    25: -21.3, 30: -22.2,
}


# =====================================================================
# Terminal Mismatch Parameters
# Bommarito et al. 2000, Nucleic Acids Res 28:1929
# Key: "XY/WZ" where X-W is terminal mismatch, Y/Z are flanking
# =====================================================================

TERMINAL_MISMATCH: dict[str, tuple[float, float]] = {
    "AA/TA": (-3.1, -7.8), "AC/TG": (-4.4, -12.5),
    "AG/TC": (-4.0, -13.2),"CA/GA": (-0.9, -4.2),
    "CC/GG": (-4.4, -12.5),"CG/GC": (-4.0, -13.2),
    "GA/CA": (0.1, -2.3),  "GC/CG": (-0.6, -1.0),
    "GG/CC": (0.5, 3.2),   "TA/AA": (-0.1, -1.7),
    "TC/AG": (-1.5, -6.1), "TG/AC": (1.0, 0.7),
}


# =====================================================================
# Dangling End Parameters
# Bommarito et al. 2000, Nucleic Acids Res 28:1929
# Key format: "-X/YZ" (5' dangle) or "XY/Z-" (3' dangle)
# =====================================================================

DANGLING_ENDS: dict[str, tuple[float, float]] = {
    # 5' dangling ends ("-X/YZ")
    "-A/AT": (-0.7, -0.8), "-A/CT": (4.4, 14.9),
    "-A/GT": (-1.6, -3.6), "-A/TT": (2.9, 10.4),
    "-C/AG": (-2.1, -3.9), "-C/CG": (-0.2, -0.1),
    "-C/GG": (-3.9, -11.2),"-C/TG": (-0.5, -0.8),
    "-G/AC": (-5.9, -16.5),"-G/CC": (-2.6, -7.4),
    "-G/GC": (-3.2, -10.4),"-G/TC": (-3.9, -10.6),
    "-T/AA": (0.2, 2.4),   "-T/CA": (-6.3, -17.1),
    "-T/GA": (-3.2, -10.4),"-T/TA": (-2.5, -6.3),
    # 3' dangling ends ("XY/Z-")
    "AA/T-": (-0.5, -1.1), "AC/G-": (-3.7, -10.0),
    "AG/C-": (-2.9, -7.6), "AT/A-": (-3.2, -8.9),
    "CA/T-": (-5.9, -16.5),"CC/G-": (-2.6, -7.4),
    "CG/C-": (-3.9, -11.2),"CT/A-": (-4.3, -10.6),
    "GA/T-": (-2.1, -3.9), "GC/G-": (-0.2, -0.1),
    "GG/C-": (-3.9, -10.4),"GT/A-": (-3.2, -10.4),
    "TA/T-": (0.2, 2.4),   "TC/G-": (-6.3, -17.1),
    "TG/C-": (-4.4, -12.5),"TT/A-": (-2.5, -6.3),
}


# =====================================================================
# Standard Genetic Code and Codon Tables
# NCBI Translation Table IDs
# =====================================================================

# Amino acid three-letter to one-letter
AA3_TO_AA1: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "STOP": "*",
}

AA1_TO_AA3: dict[str, str] = {v: k for k, v in AA3_TO_AA1.items()}

# Standard genetic code (NCBI Table 1)
CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Reverse lookup: amino acid → list of codons
AA_TO_CODONS: dict[str, list[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)


# NCBI Genetic Codes -- all 12 codes supported by NEBaseChanger
# Each entry stores only the DIFFERENCES from the standard code (Table 1).
GENETIC_CODES: dict[int, dict] = {
    1: {"name": "Standard", "diff": {}},
    2: {"name": "Vertebrate Mitochondrial", "diff": {
        "AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W",
    }},
    3: {"name": "Yeast Mitochondrial", "diff": {
        "ATA": "M", "CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T",
        "TGA": "W",
    }},
    4: {"name": "Mold/Protozoan/Coelenterate Mitochondrial", "diff": {
        "TGA": "W",
    }},
    5: {"name": "Invertebrate Mitochondrial", "diff": {
        "AGA": "S", "AGG": "S", "ATA": "M", "TGA": "W",
    }},
    6: {"name": "Ciliate/Dasycladacean/Hexamita Nuclear", "diff": {
        "TAA": "Q", "TAG": "Q",
    }},
    9: {"name": "Echinoderm/Flatworm Mitochondrial", "diff": {
        "AAA": "N", "AGA": "S", "AGG": "S", "TGA": "W",
    }},
    10: {"name": "Euplotid Nuclear", "diff": {
        "TGA": "C",
    }},
    11: {"name": "Bacterial/Archaeal/Plant Plastid", "diff": {}},
    12: {"name": "Alternative Yeast Nuclear", "diff": {
        "CTG": "S",
    }},
    13: {"name": "Ascidian Mitochondrial", "diff": {
        "AGA": "G", "AGG": "G", "ATA": "M", "TGA": "W",
    }},
    14: {"name": "Alternative Flatworm Mitochondrial", "diff": {
        "AAA": "N", "AGA": "S", "AGG": "S", "TAA": "Y", "TGA": "W",
    }},
}


def get_codon_table(code_id: int = 1) -> dict[str, str]:
    """Return the full codon->amino acid table for a given NCBI code ID."""
    if code_id not in GENETIC_CODES:
        raise ValueError(
            f"Unknown genetic code ID {code_id}. "
            f"Available: {sorted(GENETIC_CODES.keys())}"
        )
    table = dict(CODON_TABLE)  # copy standard
    table.update(GENETIC_CODES[code_id]["diff"])
    return table


def get_aa_to_codons(code_id: int = 1) -> dict[str, list[str]]:
    """Return amino acid -> list of codons for a given genetic code."""
    table = get_codon_table(code_id)
    result: dict[str, list[str]] = {}
    for codon, aa in table.items():
        result.setdefault(aa, []).append(codon)
    return result


# E. coli codon usage frequencies (per thousand)
# From the Kazusa Codon Usage Database
ECOLI_CODON_USAGE: dict[str, float] = {
    "TTT": 22.0, "TTC": 16.0, "TTA": 13.7, "TTG": 13.4,
    "CTT": 11.5, "CTC": 10.4, "CTA": 3.9,  "CTG": 52.0,
    "ATT": 30.1, "ATC": 24.2, "ATA": 4.6,  "ATG": 27.8,
    "GTT": 18.2, "GTC": 15.1, "GTA": 10.9, "GTG": 25.9,
    "TCT": 8.5,  "TCC": 8.5,  "TCA": 7.2,  "TCG": 8.8,
    "CCT": 7.1,  "CCC": 5.5,  "CCA": 8.4,  "CCG": 22.9,
    "ACT": 9.0,  "ACC": 22.9, "ACA": 7.1,  "ACG": 14.4,
    "GCT": 15.5, "GCC": 25.3, "GCA": 20.2, "GCG": 33.1,
    "TAT": 16.0, "TAC": 12.0, "TAA": 2.0,  "TAG": 0.3,
    "CAT": 12.9, "CAC": 9.5,  "CAA": 15.3, "CAG": 29.1,
    "AAT": 18.2, "AAC": 21.5, "AAA": 33.5, "AAG": 10.3,
    "GAT": 32.4, "GAC": 19.0, "GAA": 39.4, "GAG": 17.8,
    "TGT": 5.2,  "TGC": 6.3,  "TGA": 1.0,  "TGG": 15.2,
    "CGT": 20.7, "CGC": 21.4, "CGA": 3.6,  "CGG": 5.6,
    "AGT": 9.0,  "AGC": 16.0, "AGA": 2.1,  "AGG": 1.2,
    "GGT": 24.6, "GGC": 28.9, "GGA": 8.0,  "GGG": 11.0,
}


# =====================================================================
# Restriction Enzyme Recognition Sites (5'->3')
# =====================================================================

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
