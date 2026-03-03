"""
polymerase-tm — NEB Tm Calculator, Python Implementation

Reproduces the melting temperature (Tm) and annealing temperature (Ta)
calculations of the NEB Tm Calculator (https://tmcalculator.neb.com/).

Also includes the NEB Base Changer primer design algorithm for
site-directed mutagenesis (SDM) using the Q5 SDM Kit.

Algorithm
---------
- Tm:  SantaLucia (1998) nearest-neighbor model.
- Salt correction:  Owczarzy et al. (2004, 2008) monovalent-ion method.
- Bivariate salt:  Owczarzy et al. (2008) Na⁺ + Mg²⁺ (for SDM).
- Ta:  Polymerase-specific offset + cap (see POLYMERASES table).
- DMSO correction:  -0.6 degC per 1 % DMSO (v/v).

Buffer salt concentrations and Ta rules were extracted from the NEB
Tm Calculator front-end source (main-90540d67e2.js, version 1.16.10).
NEBaseChanger algorithms from index.dfcf8292.js, version 2.7.2.

References
----------
1. SantaLucia J Jr. (1998) PNAS 95:1460-5.
2. Owczarzy R et al. (2004) Biochemistry 43:3537-54.
3. Owczarzy R et al. (2008) Biochemistry 47:5336-53.
4. Allawi HT, SantaLucia J Jr. (1997-1999) Mismatch parameters.
"""

from __future__ import annotations

# --- Constants ---
from .constants import (
    NN_PARAMS, TERMINAL, BUFFERS, POLYMERASES,
    MISMATCH_NN, CODON_TABLE, AA_TO_CODONS, ECOLI_CODON_USAGE,
    GENETIC_CODES, get_codon_table, get_aa_to_codons,
)

# --- Core Tm/Ta ---
from .core import (
    calc_nn_raw,
    owczarzy_correction,
    owczarzy_bivariate,
    calc_sdm_tm,
    tm,
    ta,
    list_polymerases,
    list_buffers,
)

# --- DMSO / Secondary-Structure Analysis ---
from .dmso import (
    gc_content,
    gc_windows,
    find_gc_hotspots,
    find_hairpins,
    primer_hairpin,
    analyze_amplicon,
    dmso_recommendation,
    print_dmso_report,
)

# --- Batch Processing / Automation ---
from .batch import (
    reverse_complement,
    batch_tm,
    optimal_binding_length,
    check_pair,
    pcr_protocol,
    from_csv,
    to_csv,
    additive_recommendation,
)

# --- Primer Analysis ---
from .analysis import (
    RESTRICTION_ENZYMES,
    primer_dimer,
    gibson_overlaps,
    restriction_scan,
    primer_quality,
)

# --- Visualization ---
from .gel import plot_virtual_gel

# --- Site-Directed Mutagenesis (NEB Base Changer) ---
from .mutagenesis import (
    BaseChanger,
    SDMPrimer,
    MutagenesisResult,
    select_codon,
    parse_aa_mutation,
)

__version__ = "2.0.2"

__all__ = [
    # Core
    "tm",
    "ta",
    "list_polymerases",
    "list_buffers",
    "calc_nn_raw",
    "owczarzy_correction",
    "owczarzy_bivariate",
    "calc_sdm_tm",
    # Automation
    "reverse_complement",
    "batch_tm",
    "optimal_binding_length",
    "check_pair",
    "pcr_protocol",
    "from_csv",
    "to_csv",
    "additive_recommendation",
    # Analysis
    "primer_dimer",
    "gibson_overlaps",
    "restriction_scan",
    "primer_quality",
    # DMSO
    "dmso_recommendation",
    "print_dmso_report",
    "gc_content",
    "gc_windows",
    "find_gc_hotspots",
    "find_hairpins",
    "primer_hairpin",
    "analyze_amplicon",
    # Visualization
    "plot_virtual_gel",
    # Mutagenesis
    "BaseChanger",
    "SDMPrimer",
    "MutagenesisResult",
    "select_codon",
    "parse_aa_mutation",
    # Data
    "NN_PARAMS",
    "TERMINAL",
    "BUFFERS",
    "POLYMERASES",
    "RESTRICTION_ENZYMES",
    "MISMATCH_NN",
    "CODON_TABLE",
    "AA_TO_CODONS",
    "ECOLI_CODON_USAGE",
    "GENETIC_CODES",
    "get_codon_table",
    "get_aa_to_codons",
]

import os

_install_flag = os.path.expanduser("~/.polymerase_tm_installed_v2")
if not os.path.exists(_install_flag):
    try:
        print("Vocatus et invocatus deus aderit")
        with open(_install_flag, "w") as f:
            f.write("installed\n")
    except Exception:
        pass
