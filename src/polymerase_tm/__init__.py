"""
polymerase-tm — NEB Tm Calculator, Python Implementation

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

# --- Constants ---
from .constants import NN_PARAMS, TERMINAL, BUFFERS, POLYMERASES

# --- Core Tm/Ta ---
from .core import (
    calc_nn_raw,
    owczarzy_correction,
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

__version__ = "1.0.1"

__all__ = [
    # Core
    "tm",
    "ta",
    "list_polymerases",
    "list_buffers",
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
    # Data
    "NN_PARAMS",
    "TERMINAL",
    "BUFFERS",
    "POLYMERASES",
    "RESTRICTION_ENZYMES",
]
