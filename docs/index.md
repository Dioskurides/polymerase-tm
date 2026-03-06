# polymerase-tm

Exact Python reproduction of the [NEB Tm Calculator](https://tmcalculator.neb.com/) and [NEB Base Changer](https://nebasechanger.neb.com/) for PCR primer melting temperature (Tm), annealing temperature (Ta), and site-directed mutagenesis (SDM) primer design.

## Features

- **Exact NEB Tm Calculator reproduction** -- algorithm recreated perfectly from the NEB Tm Calculator; verified against the official tool with 0 degC deviation across all tested sequences.
- **Core functions** -- `tm()`, `ta()`, `batch_tm()`, `check_pair()`, and comprehensive primer analysis.
- **22 NEB polymerase products** with their specific buffer salt concentrations and Ta rules (Q5, Phusion, Taq, OneTaq, LongAmp, Vent, Deep Vent, and more).
- **Automatic additive recommendation** -- suggests Q5 High GC Enhancer or DMSO based on primer GC, hairpins, and amplicon analysis.
- **Batch processing** -- process hundreds of primer pairs from CSV files with full Tm/Ta/compatibility analysis.
- **PCR protocol generator** -- generates complete cycling protocols with polymerase-specific temperatures and extension times. Automatically generates **touchdown protocols** when primer Tm difference exceeds 3 degC.
- **Smart primer design** -- find the optimal binding length for a target Tm.
- **Primer dimer detection** -- checks 3'-end complementarity and self-dimer risk.
- **Gibson Assembly overlap design** -- generates full primers with overhangs for Gibson/HiFi Assembly.
- **Restriction site scanning** -- scans primers for ~120 NEB restriction enzyme sites.
- **Site-directed mutagenesis** -- reimplementation of the [NEB Base Changer v2.7.2](https://nebasechanger.neb.com/) primer design algorithm.
- **DMSO analysis** -- analyses primer hairpins, amplicon GC content, GC-rich hotspots, and template secondary structures.
- **Virtual gel visualization** -- simulated agarose gel with realistic Ferguson-plot migration physics.

## Installation

```bash
# Standard Python Installation
pip install polymerase-tm

# Or via Conda/Mamba
mamba install -c conda-forge polymerase-tm
```

## Quick Start

### Python API

```python
from polymerase_tm import tm, ta, list_buffers

# Single primer Tm (Q5, 500 nM)
print(tm("ATGTCCCTGCTCTTCTCTCGATGCAA"))          # 72

# Primer pair Ta
result_ta, tm_fwd, tm_rev = ta(
    "ATGTCCCTGCTCTTCTCTCGATGCAA",
    "GTGCCTCCGAGCCAGCACC",
)
print(f"Ta = {result_ta}, Fwd Tm = {tm_fwd}, Rev Tm = {tm_rev}")
# Ta = 72, Fwd Tm = 72, Rev Tm = 75

# Different polymerase
print(tm("ATGTCCCTGCTCTTCTCTCGATGCAA", polymerase="taq"))

# Override buffer (when not using default master mix)
print(tm("ATGTCCCTGCTCTTCTCTCGATGCAA", polymerase="taq", buffer="thermopol"))

# Direct salt concentration (mM)
print(tm("ATGTCCCTGCTCTTCTCTCGATGCAA", salt_mM=50))
```

### Automation & Batch Processing

```python
from polymerase_tm import (
    batch_tm,                 # Bulk Tm for many sequences
    optimal_binding_length,   # Find shortest binding region for target Tm
    check_pair,               # Full primer pair compatibility report
    pcr_protocol,             # Generate complete PCR cycling protocol
)

# Batch Tm for multiple primers
results = batch_tm(["ATCGATCGATCG", "GCGCGCGCGCGC", "AATTCCGGAATT"])

# Primer pair compatibility check (includes auto additive recommendation)
pair = check_pair("ATGTCCCTGCTCTTCTCTCGATGCAA", "GTGCCTCCGAGCCAGCACC")
if pair["additive"]["recommended"]:
    print(f"Use {pair['additive']['additive']} ({pair['additive']['concentration']})")

# Generate full PCR cycling protocol
protocol = pcr_protocol(
    "ATGTCCCTGCTCTTCTCTCGATGCAA",
    "GTGCCTCCGAGCCAGCACC",
    template="...full template sequence...",  # Auto-calculates amplicon_length
)
```

## Disclaimer

This package is not affiliated with New England Biolabs (NEB).
The algorithms exactly reproduce the calculations of the [NEB Tm Calculator](https://tmcalculator.neb.com/) and [NEB Base Changer](https://nebasechanger.neb.com/) for research and educational purposes.
Always verify critical calculations against the official tools.
