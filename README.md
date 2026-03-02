# polymerase-tm

Exact Python reproduction of the [NEB Tm Calculator](https://tmcalculator.neb.com/) for PCR primer melting temperature (Tm) and annealing temperature (Ta) prediction.

## Features

- **Exact NEB Tm Calculator reproduction** -- algorithm reverse-engineered from the NEB Tm Calculator front-end source; verified against the official tool with 0 degC deviation across all tested sequences.
- **22 NEB polymerase products** with their specific buffer salt concentrations and Ta rules (Q5, Phusion, Taq, OneTaq, LongAmp, Vent, Deep Vent, and more).
- **DMSO recommendation** -- analyses primer hairpins, amplicon GC content, GC-rich hotspots, and template secondary structures.
- **CLI tool** (`polymerase-tm`) for quick calculations from the terminal.
- **No dependencies** for core Tm/Ta calculations. Biopython is optional (only needed for reading GenBank template files in DMSO analysis).

## Installation

```bash
# From PyPI
pip install polymerase-tm

# With Biopython support (for GenBank template analysis)
pip install polymerase-tm[bio]

# From conda-forge (when available)
mamba install -c conda-forge polymerase-tm
```

## Quick Start

### Python API

```python
from polymerase_tm import tm, ta, dmso_recommendation

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

# Ta with 3% DMSO
ta_dmso, _, _ = ta("ATGTCCCTGCTCTTCTCTCGATGCAA", "GTGCCTCCGAGCCAGCACC", dmso_pct=3)
print(f"Ta with 3% DMSO = {ta_dmso}")  # 70

# List all available polymerases
from polymerase_tm import list_polymerases
for p in list_polymerases():
    print(f"{p['key']:25s} {p['description']}")
```

### DMSO Analysis

```python
from polymerase_tm import dmso_recommendation, print_dmso_report

report = dmso_recommendation(
    fwd_bind="ATGTCCCTGCTCTTCTCTCGATGCAA",
    rev_bind="GTGCCTCCGAGCCAGCACC",
    template_file="template.gbk",     # optional, requires biopython
)
print_dmso_report(report)
```

### Command Line

```bash
# Single primer Tm
polymerase-tm ATGTCCCTGCTCTTCTCTCGATGCAA

# Primer pair Ta
polymerase-tm ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

# Different polymerase
polymerase-tm --polymerase taq ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

# With DMSO correction
polymerase-tm --dmso 3 ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

# List all polymerases
polymerase-tm --list

# DMSO analysis with template
polymerase-tm --dmso-check --template template.gbk ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC
```

## Algorithm

| Component | Method | Reference |
|:---|:---|:---|
| Nearest-neighbor Tm | SantaLucia (1998) | PNAS 95:1460-5 |
| Salt correction | Owczarzy et al. (2004) | Biochemistry 43:3537-54 |
| Ta rules | Polymerase-specific | NEB Tm Calculator v1.16 |
| DMSO correction | -0.6 degC per 1% | NEB Tm Calculator v1.16 |

### Buffer Salt Concentrations

| Buffer | [Monovalent] (mM) | Used by |
|:---|:---|:---|
| Q5 | 150 | Q5, Q5 Hot Start, Q5 Blood Direct |
| Q5U | 170 | Q5U Hot Start |
| Q5 Master Mix | 150 | Q5 2X Master Mix |
| Phusion HF / GC | 222 | Phusion, Phusion Hot Start Flex |
| Standard Taq | 55 | Taq, Hot Start Taq, EpiMark |
| ThermoPol | 40 | Vent, Deep Vent |
| OneTaq Std | 54 | OneTaq (Standard Buffer) |
| OneTaq GC | 80 | OneTaq (GC Buffer) |
| LongAmp | 100 | LongAmp, LongAmp Hot Start |
| Crimson Taq | 55 | Crimson Taq |
| Hemo KlenTaq | 70 | Hemo KlenTaq |
| Multiplex | 90 | Multiplex PCR Master Mix |

### Ta Calculation Rules

| Polymerase family | Rule | Cap |
|:---|:---|:---|
| Q5 | min(Tm1, Tm2) + 1 | 72 degC |
| Q5U | min(Tm1, Tm2) + 2 | 72 degC |
| Phusion | 0.93 * min(Tm1, Tm2) + 7.5 | 72 degC |
| Taq / OneTaq | min(Tm1, Tm2) - 5 | 68 degC |
| LongAmp | min(Tm1, Tm2) - 5 | 65 degC |
| Vent / Deep Vent | min(Tm1, Tm2) - 2 | 72 degC |

## Disclaimer

This package is not affiliated with New England Biolabs (NEB).
The algorithm was reverse-engineered from the publicly available JavaScript source of the [NEB Tm Calculator](https://tmcalculator.neb.com/) for research and educational purposes.
Always verify critical calculations against the official tool.

## License

MIT
