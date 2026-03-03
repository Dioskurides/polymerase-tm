# polymerase-tm

[![PyPI](https://img.shields.io/pypi/v/polymerase-tm)](https://pypi.org/project/polymerase-tm/)
[![Tests](https://github.com/Dioskurides/polymerase-tm/actions/workflows/test.yml/badge.svg)](https://github.com/Dioskurides/polymerase-tm/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/pypi/pyversions/polymerase-tm)](https://pypi.org/project/polymerase-tm/)

Exact Python reproduction of the [NEB Tm Calculator](https://tmcalculator.neb.com/) for PCR primer melting temperature (Tm) and annealing temperature (Ta) prediction.

## Features

- **Exact NEB Tm Calculator reproduction** -- algorithm recreated from the NEB Tm Calculator front-end source; verified against the official tool with 0 degC deviation across all tested sequences.
- **22 NEB polymerase products** with their specific buffer salt concentrations and Ta rules (Q5, Phusion, Taq, OneTaq, LongAmp, Vent, Deep Vent, and more).
- **Automatic additive recommendation** -- suggests Q5 High GC Enhancer or DMSO based on primer GC, hairpins, and amplicon analysis.
- **Batch processing** -- process hundreds of primer pairs from CSV files with full Tm/Ta/compatibility analysis.
- **PCR protocol generator** -- generates complete cycling protocols with polymerase-specific temperatures and extension times.
- **Smart primer design** -- find the optimal binding length for a target Tm.
- **Primer dimer detection** -- checks 3'-end complementarity and self-dimer risk.
- **Gibson Assembly overlap design** -- generates full primers with overhangs for Gibson/HiFi Assembly.
- **Restriction site scanning** -- scans primers for ~120 NEB restriction enzyme sites (accepts enzyme names or custom dict).
- **Primer quality scoring** -- comprehensive 0-100 score evaluating GC clamp, runs, repeats, hairpins.
- **DMSO analysis** -- analyses primer hairpins, amplicon GC content, GC-rich hotspots, and template secondary structures.
- **CLI tool** (`polymerase-tm`) for quick calculations from the terminal.

**Dependencies:** Biopython (for GenBank template analysis in DMSO features).

## Installation

```bash
pip install polymerase-tm

# From conda-forge (when available)
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

# List all available buffers
for b in list_buffers():
    print(f"{b['name']:20s} {b['salt_mM']:>4d} mM")

# Ta with 3% DMSO
ta_dmso, _, _ = ta("ATGTCCCTGCTCTTCTCTCGATGCAA", "GTGCCTCCGAGCCAGCACC", dmso_pct=3)
print(f"Ta with 3% DMSO = {ta_dmso}")  # 70

# List all available polymerases
from polymerase_tm import list_polymerases
for p in list_polymerases():
    print(f"{p['key']:25s} {p['description']}")
```

### Automation & Batch Processing

```python
from polymerase_tm import (
    batch_tm,                 # Bulk Tm for many sequences
    optimal_binding_length,   # Find shortest binding region for target Tm
    check_pair,               # Full primer pair compatibility report
    pcr_protocol,             # Generate complete PCR cycling protocol
    reverse_complement,       # DNA reverse complement
    from_csv, to_csv,         # CSV batch I/O
)

# Batch Tm for multiple primers
results = batch_tm(["ATCGATCGATCG", "GCGCGCGCGCGC", "AATTCCGGAATT"])
for r in results:
    print(f"{r['sequence']}: Tm={r['tm']} degC, GC={r['gc_pct']}%")

# Find optimal binding length for a target Tm
result = optimal_binding_length("ATGTCCCTGCTCTTCTCTCGATGCAA", target_tm=65)
print(f"{result['binding_seq']} ({result['length']} nt, Tm={result['tm']})")
# CCTGCTCTTCTCTCGATGCAA (21 nt, Tm=67)

# Primer pair compatibility check (includes auto additive recommendation)
pair = check_pair("ATGTCCCTGCTCTTCTCTCGATGCAA", "GTGCCTCCGAGCCAGCACC")
print(f"Ta={pair['ta']}, compatible={pair['compatible']}")
if pair["additive"]["recommended"]:
    print(f"Use {pair['additive']['additive']} ({pair['additive']['concentration']})")
    # -> "Use Q5 High GC Enhancer (1x)" for Q5 with high-GC primers
    # -> "Use DMSO (3%)" for Taq with high-GC primers

# Generate full PCR cycling protocol
protocol = pcr_protocol(
    "ATGTCCCTGCTCTTCTCTCGATGCAA",
    "GTGCCTCCGAGCCAGCACC",
    template="...full template sequence...",  # Auto-calculates amplicon_length
)
for step in protocol["cycling"]:
    print(f"{step['step']:25s} {step['temp']} degC  {step['time']}")
# Initial Denaturation       98 degC  30 s
# Denaturation               98 degC  10 s
# Annealing                  72 degC  30 s
# Extension                  72 degC  1 min 15 s
# Final Extension            72 degC  2 min
# Hold                        4 degC  indefinite

# CSV pipeline: read primers, compute everything, write results
results = from_csv("primers.csv")  # expects columns: name, fwd, rev
to_csv(results, "results_with_tm.csv")
```

### DMSO Analysis

```python
from polymerase_tm import dmso_recommendation, print_dmso_report

report = dmso_recommendation(
    fwd_bind="ATGTCCCTGCTCTTCTCTCGATGCAA",
    rev_bind="GTGCCTCCGAGCCAGCACC",
    template_file="template.gbk",     # optional GenBank template
)
print_dmso_report(report)
```

### Command Line

```bash
# Single primer Tm
polymerase-tm ATGTCCCTGCTCTTCTCTCGATGCAA

# Primer pair Ta (includes auto additive recommendation)
polymerase-tm ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

# Different polymerase
polymerase-tm --polymerase taq ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

# Override buffer (e.g. when not using master mix)
polymerase-tm --buffer thermopol ATGTCCCTGCTCTTCTCTCGATGCAA

# Direct salt concentration override
polymerase-tm --salt 50 ATGTCCCTGCTCTTCTCTCGATGCAA

# With DMSO correction
polymerase-tm --dmso 3 ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

# List all polymerases / buffers
polymerase-tm --list
polymerase-tm --list-buffers

# DMSO analysis with template
polymerase-tm --dmso-check --template template.gbk FWD_SEQ REV_SEQ

# Version
polymerase-tm --version
```

## API Reference

| Function | Description |
|:---|:---|
| `tm(seq, polymerase, buffer, salt_mM)` | Melting temperature for one primer |
| `ta(seq1, seq2, polymerase, dmso_pct, buffer, salt_mM)` | Annealing temperature for a primer pair |
| `list_polymerases()` | List all 22 supported polymerases |
| `list_buffers()` | List all 15 NEB buffers with salt concentrations |
| `batch_tm(sequences, polymerase)` | Batch Tm for multiple sequences |
| `check_pair(fwd, rev, polymerase)` | Pair compatibility + additive recommendation |
| `pcr_protocol(fwd, rev, polymerase, amplicon_length)` | Full PCR cycling protocol |
| `optimal_binding_length(seq, target_tm, polymerase)` | Find shortest binding region for target Tm |
| `reverse_complement(seq)` | DNA reverse complement |
| `from_csv(path, polymerase)` | Read primer pairs from CSV, compute Tm/Ta |
| `to_csv(results, path)` | Write results to CSV |
| `primer_dimer(fwd, rev)` | Check 3' complementarity / dimer risk |
| `gibson_overlaps(fwd_bind, rev_bind, left_seq, right_seq)` | Design Gibson/HiFi Assembly primers |
| `restriction_scan(seq, enzymes)` | Scan for restriction sites (~120 NEB enzymes, accepts name list or dict) |
| `primer_quality(seq)` | Quality score 0-100 with issues list |
| `dmso_recommendation(fwd, rev, template)` | Full DMSO/additive analysis |
| `gc_content(seq)` | GC content as fraction |

> **Buffer / salt override:** `tm()` and `ta()` accept optional `buffer` (NEB buffer name) or `salt_mM` (direct mM value) to override the polymerase default. Priority: `salt_mM` > `buffer` > polymerase default.

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
The algorithm was recreated from the publicly available JavaScript source of the [NEB Tm Calculator](https://tmcalculator.neb.com/) for research and educational purposes.
Always verify critical calculations against the official tool.

## Module Structure

```
polymerase_tm/
├── __init__.py      # Re-exports (backward compatible)
├── constants.py     # NN_PARAMS, BUFFERS, POLYMERASES
├── core.py          # tm(), ta(), list_buffers(), list_polymerases()
├── dmso.py          # DMSO analysis, hairpins, GC analysis
├── batch.py         # batch_tm(), pcr_protocol(), CSV I/O
├── analysis.py      # restriction_scan(), primer_dimer(), primer_quality()
└── cli.py           # Command-line interface
```

All functions are re-exported from `__init__.py` — `from polymerase_tm import tm` works as before.

## License

MIT
