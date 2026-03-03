# polymerase-tm

[![PyPI](https://img.shields.io/pypi/v/polymerase-tm)](https://pypi.org/project/polymerase-tm/)
[![Tests](https://github.com/Dioskurides/polymerase-tm/actions/workflows/test.yml/badge.svg)](https://github.com/Dioskurides/polymerase-tm/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/pypi/pyversions/polymerase-tm)](https://pypi.org/project/polymerase-tm/)

Exact Python reproduction of the [NEB Tm Calculator](https://tmcalculator.neb.com/) and [NEB Base Changer](https://nebasechanger.neb.com/) for PCR primer melting temperature (Tm), annealing temperature (Ta), and site-directed mutagenesis (SDM) primer design.

## Features

- **Exact NEB Tm Calculator reproduction** -- algorithm reverse-engineered from the NEB Tm Calculator front-end source; verified against the official tool with 0 degC deviation across all tested sequences.
- **Zero dependencies for core functions** -- `tm()`, `ta()`, `batch_tm()`, `check_pair()`, and all primer analysis work without any external packages.
- **22 NEB polymerase products** with their specific buffer salt concentrations and Ta rules (Q5, Phusion, Taq, OneTaq, LongAmp, Vent, Deep Vent, and more).
- **Automatic additive recommendation** -- suggests Q5 High GC Enhancer or DMSO based on primer GC, hairpins, and amplicon analysis.
- **Batch processing** -- process hundreds of primer pairs from CSV files with full Tm/Ta/compatibility analysis.
- **PCR protocol generator** -- generates complete cycling protocols with polymerase-specific temperatures and extension times. Automatically generates **touchdown protocols** when primer Tm difference exceeds 3 degC.
- **Smart primer design** -- find the optimal binding length for a target Tm.
- **Primer dimer detection** -- checks 3'-end complementarity and self-dimer risk.
- **Gibson Assembly overlap design** -- generates full primers with overhangs for Gibson/HiFi Assembly.
- **Restriction site scanning** -- scans primers for ~120 NEB restriction enzyme sites with full IUPAC ambiguity code support (accepts enzyme names or custom dict).
- **Primer quality scoring** -- comprehensive 0-100 score evaluating GC clamp, runs, repeats, hairpins.
- **Hairpin detection** -- nearest-neighbor thermodynamic Tm for hairpin stems (SantaLucia 1998), G-T wobble pair tolerance, loop entropy penalty (Jacobson-Stockmayer).
- **Site-directed mutagenesis** -- reimplementation of the [NEB Base Changer v2.7.2](https://nebasechanger.neb.com/) primer design algorithm. Supports AA point mutations, nucleotide substitutions/deletions/insertions, 12 NCBI genetic codes, codon usage/parsimony selection, back-to-back primer design with Owczarzy (2008) bivariate salt correction (Na+ + Mg2+).
- **DMSO analysis** -- analyses primer hairpins, amplicon GC content, GC-rich hotspots, and template secondary structures (requires `[bio]` extra for GenBank files).
- **Virtual gel visualization** -- simulated agarose gel with realistic Ferguson-plot migration physics (requires `[viz]` extra).
- **CLI tool** (`polymerase-tm`) for quick calculations from the terminal.

**Optional dependencies:** Biopython (for GenBank template analysis), matplotlib + seaborn (for virtual gel).

## Installation

```bash
# Core package (zero dependencies)
pip install polymerase-tm

# With GenBank template support (adds biopython)
pip install polymerase-tm[bio]

# With virtual gel visualization (adds matplotlib, seaborn)
pip install polymerase-tm[viz]

# Everything
pip install polymerase-tm[all]

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

### Site-Directed Mutagenesis (NEB Base Changer)

```python
from polymerase_tm import BaseChanger

# Design SDM primers for a point mutation
template = "ATGACCATGATTACGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGG..."
bc = BaseChanger(template)
result = bc.point_mutation("T2A")
print(f"FWD: {result.forward.sequence} (Tm={result.forward.tm})")
print(f"REV: {result.reverse.sequence} (Tm={result.reverse.tm})")
print(f"Ta:  {result.ta} degC")

# Multiple mutations
results = bc.batch("T2A A3G K5R")

# Deletion / insertion
result = bc.deletion(start=10, length=3)
result = bc.insertion(position=10, insert_seq="AAAAAA")

# Parsimony codon selection (fewest base changes)
bc = BaseChanger(template, codon_mode="parsimony")

# Alternative genetic code (e.g. vertebrate mitochondrial)
bc = BaseChanger(template, genetic_code=2)
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
polymerase-tm -p taq ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

# Override buffer (e.g. when not using master mix)
polymerase-tm --buffer thermopol ATGTCCCTGCTCTTCTCTCGATGCAA

# Direct salt concentration override
polymerase-tm --salt 50 ATGTCCCTGCTCTTCTCTCGATGCAA

# With DMSO correction
polymerase-tm --dmso 3 ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC

# DMSO analysis with template
polymerase-tm --dmso-check --template plasmid.gbk FWD_SEQ REV_SEQ

# Generate PCR Cycler Protocol & Virtual Agarose Gel Plot
polymerase-tm ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC --template plasmid.gbk --plot-gel out_gel.png

# Generate a standalone virtual gel with multiple custom fragments and a 100bp ladder
polymerase-tm --plot-gel out_multi.png --ladder 100bp --plot-gel-sizes 150 200 400

# Simulate a custom gel run (e.g. 1.5% agarose, 110V, 90 minutes, 15cm gel)
polymerase-tm --plot-gel custom_physics.png --ladder 1kb_plus --plot-gel-sizes 1500 3000 --agarose 1.5 --voltage 110 --time 90 --gel-length 15.0

# Compare topological isomers of the same size plasmid (Linear vs Supercoiled vs Nicked)
polymerase-tm --plot-gel topologies.png --plot-gel-sizes 3000 3000 3000 --topology linear coiled nicked

# List all polymerases / buffers
polymerase-tm --list
polymerase-tm --list-buffers

# Version
polymerase-tm --version

# Site-directed mutagenesis (NEB Base Changer)
polymerase-tm --sdm --mutation M1A TEMPLATE_SEQ
polymerase-tm --sdm --mutation "T2A K5R" --codon-mode parsimony TEMPLATE_SEQ
polymerase-tm --sdm --mode del --mutation 10:3 TEMPLATE_SEQ
polymerase-tm --sdm --mode ins --mutation 10:AAAAAA TEMPLATE_SEQ
polymerase-tm --sdm --mutation T2A --genetic-code 2 TEMPLATE_SEQ
```

## API Reference

| Function | Description |
|:---|:---|
| `tm(seq, polymerase, buffer, salt_mM)` | Melting temperature for one primer |
| `ta(seq1, seq2, polymerase, dmso_pct, buffer, salt_mM)` | Annealing temperature for a primer pair |
| `list_polymerases()` | List all 22 supported polymerases |
| `list_buffers()` | List all 17 NEB buffers with salt concentrations |
| `batch_tm(sequences, polymerase)` | Batch Tm for multiple sequences |
| `check_pair(fwd, rev, polymerase)` | Pair compatibility + additive recommendation |
| `pcr_protocol(fwd, rev, polymerase, template, amplicon_length)` | Full PCR cycling protocol (auto-touchdown when Tm diff > 3 degC) |
| `optimal_binding_length(seq, target_tm, polymerase)` | Find shortest binding region for target Tm |
| `reverse_complement(seq)` | DNA reverse complement |
| `from_csv(path, polymerase)` | Read primer pairs from CSV, compute Tm/Ta |
| `to_csv(results, path)` | Write results to CSV |
| `primer_dimer(fwd, rev)` | Check 3' complementarity / dimer risk |
| `gibson_overlaps(fwd_bind, rev_bind, left_seq, right_seq)` | Design Gibson/HiFi Assembly primers |
| `restriction_scan(seq, enzymes)` | Scan for restriction sites (~120 NEB enzymes, IUPAC support, accepts name list or dict) |
| `primer_quality(seq)` | Quality score 0-100 with issues list |
| `find_hairpins(seq)` | Detect hairpin stem-loops with NN Tm |
| `primer_hairpin(seq)` | Hairpin scan tuned for primer-length sequences |
| `additive_recommendation(fwd, rev, polymerase)` | DMSO / GC Enhancer recommendation |
| `dmso_recommendation(fwd_bind, rev_bind, template_seq, template_file)` | Full DMSO/additive analysis |
| `print_dmso_report(report)` | Pretty-print a DMSO analysis report |
| `gc_content(seq)` | GC content as fraction |
| `plot_virtual_gel(amplicon_lengths, ladder_name, agarose_pct, ...)` | Simulated agarose gel image (requires `[viz]`) |
| `BaseChanger(template, orf_start, genetic_code, codon_mode, ...)` | SDM primer designer (NEB Base Changer v2.7.2) |
| `select_codon(target_aa, original_codon, mode, genetic_code)` | Codon selection (usage/parsimony) |
| `parse_aa_mutation(mutation_str)` | Parse AA mutation notation ("M1A", "K2R:CGC") |
| `calc_sdm_tm(seq, mono_mM, divalent_mM, primer_conc_nM)` | SDM Tm with bivariate salt correction |
| `owczarzy_bivariate(raw_tm, seq, mono_mM, divalent_mM)` | Owczarzy (2008) Na+/Mg2+ correction |
| `GENETIC_CODES` | All 12 NCBI genetic codes |
| `get_codon_table(code_id)` | Full codon table for any genetic code |
| `get_aa_to_codons(code_id)` | Amino acid → codons for any genetic code |

> **Buffer / salt override:** `tm()` and `ta()` accept optional `buffer` (NEB buffer name) or `salt_mM` (direct mM value) to override the polymerase default. Priority: `salt_mM` > `buffer` > polymerase default.

## Algorithm

| Component | Method | Reference |
|:---|:---|:---|
| Nearest-neighbor Tm | SantaLucia (1998) | PNAS 95:1460-5 |
| Salt correction (mono) | Owczarzy et al. (2004) | Biochemistry 43:3537-54 |
| Salt correction (bivariate) | Owczarzy et al. (2008) | Biochemistry 47:5336-53 |
| Ta rules | Polymerase-specific | NEB Tm Calculator v1.16 |
| SDM primer design | Back-to-back | NEB Base Changer v2.7.2 |
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
The algorithms were reverse-engineered from the publicly available JavaScript sources of the [NEB Tm Calculator](https://tmcalculator.neb.com/) and [NEB Base Changer](https://nebasechanger.neb.com/) for research and educational purposes.
Always verify critical calculations against the official tools.

## Module Structure

```
polymerase_tm/
├── __init__.py      # Re-exports (backward compatible)
├── constants.py     # NN_PARAMS, BUFFERS, POLYMERASES, GENETIC_CODES
├── core.py          # tm(), ta(), owczarzy_bivariate(), calc_sdm_tm()
├── mutagenesis.py   # BaseChanger, SDMPrimer, MutagenesisResult
├── dmso.py          # DMSO analysis, hairpins, GC analysis
├── batch.py         # batch_tm(), pcr_protocol(), CSV I/O
├── analysis.py      # restriction_scan(), primer_dimer(), primer_quality()
├── gel.py           # Virtual agarose gel visualization
└── cli.py           # Command-line interface (Tm/Ta + SDM)
```

All functions are re-exported from `__init__.py` — `from polymerase_tm import tm` and `from polymerase_tm import BaseChanger` work directly.

## License

MIT
