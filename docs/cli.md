# Command Line Interface (CLI)

The package provides a comprehensive `polymerase-tm` CLI tool for calculating Tms and performing other useful sequence analysis tasks directly from your terminal.

## Basic Usage

### Tm / Ta Calculations
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
```

### Advanced Analysis

```bash
# DMSO analysis with template
polymerase-tm --dmso-check --template plasmid.gbk FWD_SEQ REV_SEQ

# Generate PCR Cycler Protocol & Virtual Agarose Gel Plot
polymerase-tm ATGTCCCTGCTCTTCTCTCGATGCAA GTGCCTCCGAGCCAGCACC --template plasmid.gbk --plot-gel out_gel.png

# Generate a standalone virtual gel with multiple custom fragments and a 10kb ladder
polymerase-tm --plot-gel out_multi.png --ladder 10kb --plot-gel-sizes 150 200 400

# Simulate a custom gel run (e.g. 1.5% agarose, 110V, 90 minutes, 15cm gel)
polymerase-tm --plot-gel custom_physics.png --ladder 1kb_plus --plot-gel-sizes 1500 3000 --agarose 1.5 --voltage 110 --time 90 --gel-length 15.0

# Compare topological isomers of the same size plasmid (Linear vs Supercoiled vs Nicked)
polymerase-tm --plot-gel topologies.png --plot-gel-sizes 3000 3000 3000 --topology linear coiled nicked
```

### Site-Directed Mutagenesis

```bash
# Site-directed mutagenesis (NEB Base Changer)
polymerase-tm --sdm --mutation M1A TEMPLATE_SEQ
polymerase-tm --sdm --mutation "T2A K5R" --codon-mode parsimony TEMPLATE_SEQ
polymerase-tm --sdm --mode del --mutation 10:3 TEMPLATE_SEQ
polymerase-tm --sdm --mode ins --mutation 10:AAAAAA TEMPLATE_SEQ
polymerase-tm --sdm --mutation T2A --genetic-code 2 TEMPLATE_SEQ
```

## CSV Batch Automation

You can batch process hundreds of sequences natively with the CLI tool. Depending on the `action` used, the CSV file requires specific column names. You can override these defaults using `--csv-<type>-col` (e.g., `--csv-template-col DNA_Sequence`).

```bash
# General usage
polymerase-tm --csv input.csv --csv-action action_name --csv-out output.csv

# Examples
polymerase-tm --csv pairs.csv --csv-action check_pair
polymerase-tm --csv primers.csv --csv-action tm --csv-out tms_out.csv
polymerase-tm --csv mutations.csv --csv-action sdm
polymerase-tm --csv protocols.csv --csv-action protocol
```

### Table of Actions

| Action | Required Columns | Optional Columns |
| :--- | :--- | :--- |
| `check_pair` | `fwd`, `rev` | `name` |
| `tm` | `seq` | `name` |
| `protocol` | `fwd`, `rev` | `template`, `name` |
| `sdm` | `template`, `mutation` | `mode` (default: point), `name` |

For more SDM details, see the [Mutagenesis](mutagenesis.md) page.

### Utility

```bash
# List all polymerases / buffers
polymerase-tm --list
polymerase-tm --list-buffers

# Version
polymerase-tm --version
```
