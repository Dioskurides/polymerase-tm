# API Reference

This section details the Python implementation of the NEB Tm Calculator, and documents all public functions provided by `polymerase-tm`. Most user features are available directly in the top-level namespace.

```python
import polymerase_tm

# Example
polymerase_tm.tm("ATGC...")
polymerase_tm.ta("ATGC...", "CGAT...")
```

## Core

Calculates melting and annealing temperatures based on nearest-neighbor thermodynamics and Owczarzy salt corrections.

::: polymerase_tm.core.tm
::: polymerase_tm.core.ta
::: polymerase_tm.core.owczarzy_correction
::: polymerase_tm.core.owczarzy_bivariate
::: polymerase_tm.core.calc_sdm_tm
::: polymerase_tm.core.list_polymerases
::: polymerase_tm.core.list_buffers

## Site-Directed Mutagenesis

Implements the NEB Base Changer logic for primer design.

::: polymerase_tm.mutagenesis.BaseChanger
::: polymerase_tm.mutagenesis.SDMPrimer
::: polymerase_tm.mutagenesis.MutagenesisResult
::: polymerase_tm.mutagenesis.select_codon

## Batch Operations and Automation

Pipeline processing and programmatic generation of complex protocols.

::: polymerase_tm.batch.batch_tm
::: polymerase_tm.batch.check_pair
::: polymerase_tm.batch.pcr_protocol
::: polymerase_tm.batch.optimal_binding_length
::: polymerase_tm.batch.from_csv
::: polymerase_tm.batch.to_csv
::: polymerase_tm.batch.additive_recommendation

## Primer Analysis

Advanced checks for secondary structure and sequence compatibility.

::: polymerase_tm.analysis.primer_dimer
::: polymerase_tm.analysis.primer_quality
::: polymerase_tm.analysis.restriction_scan
::: polymerase_tm.analysis.gibson_overlaps

## Secondary Structure Analysis (DMSO/GC)

Evaluates the amplicon for complex features that require PCR additives.

::: polymerase_tm.dmso.find_hairpins
::: polymerase_tm.dmso.primer_hairpin
::: polymerase_tm.dmso.gc_content
::: polymerase_tm.dmso.find_gc_hotspots
::: polymerase_tm.dmso.analyze_amplicon
::: polymerase_tm.dmso.dmso_recommendation
::: polymerase_tm.dmso.print_dmso_report

## Visualization

::: polymerase_tm.gel.plot_virtual_gel
