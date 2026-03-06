# Site-Directed Mutagenesis

The `polymerase-tm` package provides a robust implementation of the **NEB Base Changer** algorithm for designing primers for the Q5 Site-Directed Mutagenesis Kit. 

The `BaseChanger` class can automatically design optimized back-to-back primer pairs with balanced Tms for four distinct modes:
1. Amino acid point mutations
2. Nucleotide substitutions
3. Nucleotide deletions
4. Nucleotide insertions

The underlying Tm engine correctly utilizes the Q5 SDM conditions (50 mM monovalent salt, 2 mM Mg2+, 500 nM primer) combined with the Owczarzy et al. (2008) bivariate salt correction model.

## Quick Python Example

```python
from polymerase_tm import BaseChanger

# Design SDM primers for a point mutation
template = "ATGACCATGATTACGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGG..."
bc = BaseChanger(template)
result = bc.point_mutation("T2A")
print(f"FWD: {result.forward.sequence} (Tm={result.forward.tm})")
print(f"REV: {result.reverse.sequence} (Tm={result.reverse.tm})")
print(f"Ta:  {result.ta} degC")

# Multiple point mutations in the same run
results = bc.batch("T2A A3G K5R")

# Deletion / insertion
result = bc.deletion(start=10, length=3)
result = bc.insertion(position=10, insert_seq="AAAAAA")

# Parsimony codon selection (fewest base changes from original)
bc = BaseChanger(template, codon_mode="parsimony")

# Alternative genetic code (e.g. vertebrate mitochondrial)
bc = BaseChanger(template, genetic_code=2)
```

## Mutation Formatting for CLI and CSV Batch

When using the `--sdm` action in the `polymerase-tm` CLI or during a `from_csv` execution with `action="sdm"`, the `mutation` string describes exactly what is to be changed in the template. The necessary string structure depends heavily on the chosen `mode`.

| Mode | Format | Example | Description |
| :--- | :--- | :--- | :--- |
| **`point`** | `<OldAA><Pos><NewAA>` | `M1A` | Mutates Methionine at Amino Acid position 1 to Alanine. Chain commands with spaces: `T2A K5R`. |
| **`sub`** | `<Pos>:<Nucleotides>` | `10:GGG` | Replaces nucleotides starting at base pair 10 with `GGG`. |
| **`del`** | `<Pos>:<Length>` | `15:3` | Deletes exactly 3 nucleotides starting at base pair 15. |
| **`ins`** | `<Pos>:<Nucleotides>` | `20:AAAA` | Inserts four Adenines (`AAAA`) at base pair 20. |

!!! note "Default Mode"
    If `mode` is omitted, it defaults to **`point`**. In this mode, mutations happen on the translated **Amino Acid sequence**, not the naked DNA. To avoid off-by-one errors or frameshifts, ensure you provide the correct `--orf-start` (0-indexed base pair relative to the start of the `template` sequence). Default is 0.

## API Methods

Reference for the `BaseChanger` operations.

- [BaseChanger](api.md#polymerase_tm.mutagenesis.BaseChanger)
- [BaseChanger.point_mutation](api.md#polymerase_tm.mutagenesis.BaseChanger.point_mutation)
- [BaseChanger.substitution](api.md#polymerase_tm.mutagenesis.BaseChanger.substitution)
- [BaseChanger.insertion](api.md#polymerase_tm.mutagenesis.BaseChanger.insertion)
- [BaseChanger.deletion](api.md#polymerase_tm.mutagenesis.BaseChanger.deletion)
- [BaseChanger.batch](api.md#polymerase_tm.mutagenesis.BaseChanger.batch)
