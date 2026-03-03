## v1.0.0 — Production Release

### Breaking Changes
- **Optional dependencies:** `biopython`, `matplotlib`, and `seaborn` are no longer required dependencies. The core `tm()`, `ta()`, `batch_tm()`, `check_pair()`, and all primer analysis functions work with zero external dependencies. Install extras as needed:
  - `pip install polymerase-tm[bio]` — GenBank template file support (biopython)
  - `pip install polymerase-tm[viz]` — Virtual gel plotting (matplotlib, seaborn)
  - `pip install polymerase-tm[all]` — Everything
- **`_additive_recommendation` renamed to `additive_recommendation`** (public API). The old underscore-prefixed name is kept as a backward-compatible alias.

### New Features
- **IUPAC ambiguity codes in restriction scan:** `restriction_scan()` now correctly matches enzymes with degenerate recognition sites (N, R, Y, S, W, etc.) via regex expansion. Previously, enzymes like BlpI (GCTNAGC), SfiI (GGCCNNNNNGGCC), HincII (GTYRAC), and SgrAI (CRCCGGYG) could never match any sequence.
- `additive_recommendation()` is now a fully documented public API function with parameter descriptions and return type documentation.
- Clear `ImportError` messages when optional features are used without their dependencies installed, pointing users to the correct extras (`pip install polymerase-tm[bio]` / `[viz]`).

### Improvements
- Development Status upgraded from Beta to Production/Stable.
- `plot_virtual_gel()` now raises `ImportError` with a helpful message instead of silently returning `False` when visualization libraries are missing.

## v0.9.1

### New Features
- Added template argument to pcr_protocol and CLI for automatic amplicon length and extension time prediction.
- Smart 3-prime matching in analyze_amplicon to support primers with 5-prime overhangs (e.g. Gibson assembly, restriction sites).
