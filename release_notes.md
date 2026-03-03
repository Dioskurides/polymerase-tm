## v1.0.2 — CLI PCR Protocol & Code Quality

### New Features
- **PCR cycling protocol always shown for primer pairs:** The CLI now automatically generates and displays a full PCR cycling protocol whenever two primers are provided — no `--template` required. Extension time defaults to 30 s without a template; with `--template`, the amplicon length is auto-detected and extension time is calculated from the polymerase's extension rate.
- **Buffer/salt override in PCR protocol:** `pcr_protocol()` now accepts `buffer` and `salt_mM` parameters. The CLI correctly passes `--buffer` and `--salt` overrides through to the protocol, ensuring the displayed Ta matches the protocol's annealing temperature.
- **Buffer/salt display in CLI output:** When `--buffer` or `--salt` is specified, the values are shown in the output for both single-primer and primer-pair modes.

### Bug Fixes
- **Critical: CLI crash with touchdown protocols** — Fixed `ValueError` on `{step['temp']:>2d}` when touchdown annealing produces a string temp like `"75 -> 72"`. Now handles both int and string temperature values.
- **Input validation in `tm()` / `calc_nn_raw()`** — Empty sequences, invalid characters (non-ATGC), unknown polymerase keys, and zero/negative salt concentrations now raise descriptive `ValueError` instead of cryptic `KeyError`/`IndexError`/`ZeroDivisionError`.
- **CLI DNA validation** — Invalid primer characters are caught before computation with a helpful error message.
- **CLI edge-case warnings** — `--template` and `--dmso-check` with a single primer now print a warning instead of being silently ignored.

### Code Quality
- Removed unused imports: `math` from constants.py, `os` from gel.py, `Bio.Seq.Seq` from dmso.py.
- Removed dead code: unused `rc` variable in `restriction_scan()`.
- Fixed misleading comment in `gibson_overlaps()` ("50 nM" → uses polymerase default).
- Fixed `list_polymerases()` docstring (`conc_nM` → `primer_conc_nM`).
- Fixed `plot_virtual_gel()` docstring (documents `ImportError` raise instead of `False` return).
- Updated README API table: `fragments` → `amplicon_lengths`, `ladder` → `ladder_name`.
- Modernized type hints in gel.py: `Dict`/`List`/`Tuple` → `dict`/`list`/`tuple`.
- 67 tests passing.

---

## v1.0.1 — Hairpin NN & Touchdown PCR

### New Features
- **Nearest-neighbor hairpin Tm:** `find_hairpins()` and `primer_hairpin()` now compute Tm using SantaLucia (1998) NN parameters with Jacobson-Stockmayer loop entropy, replacing the Wallace rule. G-T wobble pairs are tolerated in stems ≥ 6 bp. Hairpin results now include a `mismatches` field.
- **Touchdown PCR protocol:** `pcr_protocol()` auto-generates touchdown protocols when primer Tm difference > 3 °C. Starts annealing at the higher Tm and decreases by 0.5 °C/cycle. Configurable via `touchdown`, `td_step`, `td_cycles` parameters.
- **Touchdown recommendation in `check_pair()`:** Instead of "redesign primers", large Tm differences now suggest using a touchdown protocol.

### Improvements
- Hairpin deduplication now sorts by Tm (highest first) instead of stem length.
- 65 tests passing.

---

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
