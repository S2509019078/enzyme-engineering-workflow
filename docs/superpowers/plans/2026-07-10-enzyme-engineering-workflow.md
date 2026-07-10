# Enzyme Engineering Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reusable Windows-first workflow that converts sequence-evolution candidates and structures into resumable FoldX scans, standardized tables and interpretable figures.

**Architecture:** A Python package validates inputs, normalizes online-tool exports, maps sequence positions to structures, schedules FoldX safely, parses outputs and builds reports. A PowerShell launcher and YAML/TSV configuration keep local paths outside code. Reporting is independent of computation so figures can be rebuilt without rerunning FoldX.

**Tech Stack:** Python 3.10+, PyYAML, pandas, Biopython, matplotlib, seaborn, openpyxl, pytest, Windows PowerShell 5.1+, licensed FoldX 5.

## Global Constraints

- Windows PowerShell is the primary execution environment.
- FoldX binary and license are user-supplied and never committed.
- EVcouplings/ESM-like scores and FoldX ΔΔG remain separate named measurements.
- FASTA-to-PDB position mapping is mandatory before structure mutations are scheduled.
- Complete results are resumable; partial or corrupt results are recorded and retried.
- Status-log write failures must not terminate FoldX computation.
- Raw conversations, personal paths, credentials, large outputs and third-party binaries are excluded.

---

## File Map

- `pyproject.toml`, `.gitignore`, `LICENSE`: packaging and repository policy.
- `README.md`: Chinese quick start and scientific interpretation.
- `config/config.example.yaml`, `config/enzymes.example.tsv`: portable project inputs.
- `enzymeflow.ps1`: Windows launcher.
- `src/enzymeflow/models.py`: records for enzymes, candidates, mappings and FoldX results.
- `src/enzymeflow/config.py`: YAML/TSV loading and validation.
- `src/enzymeflow/normalize.py`: adapters for EVcouplings/ESM-style CSV files.
- `src/enzymeflow/mapping.py`: FASTA/PDB extraction, alignment and residue mapping.
- `src/enzymeflow/mutations.py`: candidate expansion and FoldX mutation syntax.
- `src/enzymeflow/foldx.py`: RepairPDB/BuildModel scheduling and safe resume.
- `src/enzymeflow/state.py`: robust task/status journal.
- `src/enzymeflow/report.py`: parsing, classification, XLSX and figures.
- `src/enzymeflow/cli.py`: `check`, `normalize`, `map`, `prepare`, `run`, `report`, `status`.
- `tests/fixtures/`, `tests/test_*.py`: small deterministic tests.
- `docs/online-tools.md`, `docs/input-formats.md`, `docs/interpretation.md`, `docs/troubleshooting.md`: user guidance.

### Task 1: Package, configuration and input audit

**Files:** Create `pyproject.toml`, `.gitignore`, `config/config.example.yaml`, `config/enzymes.example.tsv`, `src/enzymeflow/__init__.py`, `src/enzymeflow/models.py`, `src/enzymeflow/config.py`, `tests/test_config.py`.

**Interfaces:** Produce `WorkflowConfig.from_yaml(path)`, `load_enzymes(path) -> list[EnzymeInput]`, and `audit_inputs(config, enzymes) -> AuditReport`.

- [ ] Test relative path resolution, unique enzyme names, missing FASTA/PDB, invalid chain, FoldX run count and optional online-result paths.
- [ ] Run config tests and verify import failure.
- [ ] Implement immutable dataclasses, path resolution against project root and structured audit issues with severity and repair text.
- [ ] Run tests; expect pass.
- [ ] Commit with `git commit -m "feat: add enzyme project configuration"`.

### Task 2: Online result normalization

**Files:** Create `src/enzymeflow/normalize.py`, `tests/fixtures/ev_matrix.csv`, `tests/fixtures/ev_long.csv`, `tests/test_normalize.py`, `docs/online-tools.md`.

**Interfaces:** Produce `detect_format(frame) -> OnlineFormat`, `normalize_scores(path, enzyme, source) -> DataFrame`, with columns `enzyme,source,sequence_position,wild_type,mutant,source_score,source_direction,selected`.

- [ ] Test long and matrix CSV forms, UTF-8 BOM, alternative column names, wild-type exclusion, numeric coercion and explicit score-direction configuration.
- [ ] Run normalization tests and confirm failure.
- [ ] Implement adapters that preserve original row identifiers and reject ambiguous score direction rather than guessing.
- [ ] Run tests; expect pass.
- [ ] Document export steps and commit with `git commit -m "feat: normalize evolution score exports"`.

### Task 3: Sequence-to-structure mapping

**Files:** Create `src/enzymeflow/mapping.py`, `tests/fixtures/enzyme.fasta`, `tests/fixtures/enzyme.pdb`, `tests/test_mapping.py`.

**Interfaces:** Produce `read_fasta(path) -> str`, `pdb_chain_sequence(path, chain) -> StructureSequence`, `map_positions(fasta, structure) -> list[ResidueMapping]`, and `validate_candidate_mapping(candidates, mappings) -> MappingReport`.

- [ ] Test exact mapping, missing N-terminal residues, insertion codes, chain mismatch, nonstandard residues and wild-type mismatch.
- [ ] Run mapping tests and confirm failure.
- [ ] Implement Biopython pairwise alignment with a deterministic best-alignment rule and explicit unmapped/ambiguous states.
- [ ] Run tests; expect pass.
- [ ] Commit with `git commit -m "feat: map sequence candidates to PDB residues"`.

### Task 4: Candidate expansion and FoldX naming

**Files:** Create `src/enzymeflow/mutations.py`, `tests/test_mutations.py`.

**Interfaces:** Produce `expand_saturation(sites, sequence) -> list[Candidate]`, `to_foldx_code(candidate, mapping) -> str`, and `write_individual_list(candidates, mappings, path)`.

- [ ] Test exactly 19 substitutions per valid site, no wild-type self-substitution, deduplication, mapped chain syntax and rejection of unmapped candidates.
- [ ] Run mutation tests and confirm failure.
- [ ] Implement the canonical amino-acid alphabet and preserve both sequence-space and PDB-space mutation labels.
- [ ] Run tests; expect pass.
- [ ] Commit with `git commit -m "feat: generate mapped mutation libraries"`.

### Task 5: Durable state and FoldX scheduler

**Files:** Create `src/enzymeflow/state.py`, `src/enzymeflow/foldx.py`, `tests/test_state.py`, `tests/test_foldx.py`.

**Interfaces:** Produce `StatusJournal.append(record)`, `safe_status(record)`, `FoldXRunner.repair(enzyme)`, `FoldXRunner.build_batch(enzyme, batch)`, and `inspect_foldx_output(path, expected_count) -> OutputState`.

- [ ] Test locked primary status log fallback, atomic backup journal, command arguments with spaces, RepairPDB completeness, BuildModel completeness, partial-result retry, completed-result skip and per-enzyme failure isolation.
- [ ] Run tests and confirm failure.
- [ ] Implement `subprocess.run` without shell interpolation, per-batch directories, stdout/stderr logs and backup JSONL when CSV logging fails.
- [ ] Run tests; expect pass.
- [ ] Commit with `git commit -m "feat: run FoldX with safe resume and logging"`.

### Task 6: FoldX parsing, classification and independent reports

**Files:** Create `src/enzymeflow/report.py`, `tests/fixtures/Average_test.fxout`, `tests/test_report.py`.

**Interfaces:** Produce `parse_average_fxout(path, name_map) -> DataFrame`, `classify_ddg(value, stable_max, neutral_max) -> str`, `build_workbook(frame, path)`, and `build_figures(frame, output_dir)`.

- [ ] Test FoldX header detection, mutation-name restoration, duplicate aggregation, negative/neutral/positive ΔΔG classification, pivot orientation and missing-cell rendering.
- [ ] Run report tests and confirm failure.
- [ ] Implement CSV/XLSX output, per-enzyme summaries, full heatmaps with web-style diverging colors, score distributions and top-candidate plots.
- [ ] Run tests; expect pass.
- [ ] Verify report generation works using fixture outputs without FoldX installed and commit with `git commit -m "feat: summarize and visualize FoldX scans"`.

### Task 7: CLI and Windows launcher

**Files:** Create `src/enzymeflow/cli.py`, `enzymeflow.ps1`, `tests/test_cli.py`.

**Interfaces:** Provide `enzymeflow check`, `normalize`, `map`, `prepare`, `run [--enzyme NAME] [--force]`, `report`, `status`, and `all`.

- [ ] Test help text, dry-run task counts, stage ordering, enzyme selection, report-only operation and nonzero exit for unresolved audit errors.
- [ ] Run CLI tests and confirm failure.
- [ ] Implement argparse routing and concise Chinese progress/error messages.
- [ ] Implement a PowerShell launcher that finds the repository root and forwards all arguments to `python -m enzymeflow.cli`.
- [ ] Run tests plus PowerShell parser validation; expect pass and commit with `git commit -m "feat: add Windows workflow entry point"`.

### Task 8: Documentation and reusable scenarios

**Files:** Create `README.md`, `LICENSE`, `docs/input-formats.md`, `docs/interpretation.md`, `docs/troubleshooting.md`; add small example inputs under `examples/minimal/`.

**Interfaces:** Deliver a new-project recipe that never requires editing Python code.

- [ ] Document installation, FoldX placement, project configuration, online export import, mapping review, dry-run, overnight run, resume, report rebuild and GitHub-safe data handling.
- [ ] Explain that EV/ESM scores prioritize sequence plausibility while FoldX ΔΔG estimates stability; include a decision table for combined interpretation.
- [ ] Document locked-log recovery, partial batches, mismatched residue numbering, CPY-like sequence/PDB offsets and safe runtime estimation.
- [ ] Execute the minimal example through every non-FoldX command; expect normalized candidates, a reviewed mapping, prepared mutations and fixture-based reports.
- [ ] Commit with `git commit -m "docs: deliver reusable enzyme engineering workflow"`.

### Task 9: Release verification

**Files:** Modify only files implicated by verification failures.

**Interfaces:** Produce a clean Git repository ready for GitHub upload.

- [ ] Run `python -m compileall -q src tests` and `python -m pytest -q`.
- [ ] Parse `enzymeflow.ps1` with PowerShell and run its help command.
- [ ] Search tracked files for the original username, desktop paths, conversation IDs, credentials and private URLs; expect no matches.
- [ ] Confirm fixture-based reports rebuild without FoldX installed.
- [ ] Run `git status --short`; expect clean after fixes.
- [ ] Commit verification fixes with `git commit -m "test: harden enzyme workflow release"` only if changes were required.

