# Codebase Audit Plan

**Goal:** Map the entire codebase, identify dead code, redundant files, and candidates for compaction. Audit only — no deletions until reviewed.

## Approach

Series of subagent cron jobs, each auditing one directory. Each agent gets:
- The directory listing
- Access to read files
- Instructions to produce a structured report

Reports land in `audit/` directory for human review.

## Phases

### Phase 1: Core Model (`src/`)
- Map every module and its public API
- Trace import graph (what imports what)
- Identify unused functions/classes
- Flag any dead code paths
- Note code duplication

### Phase 2: Tests (`tests/`)
- Map test files to source modules they test
- Identify orphaned tests (testing removed/renamed code)
- Flag redundant test cases
- Check coverage gaps (modules without tests)

### Phase 3: Scripts (`scripts/`)
- Classify each script: one-off vs reusable
- Check which scripts are referenced from docs/experiments
- Identify scripts that duplicate functionality
- Flag scripts with hardcoded paths or stale assumptions (e.g., 489-node references)

### Phase 4: Experiments (`experiments/`)
- Map experiment files and their dependencies
- Identify completed/obsolete experiments
- Check for hardcoded values that should be config
- Flag experiments still referencing old network sizes

### Phase 5: Conservation Module (`conservation/`)
- Already reviewed (REVIEW.md exists) — light audit
- Check for stale imports or unused utilities
- Verify test coverage

### Phase 6: Specs & Docs (`specs/`, `paper/`)
- Identify outdated specs superseded by implementation
- Flag specs that reference old architecture (pre-three-trait, pre-896-node)
- Check paper .tex files for stale content

### Phase 7: Data & Config (`data/`, `scenarios/`, config files)
- Identify orphaned data files
- Check scenario configs for stale parameters
- Map which data files are actually loaded by code

### Phase 8: Synthesis
- Merge all phase reports
- Produce single `audit/SUMMARY.md` with:
  - Total file count and lines of code
  - Dead code candidates (with reasoning)
  - Duplication candidates
  - Compaction recommendations
  - Dependency graph

## Output

Each phase writes to `audit/phase_N_report.md`. Final synthesis in `audit/SUMMARY.md`.

## When to Run

After calibration stabilizes. Not urgent — correctness > tidiness right now.
