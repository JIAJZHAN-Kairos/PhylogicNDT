# Source identity and validation

## Source identity

- This branch builds on the accepted sex-aware V2 snapshot dated 2026-10-27,
  with the [2026-11-02 review repairs](REVIEW_FIXES_20261102.md).
- Baseline lineage: `https://github.com/phylyc/PhylogicNDT.git`, commit
  `4564634affb7a9ce6f1f2297979f80c501ef874a`, plus the local changes described
  in [change categories and limitations](REPAIRS_AND_LIMITATIONS.md).
  That baseline commit alone does not identify this branch's implementation.
- The branch was created in the existing personal fork without changing
  `master`. It retains the upstream LICENSE and attribution.
- `SOURCE_COMMIT` is the legacy provenance record from the accepted source.
  Its local-patch list predates the full V2 changes; the complete file hashes
  in `provenance/accepted_source_files.json` identify that historical version.
  `provenance/repaired_source_files_20261102.json` identifies the seven-fix
  snapshot. The current deployment, including the Nextflow wrapper and audited
  negative-CN-to-zero input rule, is identified by
  `provenance/nextflow_source_files_20261102_v021.json`.

```bash
python provenance/verify_source.py --phylogic-root .
```

This checks the current repair snapshot. It is no longer byte-identical to
the 2026-10-27 inference code; the historical manifest and smoke receipt remain
unchanged for provenance. Hash verification identifies source, not correctness.

## Regression tests

```bash
python -m unittest discover -s tests -v
```

The 22 core regression tests cover supported sex-specific normal baselines,
true haploid changes, absent-homolog pseudo-losses, region-specific normal CN
in CCF/multiplicity, preservation of autosomal CCF calculations, WGD clock
evidence pooling and eligibility, unknown event times, driver/CN deduplication,
CN-cluster reconciliation and timing-graph size guards.

The original tests' local-path dependencies were replaced with versioned
historical reference functions under `tests/reference/`. These reference
files are test baselines only, not runtime inference modules. The original
comparison formula and pre-sex-aware CCF implementation provide checks for
the specific invariance claims. The root PhylogicNDT LICENSE applies to these
reference functions as well.

GitHub Actions installs the pinned dependencies on Linux / Python 3.12,
checks source hashes, runs the regression tests and checks the three CLI
entry points, plus the 14 input-adapter tests and the synthetic integration smoke test.

The review regression suite additionally covers measured arm-CN CCFs, WGD
rebasing, gain uncertainty and segment merging, PCGR indel allele mapping,
focal mutation-posterior propagation/conflicts, multiallelic read depth,
inclusive segment endpoints/PAR boundaries, and repeated serial/parallel
LeagueModel runs with a fixed seed.

## Synthetic integration check

The historical workflow ran the 2026-10-27 core snapshot locally on
fabricated male/female samples through Cluster (250 iterations), Timing
(`original`) and LeagueModel (2 permutations / 10 seasons). Its compact
[patient-free validation receipt](../provenance/synthetic_smoke.json)
records the actual test runtime and output counts.

The synthetic smoke test and regression tests are not a biological benchmark
or a full production cohort acceptance. New datasets still require input/QC
review, appropriate driver/peak definitions, explicit population/exclusion
manifests and output-level validation. Follow the
[SASH/PURPLE tutorial](SASH_PURPLE_TO_PHYLOGIC.md).
