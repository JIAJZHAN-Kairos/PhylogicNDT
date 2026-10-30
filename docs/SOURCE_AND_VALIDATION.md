# Source identity and validation

## Source identity

- This branch packages the accepted sex-aware repair V2 snapshot dated
  2026-10-27, including the earlier timing repairs.
- Baseline lineage: `https://github.com/phylyc/PhylogicNDT.git`, commit
  `4564634affb7a9ce6f1f2297979f80c501ef874a`, plus the local changes described
  in this branch. That baseline commit alone is not the repaired version.
- The branch was created in the existing personal fork without changing
  `master`. It retains the upstream LICENSE and attribution.
- `SOURCE_COMMIT` is the legacy provenance record from the accepted source.
  Its local-patch list predates the full V2 repairs; the complete file hashes
  in `provenance/accepted_source_files.json` identify the actual version.

```bash
python provenance/verify_source.py --phylogic-root .
```

This checks 47 source/reference files. Inference code is byte-identical to
the accepted snapshot. The new README, portable regression tests, historical
test fixtures and CI configuration are packaging additions, not changes to
the inference implementation.

## Regression tests

```bash
python -m unittest discover -s tests -v
```

The 22 repair regression tests cover sex-normal baselines/PAR boundaries,
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

## Synthetic integration check

The included workflow has already run this same core snapshot locally on
fabricated male/female samples through Cluster (250 iterations), Timing
(`original`) and LeagueModel (2 permutations / 10 seasons). Its compact
[patient-free validation receipt](../provenance/synthetic_smoke.json)
records the actual test runtime and output counts.

The synthetic smoke test and regression tests are not a biological benchmark
or a full production cohort acceptance. New datasets still require input/QC
review, appropriate driver/peak definitions, explicit population/exclusion
manifests and output-level validation. Follow the
[SASH/PURPLE tutorial](SASH_PURPLE_TO_PHYLOGIC.md).
