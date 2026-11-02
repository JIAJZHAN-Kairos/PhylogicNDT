# SASH/PURPLE review fixes — 2026-11-02

These changes repair seven reproducible issues found at `a70a1db`.
The previous `accepted_source_files.json` and `synthetic_smoke.json` are
historical records and remain unchanged. The current implementation is
identified by `provenance/repaired_source_files_20261102.json`.

| Issue | Repair | Regression coverage |
|---|---|---|
| Arm CN CCF became 1 after discretization | Retain continuous absolute allele CN, average it by covered bases within the selected arm state, then compute `(observed CN - event baseline)/(modal CN - event baseline)`. Use the doubled baseline for WGD-relative events; incompatible directions remain unresolved | Autosomal gain/loss, male X, segment-length weighting, WGD rebasing |
| Gain intervals used a different scale from the point estimate | Divide both estimate and uncertainty by the CN change from baseline | Autosomal/X/Y and higher-CN intervals; identical fragmented gains merge into an arm event |
| PCGR indels lost their annotations | Prefer `ALLELE_NUM`; otherwise accept only uniquely mapped raw, anchor-trimmed or minimal VEP alleles | Insertions/deletions, multiallelic deletions, ambiguous aliases retained as unknown |
| Focal mutation estimates did not reach output | Reconcile event-specific mutation posteriors into canonical mutations; discordant local estimates remain unknown with `conflicting_local_gain_estimates` | Actual timing TSV and event comparisons; reversed conflicting-event order |
| Multiallelic depth was reduced | Use all non-target-ALT reads as binomial failures; preserve literal VCF counts separately | Full AD depth and CCF in the core MAF reader; malformed negative AD rejected |
| Closed segment endpoints were lost | Adapter converts `[start,end]` to `[start,end+1)` once; core interval semantics stay unchanged | Adjacent segments, a single-base segment, PAR/non-PAR boundary |
| LeagueModel ignored serial seed | Seed Python and NumPy before model initialization; retain chunk-index order when collecting parallel results; seed 0 is valid | Repeat complete serial and parallel CLI runs and compare odds and positions |

## Input and reproducibility contracts

For a multiallelic `AD=70,20,10`, the prepared rows now contain:

| Target ALT | t_alt_count | t_ref_count (effective non-target count) | vcf_ref_count | vcf_other_alt_count |
|---|---:|---:|---:|---:|
| ALT 1 | 20 | 80 | 70 | 10 |
| ALT 2 | 10 | 90 | 70 | 20 |

The core uses `t_ref_count + t_alt_count` as binomial depth, so both rows use
100 reads. Original REF and other-ALT counts remain in the prepared MAF for
audit. Core Cluster/Timing outputs carry the effective counts. Biallelic
input counts are unchanged. Split alleles share a locus and are not thereby
independent biological events.

The prepared segment bridge is one-based half-open; the input PURPLE file
remains one-based closed. Old prepared files do not gain the new endpoint or
count fixes merely by rerunning the core. Regenerate preparations in a new
directory and run all stages into new result directories. Input provenance
records the new coordinate and count contracts.

For repeatable LeagueModel results, hold the runtime, inputs, `--seed`,
`PYTHONHASHSEED`, job count and chunk size fixed. Serial and parallel runs use
different random-stream schedules; cross-mode numerical equality is not
promised. Worker completion order no longer changes permutation IDs.

The uncertainty remains a support-based proxy, not a calibrated confidence
interval. CN-state selection and WGD presence/ordering models are otherwise
unchanged. In particular, the inherited high-copy/WGD rules are model
assumptions and do not establish ecDNA formation time.

## Validation

`tests/test_review_regressions.py` adds 13 tests to the 36 existing tests.
The reviewed source fails the new defect checks; the repair passes all 49
tests with the versions pinned in `requirements.txt` on Python 3.12.

```bash
python -m unittest discover -s tests -v
python provenance/verify_source.py --phylogic-root .
DEMO_ROOT=/absolute/path/to/new-demo bash examples/run_demo.sh
```

The full demo runs preparation, Cluster, Timing with `--comparison_method
original`, comparison aggregation and LeagueModel. It checks posterior
normalization, retained sex/PAR behavior, comparison probabilities, reports,
and saved LeagueModel output. This is synthetic implementation validation;
existing patient/cohort results are not overwritten or revalidated by it.
