# Repaired core: scope, provenance and limits

## Repository and license scope

This branch contains both the repaired inference engine and the independent
SASH/PURPLE input adapters. The engine, inherited code and core regression
fixtures retain the [upstream PhylogicNDT LICENSE](../LICENSE), including
its non-commercial research and redistribution terms. The independently
written `workflow/` and `examples/` scripts, `tests/test_workflow.py` and the
new SASH/PURPLE workflow documentation are covered by [workflow/LICENSE](../workflow/LICENSE).
That MIT license does not relicense the inference engine. Original authorship
and attribution are retained; this fork does not claim upstream endorsement.

The published matching core packages the **sex-aware repair V2, 2026-10-27**,
including the preceding timing repairs. Its baseline was
`https://github.com/phylyc/PhylogicNDT.git` at
`4564634affb7a9ce6f1f2297979f80c501ef874a`, plus local adaptations.
This commit alone is **not** the repaired implementation. The upstream
`broadinstitute/PhylogicNDT` repository describes the original tool family.

The exact accepted snapshot is identified by
[`accepted_source_files.json`](../provenance/accepted_source_files.json), covering
47 source/reference files. No clinical identifiers, cohort manifests or
mutation/copy-number records are in that hash list. The legacy core
`SOURCE_COMMIT` file lists older adaptations only; the full file hashes take
precedence for version identity. Check this repository's core using:

```bash
python3 workflow/check_core.py --phylogic-root /path/to/repaired/core
```

## What changed beyond input conversion

These are implementation adaptations/repairs to a published method, not a
new independently benchmarked clonal-evolution algorithm.

| Area | Accepted repair/behaviour |
|---|---|
| hg38 compatibility | Normalise input chromosome labels to the selected genome; retain local allele-specific CN through Cluster/Timing |
| Purity and ploidy | Carry validated values through the input/CCF pipeline; retain the existing invertible PURPLE-to-allelic-capseg bridge |
| CN/cluster reconciliation | Use saved Cluster posterior modes/support rather than silently deriving clonal identity from a single CCF point estimate; conflicting/unsupported reconstructions remain explicit |
| Driver annotation | Handle protein-change syntax and overlapping driver representations without counting equivalent records as independent evolutionary events |
| Uninformative event times | Keep presence with unknown ordering; do not manufacture timing from flat prior-only distributions, unsupported high-copy states or no matching CN reconstruction |
| Sex-specific normal CN | Male non-PAR X/Y have one normal homolog; female X has two and female Y has none. Apply region-specific normal CN in CCF, multiplicity and event calibration |
| True X/Y changes | Do not call loss of a normally absent allele; retain observed gain/loss of existing homologs. There is no blanket X/Y filter |
| Negative X/Y estimates | Bound physically impossible estimates at zero in the corrected core and retain original records/reasons in QC |
| PAR | Split hg38 male X/Y intervals at PAR boundaries. Unresolved reference/alignment provenance gives no assigned CCF/time, while original records remain auditable |
| WGD clock evidence | Preserve WGD-neutral doubled states as clock candidates, pool matched clonal mutations once per arm/region, and avoid treating the two allele labels as independent clocks |
| Uninformative WGD | Retain the event with unknown time if no informative clock exists; it must not create quantitative downstream timing via an artificial uniform fallback |
| WGD definitions | Keep PURPLE WGD grouping separate from native Timing calls; do not force agreement or an earlier WGD time |
| Comparisons | The accepted analysis explicitly uses `--comparison_method original`; eligible-event comparison calculations retain the original method, with uninformative events left unknown |
| LeagueModel | Preserve native automatic selection, arm/focal separation and true sex-chromosome candidates; retain roundoff safeguards and deterministic seed settings |
| Reports | Distinguish presence from estimable time, include Yp/Yq axes, and prevent large timing comparison graphs from overwhelming report generation |

The input-adapter changes are narrower: explicit VCF tumour-column
selection; manifest-relative path resolution; validation of sample mapping,
sex and CN files; per-specimen SIFs; immutable preparation directories;
source hashes and preparation/comparison QC. **The adapters prepare inputs; the core modules in this same repository
implement the CCF/Timing repairs.**

## Sex/PAR assumptions that remain important

- Supported normal states are validated MALE/FEMALE baselines; unusual
  germline CN/sex karyotypes require separate review.
- `par_mode=unresolved` is the default. The analysis did not have the original
  DRAGEN/SASH reference/PAR configuration and matched-normal CN evidence
  needed to resolve male PAR. Do not infer that information from a filename.
- `x_diploid_y_masked` and `separate_haploid` are explicit input contracts,
  not modes to select by trial and error. Use them only when the reference
  and matched-normal interpretation are established.
- Retained X/Y candidates are not automatically validated biological
  sex-chromosome alterations. Review coverage, alignment/reference and
  specimen metadata, especially for low-support focal/PAR candidates.
- A complete loss of the only homolog does not by itself prove a pre-WGD
  event. A no-time/presence-only result is valid where ordering evidence
  is missing.

## Reproducibility boundary

The independent adapter can be tested publicly with fabricated data. Exact
reproduction of the historical analysis additionally requires the accepted
core under its upstream license, matching runtime, the original permitted study inputs,
driver/peak definitions, cohort/exclusion manifests and analysis parameters.
The repository does not supply patient data, a study cohort, an ecDNA
classifier dataset or the study's GISTIC peak reference.

The original mutation preparation used generic variant-class labels and a
support-based segment uncertainty proxy. These are documented compatibility
choices, not new measured CN confidence estimates. New downstream users
must state these assumptions rather than imply that PURPLE directly
provided the proxy uncertainty.
