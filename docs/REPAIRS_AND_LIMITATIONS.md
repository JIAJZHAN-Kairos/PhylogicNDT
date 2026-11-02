# Changes: upstream bugs, input adaptations and project choices

The [2026-11-02 review fixes](REVIEW_FIXES_20261102.md) update CN CCF
propagation/uncertainty, indel annotation, focal mutation output, multiallelic
depth, segment endpoints and RNG handling. The snapshot description below
records the earlier baseline; current hashes are in
`provenance/repaired_source_files_20261102.json`.

This branch combines a confirmed upstream code-bug fix, SASH/PURPLE and hg38
input adaptations, and project-specific extensions and conservative rules.
These categories should be distinguished in Methods and when sharing the
workflow. Not every difference from upstream is a bug fix, and the changes
are not an independently benchmarked clonal-evolution algorithm.

## 1. Confirmed upstream code bug

`SomMutation.clean_local_cn()` in the
[baseline source](https://github.com/phylyc/PhylogicNDT/blob/4564634affb7a9ce6f1f2297979f80c501ef874a/data/SomaticEvents.py#L198)
used:

```python
self.local_cn_a1 = float(cn1) if not np.nan else np.nan
```

The condition tests the constant `np.nan`, not `cn1`. It is always false,
so calling this function replaces even a valid CN with NaN. The same error
affected the second allele. The current implementation converts the supplied
values and handles missing/invalid input explicitly.

This is a directly verified implementation defect. Its presence does not
establish that every upstream analysis invoked this path or was affected.

## 2. SASH/PURPLE and hg38 adaptations

| Area | Adaptation and attribution |
|---|---|
| SASH variant input | Select the exact tumour VCF column, map each alternate allele to its AD count, and extract SASH/PCGR gene/protein annotations. This is an input adapter, not a new mutation caller |
| PURPLE CN input | Attach absolute minor/major allele CN to mutations and encode segments through the existing invertible allelic-capseg bridge. The support-based CN uncertainty proxy is an adapter convention, not a measured PURPLE confidence interval |
| Chromosome labels | Match `1`/`chr1` and X/Y aliases to the selected hg38 genome dictionary. Exact-key mismatches in our input path could skip records; this is compatibility handling, not evidence that correctly formatted upstream inputs fail |
| Local CN between stages | Preserve `Allelic_CN_minor`/`Allelic_CN_major` when reading Cluster output into Timing. The baseline reader did not pass these fields into the mutation object; the change avoids dependence on reconstructing them later from segments |
| Purity and ploidy | Pass actual sample values to CN decoding. The baseline call omitted them and could use defaults of purity 1 and ploidy 2 in this input path. This is a parameter-propagation fix for the selected integration, not a claim about every upstream input type |
| Sample and file contracts | Generate independent per-specimen SIFs, resolve manifest-relative paths, validate sample columns and CN segments, preserve existing preparations, and record input hashes |
| Group comparisons | Aggregate explicit specimen lists, normalise comparison-table headers, and reject specimen mismatches, duplicate event pairs and invalid probabilities |

The baseline
[Cluster-output reader](https://github.com/phylyc/PhylogicNDT/blob/4564634affb7a9ce6f1f2297979f80c501ef874a/data/Sample.py#L214)
and [CN-parser call](https://github.com/phylyc/PhylogicNDT/blob/4564634affb7a9ce6f1f2297979f80c501ef874a/data/Sample.py#L99)
provide the implementation evidence for these interface changes. hg38 label
normalisation does not perform a genome-build liftover; mutation, segment,
centromere, cytoband and focal-region coordinates must already use hg38.

## 3. Project extensions and conservative handling

These are local or inherited workflow choices. They can change event
eligibility, CCFs and inferred times; they are not all established upstream
bugs. The original tool also accepts precomputed CCFs, so its input
assumptions should not be inferred from our SASH/PURPLE path alone.

| Area | Project behaviour |
|---|---|
| CN/cluster reconciliation | Read saved Cluster posterior distributions and use their modes/support for CN assignment. Keep conflicting or unsupported assignments explicit |
| Driver representations | Handle protein-change syntax and equivalent/overlapping representations to avoid treating the same biological event as independent events |
| Unknown event ordering | Retain event presence but leave ordering unknown for flat prior-only distributions, unsupported high-copy states or missing CN reconstruction. Presence alone is not timing evidence |
| Sex-specific normal CN | Use validated sex and supported regional normal CN in CCF, multiplicity and CN-event calibration. Retain true X/Y candidates without a blanket chromosome filter; their biological validity still requires QC |
| Negative X/Y estimates | Bound physically impossible estimates at zero and retain original records/reasons in QC. This is an explicit preprocessing convention |
| WGD clock evidence | Include WGD-neutral doubled states as clock candidates, pool matched clonal mutations once per arm/region, and avoid counting two allele labels as independent clocks. This is a change to evidence handling, not merely input formatting |
| WGD without an informative clock | Retain WGD presence with unknown time rather than create quantitative downstream timing from a uniform fallback. Do not force an earlier WGD or agreement with PURPLE |
| Reports and numerical safeguards | Distinguish presence from estimable time, include Yp/Yq axes, bound comparison-graph size, and retain probability-roundoff/solver safeguards. These combine project reporting and numerical-compatibility changes |

The accepted workflow uses `--comparison_method original` for eligible-event
comparisons and native LeagueModel automatic event selection. Seeds,
parallelism and production permutation/season counts are execution settings,
not a new LeagueModel algorithm or a formal between-group significance test.
PURPLE WGD grouping and native Timing WGD calls remain separate definitions.
A complete loss of a sole homolog does not by itself establish pre-WGD
ordering. Unsupported high-copy/ecDNA events do not acquire a formation time
from these changes.

## Source identity, license and validation limits

This branch builds on the accepted **sex-aware V2 snapshot, 2026-10-27**.
Its lineage is `phylyc/PhylogicNDT` commit
`4564634affb7a9ce6f1f2297979f80c501ef874a`, plus inherited and local changes.
That commit alone is not this branch's implementation. Some adaptations
predate this publication; they should not all be attributed as newly developed
here. That historical snapshot's 47 source/reference files are identified by
[`accepted_source_files.json`](../provenance/accepted_source_files.json).

```bash
python3 workflow/check_core.py --phylogic-root .
```

The engine, inherited code and core regression fixtures retain the
[upstream LICENSE](../LICENSE). Independently written adapters, examples,
`tests/test_workflow.py` and new workflow documentation use
[workflow/LICENSE](../workflow/LICENSE); this does not relicense the core.
Original authorship is retained.

Regression tests and synthetic examples check implementation behaviour and
compatibility, not biological accuracy. Exact study reproduction also needs
the matching runtime, permitted inputs, driver/peak definitions and explicit
cohort/exclusion manifests. Patient data, clinical cohort lists and study
results are not supplied. See [source validation](SOURCE_AND_VALIDATION.md)
and the [SASH/PURPLE tutorial](SASH_PURPLE_TO_PHYLOGIC.md).
