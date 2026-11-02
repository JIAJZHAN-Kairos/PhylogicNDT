# PhylogicNDT — from SASH/PURPLE to clonal-event timing

This branch integrates PhylogicNDT with SASH/PURPLE input preparation,
synthetic examples and regression tests.

**[Full tutorial: SASH VCF + PURPLE → Cluster → Timing → LeagueModel](docs/SASH_PURPLE_TO_PHYLOGIC.md)**

## 1. Clone and install

```bash
git clone --branch sash-purple-sexaware \
  https://github.com/JIAJZHAN-Kairos/PhylogicNDT.git
cd PhylogicNDT
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python provenance/verify_source.py --phylogic-root .
```

The hash check confirms 47 source/reference files match the accepted
2026-10-27 sex-aware V2 snapshot.
Inference files are unchanged from that accepted version. The existing
`master` branch is retained.

## 2. Prepare inputs from real SASH/PURPLE files

For each selected tumour specimen, provide:

- A SASH somatic `VCF.gz` with tumour `AD` counts; use the exact tumour sample
  name from the VCF header.
- Its PURPLE `*.purple.cnv.somatic.tsv`, with major/minor allele copy number.
- Its PURPLE `*.purple.purity.tsv`, with purity, ploidy and validated sex.

Create a tab-separated manifest with these columns:

| Column | Example |
|---|---|
| `sample_id` | `TUMOR_A` |
| `vcf_path` | `/data/A.somatic.vcf.gz` |
| `cn_path` | `/data/A.purple.cnv.somatic.tsv` |
| `purity_path` | `/data/A.purple.purity.tsv` |
| `vcf_sample` | Exact VCF tumour column, e.g. `A_TUMOR` |

```bash
export REPO="$PWD"
export CORE="$REPO"
export WORK=/absolute/path/to/new-analysis
python workflow/prepare_sash_purple.py \
  --manifest /absolute/path/to/selected_samples.tsv \
  --outdir "$WORK/input" \
  --results-root "$WORK/results"
```

This generates MAF-like mutation tables with ref/alt counts and local
minor/major CN, the allelic-capseg segment bridge, **one Cluster SIF and one
Timing SIF per specimen**, QC and input SHA-256 records. Different patients
are not combined into a joint-patient SIF. The output directory must be new
or empty.

## 3. Run Cluster, Timing and group-level LeagueModel

Follow the [complete commands in the tutorial](docs/SASH_PURPLE_TO_PHYLOGIC.md):

```text
SASH somatic VCF + PURPLE CN / purity / ploidy / sex
  → workflow/prepare_sash_purple.py
  → Cluster: mutation and cluster CCF posteriors
  → Timing: saved Cluster posterior + original event comparisons
  → explicit group sample lists + workflow/prepare_comparisons.py
  → LeagueModel: cohort-level relative event ordering
```

The tutorial covers gene annotation, driver/focal-peak files, container path
mapping, per-sample QC, hypermutation exclusions, group denominators,
preflight and production settings. Production LeagueModel uses 500
permutations / 200 seasons / 4 workers with native automatic event selection.

## Synthetic example and tests

All synthetic input generation is included; no patient data is needed:

```bash
python -m unittest discover -s tests -v
bash examples/run_demo.sh
```

The demo prepares fabricated male/female SASH/PURPLE files, runs Cluster
(250 iterations), Timing (`--comparison_method original`) and a small native
LeagueModel test (**2 permutations / 10 seasons**), then verifies outputs.
The test suite contains **22 core regression tests and 14 input/comparison
adapter tests**. GitHub Actions runs both and the synthetic integration test
on Linux / Python 3.12. See [validation](docs/VALIDATION.md).

## Changes and attribution

| Category | Scope |
|---|---|
| Confirmed upstream code bug | `clean_local_cn()` tested `not np.nan` instead of the supplied CN, replacing valid values with NaN when called |
| SASH/PURPLE and hg38 adaptations | Read tumour AD and PCGR annotations, attach PURPLE allele CN, normalise chromosome labels, preserve local CN between stages, and pass actual purity/ploidy through the input path |
| Project extensions and conservative handling | Use saved Cluster posteriors for CN reconciliation, handle equivalent driver representations and WGD clock evidence, apply sex-specific normal CN, and leave unsupported event ordering unknown |

The last two categories are not blanket claims of upstream bugs. Some are
inherited adaptations, and others are local analysis choices that can affect
CCFs, eligible events and timing. The accepted workflow retains the original
eligible-event comparison method and native LeagueModel automatic selection.

- [Change categories, evidence and limitations](docs/REPAIRS_AND_LIMITATIONS.md)
- [Source identity and regression validation](docs/SOURCE_AND_VALIDATION.md)
- [Original fork README](docs/UPSTREAM_README.md)

Original PhylogicNDT authorship and [upstream LICENSE](LICENSE) are retained
for the core, inherited code and core regression fixtures. The independently
written `workflow/` and `examples/` scripts, `tests/test_workflow.py` and new
SASH/PURPLE workflow documentation use [workflow/LICENSE](workflow/LICENSE).
That MIT license does not relicense the core. No upstream endorsement is
claimed. Study patient data, private cohort lists and study results are not
added by this branch.
