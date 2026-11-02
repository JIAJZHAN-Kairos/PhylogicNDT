# From our SASH VCF and PURPLE results to PhylogicNDT

This describes the **hg38 / GRCh38, one selected tumour specimen per individual**
workflow. Multiple samples from the same patient require a deliberately built
joint-patient SIF and validation; this adapter does not construct that design.
Never pass an entire cohort as one patient's joint SIF.

Use branch `sash-purple-sexaware` for the commands and source hash check below.
See [change categories and limitations](REPAIRS_AND_LIMITATIONS.md) for the
distinction between upstream bugs, input adaptations and project choices.

## 1. Select specimens and make an explicit input manifest

For every selected specimen, collect these files from the same SASH/PURPLE
run and the same genome build:

| Input | Required content |
|---|---|
| SASH somatic VCF, gzip-compressed | `FILTER=PASS`; tumour `FORMAT/AD` ref/alt counts; exact tumour sample name in `#CHROM` |
| PURPLE `*.purple.cnv.somatic.tsv` | `chromosome`, `start`, `end`, `minorAlleleCopyNumber`, `majorAlleleCopyNumber`; preferably `bafCount`, `depthWindowCount` |
| PURPLE `*.purple.purity.tsv` | Exactly one row with `purity`, `ploidy`, validated `gender` (`MALE` or `FEMALE`) |
| Driver list | One gene symbol per line, appropriate to the intended analysis |
| Focal CN regions | Optional, build-matched GISTIC-like candidate peaks; use the same file in Cluster and Timing |

Use the allele-specific major/minor CN estimates, **not just total CN**.
The purity is a fraction in `(0, 1]`, not a percentage. Ploidy must be positive.
This workflow consumes existing calls; it does not rerun SASH or PURPLE.

Example manifest schema (tab-separated):

| Column | Meaning |
|---|---|
| `sample_id` | Unique output specimen ID; safe letters/numbers/underscore/dot/hyphen |
| `vcf_path` | Path to that specimen's SASH somatic VCF.gz |
| `cn_path` | Corresponding PURPLE somatic CN table |
| `purity_path` | Corresponding PURPLE purity table |
| `vcf_sample` | Optional if exactly equal to `sample_id`; otherwise mandatory exact tumour VCF column |
| `sex` | Optional validated MALE/FEMALE value; must agree with PURPLE if PURPLE has a recognised value |
| `timepoint` | Optional numeric value; default 0, not an inferred evolutionary time |

See the generated `demo/raw/manifest.tsv` for an actual tab-separated file.
Relative file paths are resolved from the manifest directory. The script
checks files and selected VCF names, but does not prove biological donor
identity or cross-run compatibility. Perform those checks and upstream QC
before adding a specimen. Upstream QC acceptance is **not inferred** from
the existence of these files.

The parser checks `germlineAberrations` if present and rejects non-standard
reported germline/sex-chromosome states for individual review. Missing or
uncertain sex must be resolved rather than guessed. The two-normal-sex
baseline is not a validated model for every karyotype.

## 2. Convert SASH/PURPLE to mutation, segment and SIF files

Run from this repository with Python 3.12:

```bash
export REPO="$PWD"
export WORK=/absolute/path/to/new-analysis
mkdir -p "$WORK"
python3 workflow/prepare_sash_purple.py \
  --manifest /absolute/path/to/selected_samples.tsv \
  --outdir "$WORK/input" \
  --results-root "$WORK/results"
```

Use a **new or empty output directory**. For a rerun, select another directory
so the previous inputs and results remain intact. An interrupted preparation
is not accepted; a successful preparation ends with `preparation_qc.tsv`,
`samples.txt` and `input_provenance.json`. Check those before inference.

The adapter streams the VCF and:

1. Selects the explicitly named tumour column and retains PASS sequence
   variants with usable `AD` counts.
2. Splits multiallelic records into one row per alternate allele, using its
   corresponding `AD` entry. The reference count remains `AD[0]`.
3. Maps each locus to the overlapping PURPLE segment and attaches local
   minor/major allele CN. It omits sites with no usable CN segment, total
   allele CN ≤ 0, malformed/missing counts, zero depth or symbolic alleles.
4. Extracts gene symbols and protein changes from the SASH/PCGR
   `INFO/PCGR_CSQ` annotation. Standard VEP `CSQ` alone is **not parsed** by
   this adapter. Without compatible annotation, the variant is kept with
   `Hugo_Symbol=unknown` and loses gene-specific driver interpretation.
5. Writes a MAF-like table for `calc_ccf` plus an allelic-capseg bridge and
   independent per-specimen SIFs.

The MAF-like table is an input adapter, not a complete functional-annotation
product. Generic indel/type labels do not establish coding consequence;
protein changes are taken from PCGR annotations. Review annotation coverage
in the QC table before interpreting driver events.

Outputs:

```text
input/
  maf_with_cn/SAMPLE.withCN.maf
  segments/SAMPLE.alleliccapseg.tsv
  cluster_sif/SAMPLE.sif
  timing_sif/SAMPLE.sif
  preparation_qc.tsv
  samples.txt
  input_provenance.json
```

Mutation columns include `Hugo_Symbol`, chromosome/position/ref/alt,
`Protein_Change`, `t_ref_count`, `t_alt_count`, `local_cn_a1` (minor), and
`local_cn_a2` (major).

The first six Cluster SIF columns must remain in this order, because the
accepted core's Cluster reader reads them positionally:

```text
sample_id  maf_fn  seg_fn  purity  timepoint  ploidy
```

The adapter adds the metadata columns used by the core. Keep the generated
SIF intact.

Timing SIF adds `cluster_ccf_fn` and points `maf_fn` to Cluster's
`SAMPLE.mut_ccfs.txt`. Both are generated in advance with the expected result
paths; **Timing must run only after Cluster succeeds**.

### CN encoding and uncertainty

The mutation table stores PURPLE **absolute** allele CN. The segment table
stores the historical allelic-capseg bridge expected by this core:

```text
mu(C) = 2 * ((1 - purity) + purity * C) / ploidy
C     = 1 + (mu * ploidy / 2 - 1) / purity
```

Keep this inverse pair intact, including on X/Y. This core applies
sex- and region-specific normal baselines **after decoding**. Re-encoding
male X/Y using another normal baseline would double-correct the input.

PURPLE's CN TSV does not provide the required measured allele-CN uncertainty.
The adapter preserves the established support-based proxy:

```text
baf_se = 0.5 / sqrt(max(1, bafCount))
depth_se = 0.5 / sqrt(max(1, depthWindowCount))
cn_sigma = max(0.05, baf_se + depth_se)
sigma_mu = max(0.01, 2 * purity * cn_sigma / ploidy)
```

This is explicitly recorded as `sigma_method=proxy_support`. It is not a
measured PURPLE confidence interval. Support counts default to zero if absent,
giving a broader proxy. Keep this assumption in Methods and sensitivity
interpretation. Coordinates are copied from the existing PURPLE bridge
contract; do not add an ad hoc BED shift to these segment files.

### Container path mapping

For a container with host input mounted at `/input` and results at `/results`,
generate runtime rather than host SIF paths:

```bash
python3 workflow/prepare_sash_purple.py \
  --manifest /host/selected_samples.tsv --outdir /host/new-run/input \
  --input-root /input --results-root /results
```

These runtime paths must exactly match the later binds. Do not mix host
absolute paths and container-only paths in the same SIF.

## 3. Verify the source and install dependencies

From the repository root:

```bash
export CORE="$REPO"
python3 workflow/check_core.py --phylogic-root "$CORE"
```

The checker requires the exact 47 source/reference file hashes from the
accepted sex-aware V2 snapshot. An unpatched upstream clone will fail this check.
Version identity does not itself grant redistribution rights or validate a
new dataset. The [core dependency reference](../provenance/core-runtime.freeze.txt)
records the existing Linux analysis container. The adapter scripts alone
use the standard library; inference requires the pinned dependencies.

Install the pinned inference dependencies, for example:

```bash
python3.12 -m venv /path/to/venv-phylogic
source /path/to/venv-phylogic/bin/activate
python -m pip install -r "$REPO/requirements.txt"
```

Run this repository's core directly by file path; the commands above
install its Python dependencies, not another engine version. The observed Linux container itself is not
distributed by this repository. Package/platform changes require the same smoke tests;
the recorded local macOS test environment is listed separately in the
validation report.

## 4. Run Cluster, separately for each unrelated specimen

The following examples use native Python and host SIF paths from step 2.
Use a Python 3.12 environment containing the core's dependencies. `PYTHONHASHSEED`
and the application seed make the intended seed settings explicit; platform
and dependency differences can still affect floating-point results.

```bash
export MPLBACKEND=Agg PYTHONHASHSEED=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export SAMPLE=TUMOR_A
export DRIVER_LIST=/absolute/path/to/driver_genes.txt
export PEAKS=/absolute/path/to/hg38_focal_peaks.bed
mkdir -p "$WORK/results/$SAMPLE"
cd "$WORK/results/$SAMPLE"
python3 "$CORE/PhylogicNDT.py" Cluster \
  -i "$SAMPLE" -sif "$WORK/input/cluster_sif/$SAMPLE.sif" \
  --maf_input_type calc_ccf --seg_input_type alleliccapseg \
  --use_indels --n_iter 250 --genome_build hg38 --seed 20260813 \
  --driver_genes_file "$DRIVER_LIST" --cn_peaks "$PEAKS" \
  > cluster.log 2>&1
```

The focal-region input is a whitespace-delimited GISTIC-like candidate list,
for example `chr8 127188796 127328999 AMP`. Do not treat a candidate list from
another cohort as evidence of recurrence in the new cohort. If analysing only
arm-level CN and driver mutations, omit `--cn_peaks` from **both** stages.

Check exit status and log, then verify nonempty `SAMPLE.mut_ccfs.txt`,
`SAMPLE.cluster_ccfs.txt`, CN annotations with finite local CN, and
`SAMPLE.Cluster.sex_cn_qc.tsv`. The generated SIF must have the correct
purity/ploidy, sex and specimen paths. A rendered HTML alone does not certify
a valid analysis.

The driver list annotates events; it does not mean clustering is restricted
to those genes. Use the full intended mutation input and document any
downsampling separately. This public adapter performs no 10k mutation
downsampling and no automatic hypermutation filtering.

## 5. Run Timing using the same source and Cluster outputs

Still in that specimen's result directory:

```bash
python3 "$CORE/PhylogicNDT.py" Timing \
  -i "$SAMPLE" -sif "$WORK/input/timing_sif/$SAMPLE.sif" \
  --comparison_method original --genome_build hg38 --seed 20260813 \
  --driver_genes_file "$DRIVER_LIST" --cn_peaks "$PEAKS" \
  > timing.log 2>&1
```

This core requires the **joint Cluster posterior file** via
`cluster_ccf_fn`; do not reconstruct this from a single mean CCF or mix source
versions. Outputs include `SAMPLE.timing.tsv`, `SAMPLE.comp.tsv`,
`SAMPLE.event_qc.tsv`, CN/driver provenance QC and Timing HTML.

Some losses/high-CN events legitimately have no estimable time. Inspect event
status and evidence rather than interpreting a blank timing panel as proof
of a failed run. See [conservative handling and limitations](REPAIRS_AND_LIMITATIONS.md#3-project-extensions-and-conservative-handling).

## 6. Make explicit sample lists for each group

Create one plain-text file per group, one selected specimen ID per line.
Use the current study's fixed sample/donor manifest, not a historical
patient-summary table. Resolve repeated-patient membership and missing group
annotations first.

Apply the chosen hypermutation exclusion **before** making these lists.
Intersect the exclusion mask with the current analysis universe; an old mask
may include specimens already outside the cohort. Preserve the before/after
counts and excluded IDs in your local audit. The script does not choose a
TMB cutoff or infer group membership for you.

For ecDNA groups, use the selected specimen's native AmpliconClassifier
`Classification`: at least one ecDNA feature means positive under our
definition. Explicit non-ecDNA/NA classifier records may be negative;
**a missing classifier result must not be converted to negative**. Freeze
and report group counts, cohort/sex composition and PURPLE/native WGD
composition. Do not publish identifiable clinical manifests in this repo.

```bash
python3 "$REPO/workflow/prepare_comparisons.py" \
  --samples "$WORK/groups/ecDNA_pos.samples.txt" \
  --results-root "$WORK/results" \
  --output "$WORK/groups/ecDNA_pos.comparisons.tsv"
```

Repeat for the other group(s). The script refuses empty/duplicate sample
lists, missing/header-only outputs, specimen mismatches, duplicate event pairs
and invalid probabilities. It normalises either native header layout to:

```text
sample  event1  event2  p_event1_win  p_event2_win  unknown
```

Tabs separate columns. Probabilities sum to one. Fully unknown `(0,0,1)`
comparisons are retained for event presence, without claiming ordering
evidence. An adjacent `.audit.json` records source hashes and counts. A failed
aggregation can leave a partial output; fix the input and use a fresh output
path. Do not run LeagueModel on an aggregation that failed.

## 7. Preflight, then production LeagueModel

Run in a new, group-specific output directory. First substitute
`--n_perms 2 --n_seasons 10 --n_jobs 1 --perm_chunk_size 1` and use a separate
preflight directory. Inspect event maps, probabilities, sample count and
native model parameters before the full run.

```bash
export GROUP=ecDNA_pos
mkdir -p "$WORK/league/$GROUP"
cd "$WORK/league/$GROUP"
python3 "$CORE/PhylogicNDT.py" LeagueModel \
  --cohort "$GROUP" \
  --comparison_fn "$WORK/groups/$GROUP.comparisons.tsv" \
  --final_sample_list "$WORK/groups/$GROUP.samples.txt" \
  --n_perms 500 --n_seasons 200 --n_jobs 4 --perm_chunk_size 8 \
  --percent_subset 0.8 --min_event_prevalence 0.05 --seed 20260813 \
  > league.log 2>&1
```

No explicit event list is supplied: native automatic selection is retained,
including eligible X/Y candidates. Do not silently replace these study
settings with upstream defaults. These bootstrap/permutation model settings
do not constitute a formal between-group significance test.

Accept a production run only after successful exit **and** output checks:
exact sample count; saved model `n_perms=500`, `num_seasons=200`; complete
finite score arrays of length 500; nonempty PDFs/tables; correctly classified
arm/focal labels; frozen input/source hashes; explicit evidence counts for
selected events. Use the native event label map to reconcile comparison
labels. Pickles are unsafe to load from untrusted sources; inspect only your
own generated model files in the matching environment.

Expected outputs include `GROUP.positions.pdf`, `GROUP.log_odds.pdf`,
`GROUP.full_prevalence.tsv`, event maps and `GROUPfig_pkl.pkl`. Keep the large
model and raw inputs on the analysis host; transfer only needed figures,
compact tables and validation records, with SHA-256 checks.

## 8. Interpretation

- Position is a relative season rank within that group's selected event set,
  not molecular time, years or a phylogenetic branch length.
- Stored early log10 odds is positive for earlier events. The native PDF
  plots its **negative**. With 200 seasons, ±log10(200) ≈ ±2.301 is the score
  bound; reaching it does not specify an absolute biological time.
- Check direct informative comparisons. Presence-only/homdel events can
  appear in native plots without an estimable time and must be labelled.
- Groups can select different event sets. Compare shared events and evidence
  alongside cohort/sex/WGD composition; do not compare raw positions as if
  the scales and populations were identical.
- PURPLE WGD grouping and native Timing WGD calls are independent definitions.
  Do not force concordance or force WGD to be early.
- Cohort event ordering is not a clone tree connecting unrelated patients,
  and ecDNA-group differences do not establish ecDNA causation/formation time.
