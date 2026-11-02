# From our SASH VCF and PURPLE results to PhylogicNDT

Workflow: **hg38 / GRCh38, one selected tumour specimen per individual**,
on branch `sash-purple-sexaware`.

See [change categories](REPAIRS_AND_LIMITATIONS.md) for upstream bug fixes,
SASH/PURPLE and hg38 adaptations, and project-specific handling.

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

Example manifest schema (tab-separated):

| Column | Meaning |
|---|---|
| `sample_id` | Unique output specimen ID; safe letters/numbers/underscore/dot/hyphen |
| `vcf_path` | Path to that specimen's SASH somatic VCF.gz |
| `cn_path` | Corresponding PURPLE somatic CN table |
| `purity_path` | Corresponding PURPLE purity table |
| `vcf_sample` | Optional if exactly equal to `sample_id`; otherwise mandatory exact tumour VCF column |
| `sex` | Optional validated MALE/FEMALE value; must agree with PURPLE if PURPLE has a recognised value |
| `timepoint` | Optional numeric value; default 0 |

See the generated `demo/raw/manifest.tsv` for an actual tab-separated file.
Relative file paths are resolved from the manifest directory. Select one
QC-passed specimen per individual and use the exact tumour VCF column name.

The parser checks `germlineAberrations` if present and rejects non-standard
reported germline/sex-chromosome states for individual review. Validated
MALE/FEMALE metadata is required.

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

Use a **new or empty output directory**. After preparation, check
`preparation_qc.tsv`, `samples.txt` and `input_provenance.json`.

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

Protein changes come from PCGR annotations; indel/type labels are generic.
Check annotation coverage in `preparation_qc.tsv`.

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

This is recorded as `sigma_method=proxy_support`, a support-based proxy
rather than a measured PURPLE confidence interval. Missing support counts
default to zero. Preserve the PURPLE bridge coordinates without a BED shift.

### Container path mapping

For a container with host input mounted at `/input` and results at `/results`,
generate runtime rather than host SIF paths:

```bash
python3 workflow/prepare_sash_purple.py \
  --manifest /host/selected_samples.tsv --outdir /host/new-run/input \
  --input-root /input --results-root /results
```

Match these runtime paths to the container bind mounts.

## 3. Verify the source and install dependencies

From the repository root:

```bash
export CORE="$REPO"
python3 workflow/check_core.py --phylogic-root "$CORE"
```

The checker verifies 47 source/reference file hashes. The
[runtime reference](../provenance/core-runtime.freeze.txt) records the Linux
analysis environment.

Install the pinned inference dependencies, for example:

```bash
python3.12 -m venv /path/to/venv-phylogic
source /path/to/venv-phylogic/bin/activate
python -m pip install -r "$REPO/requirements.txt"
```

## 4. Run Cluster, separately for each unrelated specimen

The examples use Python 3.12, host SIF paths from step 2 and fixed seeds.

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
for example `chr8 127188796 127328999 AMP`. For arm-level CN and driver
mutations only, omit `--cn_peaks` from **both** Cluster and Timing.

Check exit status and log, then verify nonempty `SAMPLE.mut_ccfs.txt`,
`SAMPLE.cluster_ccfs.txt`, CN annotations with finite local CN, and
`SAMPLE.Cluster.sex_cn_qc.tsv`. The generated SIF must have the correct
purity/ploidy, sex and specimen paths.

Clustering uses the supplied mutation input; the driver list annotates events.
Apply any mutation downsampling or hypermutation exclusion separately.

## 5. Run Timing using the same source and Cluster outputs

Still in that specimen's result directory:

```bash
python3 "$CORE/PhylogicNDT.py" Timing \
  -i "$SAMPLE" -sif "$WORK/input/timing_sif/$SAMPLE.sif" \
  --comparison_method original --genome_build hg38 --seed 20260813 \
  --driver_genes_file "$DRIVER_LIST" --cn_peaks "$PEAKS" \
  > timing.log 2>&1
```

Timing reads the **joint Cluster posterior file** via `cluster_ccf_fn`.
Outputs include `SAMPLE.timing.tsv`, `SAMPLE.comp.tsv`,
`SAMPLE.event_qc.tsv`, CN/driver provenance QC and Timing HTML.

Check event presence and timing status in `SAMPLE.event_qc.tsv`.

## 6. Make explicit sample lists for each group

Create one plain-text file per group, one selected specimen ID per line.
Use the current study's sample/donor manifest, with repeated patients and
missing group annotations resolved.

Apply the chosen hypermutation exclusion **before** making these lists.
Intersect the exclusion mask with the current cohort and record before/after
counts and excluded IDs.

For ecDNA groups, use the selected specimen's native AmpliconClassifier
`Classification`: at least one ecDNA feature means positive under our
definition. Explicit non-ecDNA/NA classifier records may be negative;
**a missing classifier result must not be converted to negative**. Record
group counts and cohort/sex/WGD composition.

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
comparisons represent unknown ordering. An adjacent `.audit.json` records
source hashes and counts. After a failed aggregation, correct the input and
rerun to a fresh output path before LeagueModel.

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

The command uses native automatic event selection, including eligible X/Y
candidates.

Accept a production run only after successful exit **and** output checks:
exact sample count; saved model `n_perms=500`, `num_seasons=200`; complete
finite score arrays of length 500; nonempty PDFs/tables; correctly classified
arm/focal labels; frozen input/source hashes; explicit evidence counts for
selected events. Use the native event label map to reconcile comparison
labels.

Expected outputs include `GROUP.positions.pdf`, `GROUP.log_odds.pdf`,
`GROUP.full_prevalence.tsv`, event maps and `GROUPfig_pkl.pkl`.
