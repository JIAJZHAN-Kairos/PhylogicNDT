# PhylogicNDT

## Installation 
First:  Clone this repository

    git clone https://github.com/broadinstitute/PhylogicNDT.git
    cd PhylogicNDT

### Docker Install
Install docker from https://www.docker.com/community-edition#/download

	docker build --tag phylogicndt . 

## Using the Package

    ./PhylogicNDT.py --help

If running from the docker, first run:

	docker run -i -t phylogicndt
	cd phylogicndt


## Clustering 

To run clustering on the provided sample input data:

 To specify inputs: 

	./PhylogicNDT.py Cluster -i Patient_ID  -s Sample1_id:Sample1_maf:Sample1_CN_seg:Sample1_Purity:Sample1_Timepoint -s Sample2_id:Sample2_maf:Sample2_CN_seg:Sample2_Purity:Sample2_Timepoint ... SampleN_info 

alternatively - provide a tsv sample_information_file (.sif) 

with headers: sample_id maf_fn seg_fn purity timepoint

    ./PhylogicNDT.py Cluster -i Patient_ID  -sif Patient.sif

the .maf should contain pre-computed raw ccf histograms based on mutations alt/ref count 
(Absolute annotated mafs or .Rdata files are also supported)
if the ccf histograms are absent - the `--maf_input_type` flag must be set to
`calc_ccf` and sample purity must be provided. Also local copy number must
be attached to each mutation in the maf with columns named `local_cn_a1` and `local_cn_a2`

CN_seg is optional to annotate copy-number information on the trees
	
<sub>Acknowledgment: Clustering Module is partially inspired (primary 1D clustering) by earlier work of Carter & Getz (Landau D, Carter S , Stojanov P et al. Cell 152, 714–726, 2013)</sub>

### Clustering Supplement
- **`local_cn_a1`**: Minor allele local copy number  
- **`local_cn_a2`**: Major allele local copy number 
> Required when using `--maf_input_type calc_ccf`.  
> The MAF must also include `t_alt_count` and `t_ref_count` (or equivalent), along with sample purity provided via `-s` or `.sif`.  
> If CCF histograms are already included (e.g., Absolute annotated MAFs or `.Rdata` files), then this flag is not necessary.

### Typical Usage

1. **Directly specify sample information:**
```bash
./PhylogicNDT.py Cluster \
  -i Patient_ID \
  -s Sample1_id:Sample1.maf:Sample1.seg:0.62:1 \
  -s Sample2_id:Sample2.maf:Sample2.seg:0.58:2 \
  -ni 1000
```
2. **Use a .sif file (tab-delimited with headers):**

| sample_id | maf_fn      | seg_fn      | purity | timepoint |
|-----------|-------------|-------------|--------|-----------|
| S1        | Sample1.maf | Sample1.seg | 0.62   | 0         |
| S2        | Sample2.maf | Sample2.seg | 0.58   | 0         |
```bash
./PhylogicNDT.py Cluster -i Patient_ID -sif Patient.sif
```
3. **When MAF does not contain pre-computed CCF histograms:**
```bash
./PhylogicNDT.py Cluster \
  -i Patient_ID \
  -sif Patient.sif \
  -mt calc_ccf 
```
### Full Option Reference — `PhylogicNDT.py Cluster`

| Parameter | Description | Usage / Notes |
|-----------|-------------|---------------|
| `-h, --help` | Show help message and exit | Quick reference |
| `-i, --indiv_id INDIV_ID` | Patient/Case ID | **Required** |
| `-s, --sample SAMPLE_DATA` | Sample data in format `sample_id:maf_fn:seg_fn:purity:timepoint`; repeatable | Use when not using `.sif` |
| `-sif, --sample_information_file SIF` | Tab-separated file with headers `sample_id maf_fn seg_fn purity timepoint` | Recommended for multiple samples |
| `-bl, --artifact_blacklist` | Path to artifact blacklist | Filter recurrent artifacts |
| `-wl, --artifact_whitelist` | Path to artifact whitelist (takes precedence over blacklist) | Use curated real calls |
| `-drv, --driver_genes_file` | Driver list file | Highlight drivers in report |
| `-tr, --treatment_data` | Path to treatment data file | For annotated timelines |
| `-ts, --tumor_size` | Tumor size | Optional metadata |
| `-bt, --blacklist_threshold` | CCF threshold for blacklisting clusters (BuildTree/Cell Population) | Default usually fine |
| `--seed SEED` | Random seed | For reproducibility |
| `--cluster_ccf_trace` | Save MCMC trace for constrained CCF values | Extra diagnostics |
| `--cluster_order CLUSTER_ORDER` | Comma-separated order of clusters | Consistent coloring across trees |
| `--tree_number TREE_NUMBER` | Select tree from ranked posteriors (1-indexed, default=1) | For reproducibility |
| `--select_tree` | Interactively select a tree | Manual inspection |
| `--intersect_cn_trees` | Reconcile CN segments/focal events across samples | Ignores `--cn_peaks` if set |
| `--cn_peaks GISTIC_FN` | Interval file with focal amp/del regions | Use GISTIC-style input |
| `-rb, --run_with_BuildTree` | Run BuildTree right after clustering | One-step workflow |
| `--PoN PON` | Panel of Normals to use (`false` to skip) | For artifact filtering |
| `--Delete_Blacklist` | Generate new blacklist from PoN | Refresh artifact filter |
| `-g, --grid_size GRID_SIZE` | Number of CCF bins | Must align with input grid size |
| `-ni, --n_iter ITER` | Number of iterations | Increase (e.g., 2000–5000) for stability |
| `--use_indels` | Include indels in clustering | Default: added after clustering |
| `--impute` | Assume 0 CCF for missing mutations | For sparse data |
| `-mc, --min_coverage MIN_COV` | Minimum coverage threshold | Variants below are excluded then reassigned |
| `-ct, --cancer_type CANCER_TYPE` | Cancer type string | Affects focal event calling |
| `--Pi_k_r PI_K_R` | *r* parameter of negative binomial prior (#clusters) | Advanced |
| `--Pi_k_mu PI_K_MU` | *mu* parameter of negative binomial prior (#clusters) | Advanced |
| `--order_by_timepoint` | Order samples by `timepoint` field | For longitudinal samples |
| `--maf` | Output clustered MAF | For downstream use |
| `--no_html` | Suppress HTML report | For batch CLI runs |
| `--time_points TIME_POINTS` | Manually supply timepoints | Rarely needed |
| `--scale` | Scaling option (undocumented) | Advanced |
| `-ns, --use_first_n_samples N_SAMPLES` | Use only first N samples (0 = all) | Useful for quick tests |
| `-mt, --maf_input_type MAF_INPUT_TYPE` | Specify MAF input type (`calc_ccf`) | Required if no CCF histograms are present |

**In BrCa project:**
```bash
PhylogicNDT.py Cluster \
    -i "$sample_id" \
    -sif "$indiv_sif" \
    --driver_genes_file "$DRIVER" \ #from intogen BRCA DB
    --cn_peaks "$PEAKS" \ #from Gistc2.0 output
    --maf_input_type calc_ccf \ 
    --use_indels \
    -rb \
    -seg_input_type='timing_format'
```

## BuildTree (and GrowthKinetics) 
The GrowthKinetics module fully incorporates the BuildTree libraries, so when rates are desired, there is no need to run both. 

 - The -w flag should provide a measure of tumor burden, with one value per input sample maf in clustering. **When ommited, stable tumor burden is assumed.** 
  - The -t flag should provide relative time for spacing the samples.
    **When omitted, equal spacing is assumed.** 

Just BuildTree

	./PhylogicNDT.py BuildTree -i Indiv_ID -sif Patient.sif  -m mutation_ccf_file -c cluster_ccf_file 

GrowthKinetics

	./PhylogicNDT.py GrowthKinetics -i Indiv_ID -sif Patient.sif -ab cell_population_abundance_mcmc_trace -w 10 10 10 10 10 -t 1 2 3 4 5 

Run Cluster together with BuildTree

	./PhylogicNDT.py Cluster -i Patient_ID  -sif Patient.sif -rb

## SinglePatientTiming

SinglePatientTiming requires a maf input and a seg file input for each sample.
The maf file should be the output of PhylogicNDT Clustering module (`*.mut_ccfs.txt`).
The seg file should have the following columns:

    Chromosome	Start.bp	End.bp	mu.minor	sigma.minor	mu.major	sigma.major

### SEG file format (for PhylogicNDT)

The segmentation (`seg_fn`) file provided to **PhylogicNDT** should contain the following columns:

| Column       | Meaning in PhylogicNDT | Corresponds to PURPLE output |
|--------------|-------------------------|-------------------------------|
| `Chromosome` | Chromosome number / label | `chromosome` column in PURPLE `*.purple.segment.tsv` |
| `Start.bp`   | Start position (bp) of the segment | `start` |
| `End.bp`     | End position (bp) of the segment | `end` |
| `mu.minor`   | Estimated **mean copy number** of the *minor allele* in this segment | `minorAlleleCopyNumber` |
| `sigma.minor`| **Uncertainty / standard deviation** of the minor allele CN estimate | `minorAlleleCopyNumberDeviation` |
| `mu.major`   | Estimated **mean copy number** of the *major allele* in this segment | `majorAlleleCopyNumber` |
| `sigma.major`| **Uncertainty / standard deviation** of the major allele CN estimate | `majorAlleleCopyNumberDeviation` |

---

To run SinglePatientTiming:

    ./PhylogicNDT.py Timing -i Indiv_ID -sif Patient.sif

**In BrCa project:**

```bash
PhylogicNDT.py Timing \
      -i "$sample_id" \
      -sif "$indiv_sif" \
      --driver_genes_file "$DRIVER" \
      --cn_peaks "$PEAKS"
```

### LeagueModel
LeagueModel requires an input of comparison tables.
The comparison tables should be the output of SinglePatientTiming ending in ".comp.tsv"

To run LeagueModel:

    ./PhylogicNDT.py LeagueModel -cohort Cohort -comps comp1 comp2 ... compN

Alternatively, one can use a single aggregated table. The table should have the following columns:

    sample  event1  event2  p_event1_win    p_event2_win    unknown

To run with the aggregated table:

    ./PhylogicNDT.py LeagueModel -cohort Cohort -comparison_cn comps

**In BrCa project:**

```bash
PhylogicNDT.py LeagueModel \
      -cohort "$cohort" \
      -comps $COMP_FILES \
      -drv "$DRIVER" \
      --cn_peaks "$PEAKS" \
      --seed 123
```
**✅ For detailed output explanation, see `output/README.md`**

## Notes on Code Modifications

This fork contains several modifications compared to the original PhylogicNDT:

1. **Stirling number calculation (`lstirling`)**  
   - Original implementation: used `root_scalar(phiprime, bracket=[0.1, n*m])`, which sometimes failed with  
     error *"f(a) and f(b) must have different signs"* for some `(n,m)` values.  
   - Modification: replaced with a log-domain DP fallback (`logaddexp`) in `safe_log_stirling_coef`.  
   - Impact: more robust for large `(n,m)`; slightly slower but avoids crashes.

2. **Cytogenetic band reference (`cytoBand.txt`)**  
   - Original: `PhylogicNDT/data/supplement_data/cytoBand.txt` was provided in **hg19/GRCh37** format.  
   - Modification: updated to the **hg38/GRCh38** version of cytoband annotations.  
   - Impact: ensures consistency when running analyses on hg38-aligned data (e.g., BrCa project).

3. **Driver gene list integration (IntOGen BRCA)**  
   - Added IntOGen-derived driver gene list specific to **breast cancer (BRCA)**.  
   - Location: `PhylogicNDT/data/supplement_data/Driver_genes_v1.0.txt` (custom addition).  
   - Impact: enables downstream clustering/BuildTree/Timing modules to highlight BRCA driver events consistently in project analyses.

4. **GISTIC CNV peaks (`gistic_cn_peaks.bed`)**  
   - Added a GISTIC focal amplification/deletion peaks file in `PhylogicNDT/data/supplement_data/gistic_cn_peaks.bed`.  
   - Impact: allows the `--cn_peaks` option in Clustering/Timing modules to use pre-defined focal CNV regions, improving consistency of CNV annotation in hg38 analyses.