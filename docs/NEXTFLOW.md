# Nextflow / Seqera

The `sash-purple-sexaware` branch provides a DSL2 workflow: one isolated task
per specimen runs SASH/PURPLE preparation, Cluster and Timing. A cohort task
requires all expected specimens, validates comparison probabilities, and runs
LeagueModel for explicitly assigned groups. Process errors terminate the run;
failed samples are never silently dropped. Infrastructure exit codes can retry.

Use a new output prefix for every independent analysis. Published files use
`overwrite: false`; preserve the work directory to resume interrupted runs.
Use immutable Git commits and container digests in production Launchpad versions.
The image contains public source and dependencies only; input manifests and
patient files belong in private storage.

The CSV samplesheet has `sample_id,vcf_path,cn_path,purity_path,vcf_sample,sex,par_mode,league_group`.
File fields accept local or S3 paths. Optional `vcf_sha256,cn_sha256,purity_sha256`
are checked against staged bytes. IDs must be unique. Set `expected_samples`
to the explicit denominator. Each row is an independent specimen; longitudinal
joint clustering is not performed. The implementation currently supports hg38.
`league_group=EXCLUDED` retains Cluster/Timing results but omits that sample
from LeagueModel. Assign groups and exclusions explicitly in the private manifest.

Use SASH `<sample_pair>/smlv_somatic/filter/<tumor_id>.pass.vcf.gz` directly.
This file retains `PCGR_CSQ`, tumour FORMAT/AD and the named sample columns.
The report-only PCGR VCF omits genotypes and is not an interchangeable input.

```bash
nextflow run JIAJZHAN-Kairos/PhylogicNDT -r sash-purple-sexaware -profile aws \
  --input s3://PRIVATE/manifest.csv --expected_samples 380 \
  --driver_genes s3://PRIVATE/driver_list.txt \
  --cn_peaks s3://PRIVATE/gistic_cn_peaks.bed \
  --container public.ecr.aws/REGISTRY/phylogicndt-sexaware@sha256:DIGEST \
  --outdir s3://PRIVATE/new-results -work-dir s3://PRIVATE/new-work
```

Defaults match the established run settings: Cluster 250 iterations,
Timing `original`, seed 20260813, LeagueModel 500 permutations / 200 seasons,
subset 0.8, minimum prevalence 0.05. LeagueModel uses four workers with chunk
size 10. Hold these settings, container digest and PYTHONHASHSEED fixed for
repeatability. Compare ordering estimates cautiously across changed settings.

Outputs: `samples/SAMPLE.run/` contains preparation QC, provenance, original
Cluster/Timing outputs and logs. `cohort/` contains complete-cohort comparisons,
per-group LeagueModel outputs and a denominator validation receipt.
Generated SIF paths refer to the task filesystem; the prepared tables and
provenance remain available in the published sample directory.

The deployment clips negative PURPLE allele CN estimates to zero at adaptation,
per the cohort owner's input rule. Each affected segment's original and clipped
values are retained in input provenance; QC reports the number clipped. Original
PURPLE files are not edited. A segment with total CN zero contributes no valid
mutation CCF. `provenance/nextflow_source_files_20261102.json` identifies this
version; older repair snapshots and validation receipts remain historical.
