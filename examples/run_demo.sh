#!/usr/bin/env bash
# Run fabricated SASH/PURPLE inputs through this repository's repaired core.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
PHYLOGIC_ROOT=${PHYLOGIC_ROOT:-"${repo}"}
python_bin=${PYTHON:-python3}
core=$(cd "${PHYLOGIC_ROOT}" && pwd)
demo=${DEMO_ROOT:-"${repo}/demo"}
mkdir -p "${demo}"
demo=$(cd "${demo}" && pwd)
"${python_bin}" "${repo}/workflow/check_core.py" --phylogic-root "${core}"
"${python_bin}" "${repo}/examples/make_demo.py" --outdir "${demo}/raw"
"${python_bin}" "${repo}/workflow/prepare_sash_purple.py" \
  --manifest "${demo}/raw/manifest.tsv" --outdir "${demo}/input" --results-root "${demo}/results"
export MPLBACKEND=Agg PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for sample in DEMO_MALE DEMO_FEMALE; do
  mkdir -p "${demo}/results/${sample}"
  (
    cd "${demo}/results/${sample}"
    "${python_bin}" "${core}/PhylogicNDT.py" Cluster -i "${sample}" \
      -sif "${demo}/input/cluster_sif/${sample}.sif" \
      --maf_input_type calc_ccf --seg_input_type alleliccapseg \
      --use_indels --n_iter 250 --genome_build hg38 --seed 20260813 \
      --driver_genes_file "${demo}/raw/driver_genes.txt" --cn_peaks "${demo}/raw/focal_peaks.bed" > cluster.log 2>&1
    "${python_bin}" "${core}/PhylogicNDT.py" Timing -i "${sample}" \
      -sif "${demo}/input/timing_sif/${sample}.sif" \
      --comparison_method original --genome_build hg38 --seed 20260813 \
      --driver_genes_file "${demo}/raw/driver_genes.txt" --cn_peaks "${demo}/raw/focal_peaks.bed" > timing.log 2>&1
  )
done
"${python_bin}" "${repo}/workflow/prepare_comparisons.py" \
  --samples "${demo}/input/samples.txt" --results-root "${demo}/results" --output "${demo}/demo.comparisons.tsv"
mkdir -p "${demo}/league"
(
  cd "${demo}/league"
  "${python_bin}" "${core}/PhylogicNDT.py" LeagueModel --cohort SYNTHETIC_DEMO \
    --comparison_fn "${demo}/demo.comparisons.tsv" --final_sample_list "${demo}/input/samples.txt" \
    --n_perms 2 --n_seasons 10 --n_jobs 1 --perm_chunk_size 1 \
    --percent_subset 0.8 --min_event_prevalence 0.05 --seed 20260813 > league.log 2>&1
)
"${python_bin}" "${repo}/workflow/verify_demo.py" --demo-root "${demo}" \
  --phylogic-root "${core}" --report "${demo}/smoke_validation.json"
echo "Synthetic smoke run finished: ${demo}. LeagueModel uses 2/10 for testing, not production."
