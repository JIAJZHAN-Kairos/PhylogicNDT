# PhylogicNDT — SASH/PURPLE hg38 sex-aware repair

This branch, **`sash-purple-sexaware`**, contains the accepted **2026-10-27
sex-aware V2** source snapshot used for Cluster, Timing and LeagueModel, with
the preceding timing repairs included. The `master` branch is unchanged.

The implementation builds on PhylogicNDT and retains its original authorship,
copyright notices and [upstream LICENSE](LICENSE). This fork does not
relicense PhylogicNDT as MIT or claim endorsement by the original authors.
The independently maintained input adapters have their own repository and
license.

## Download and verify this version

```bash
git clone --branch sash-purple-sexaware \
  https://github.com/JIAJZHAN-Kairos/PhylogicNDT.git
cd PhylogicNDT
python3 provenance/verify_source.py --phylogic-root .
```

The verification checks **47 source/reference file SHA-256 values** against
the accepted analysis snapshot. The inference files are byte-identical to
that snapshot; they have not been rewritten for publication.

Install dependencies in a Python 3.12 environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python PhylogicNDT.py --help
```

## Start from SASH somatic VCF + PURPLE results

**[Full step-by-step tutorial →](https://github.com/JIAJZHAN-Kairos/sash-purple-phylogic-workflow/blob/main/docs/SASH_PURPLE_TO_PHYLOGIC.md)**

The companion repository provides:

- Explicit specimen-to-VCF/PURPLE mapping and exact VCF tumour-column selection.
- Ref/alternate read counts and local minor/major CN in MAF-like inputs.
- PURPLE purity/ploidy/sex metadata, the allelic-capseg bridge and one SIF per unrelated specimen.
- Cluster → Timing commands, including the saved Cluster posterior and `--comparison_method original`.
- Group sample-list validation, comparison-table aggregation and native LeagueModel settings.
- Fabricated male/female inputs and a runnable integration smoke test.

```bash
export PHYLOGIC_ROOT="$PWD"
cd ..
git clone https://github.com/JIAJZHAN-Kairos/sash-purple-phylogic-workflow.git
cd sash-purple-phylogic-workflow
export PYTHON="$PHYLOGIC_ROOT/.venv/bin/python"
bash examples/run_demo.sh
```

The demo uses 250 Cluster iterations and original Timing comparisons, followed
by a small LeagueModel test with **2 permutations / 10 seasons**. The tutorial
separately documents the **500 permutations / 200 seasons / 4 workers** cohort
settings, with native automatic event selection and no blanket X/Y filter.

## Main repairs and interpretation

- hg38 chromosome normalisation and preservation of local allele-specific CN.
- CN/Cluster posterior reconciliation, driver-representation deduplication,
  and explicit unknown timing for unsupported/presence-only events.
- Sex- and region-specific normal CN in CCF/multiplicity/event calibration;
  true X/Y candidates are retained without calling loss of absent homologs.
- hg38 PAR splitting and explicit non-estimable male PAR status when the
  reference/alignment provenance is unresolved.
- WGD clocks use matched clonal evidence pooled once per region, including
  neutral doubled states; uninformative WGD cannot generate artificial
  downstream timing. PURPLE/native WGD calls remain distinct.
- Native LeagueModel event selection and arm/focal classification are retained.

See [source and validation](docs/SOURCE_AND_VALIDATION.md), the
[repair/limitation reference](https://github.com/JIAJZHAN-Kairos/sash-purple-phylogic-workflow/blob/main/docs/REPAIRS_AND_LIMITATIONS.md),
and the [original fork README](docs/UPSTREAM_README.md).

These are implementation adaptations and repairs to the published method,
not a new independently benchmarked evolutionary-inference algorithm.
Unknown male PAR remains unresolved; a candidate event's presence is not
proof of an estimable time. LeagueModel infers cohort-level relative event
ordering, not a clone tree connecting unrelated patients or ecDNA formation
time. No study patient data, private cohort lists or study results are added
by this repaired branch.
