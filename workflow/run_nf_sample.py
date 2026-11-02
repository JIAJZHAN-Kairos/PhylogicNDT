"""Run one explicitly mapped specimen in an isolated Nextflow task directory."""
import argparse
import base64
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from prepare_sash_purple import prepare


CORE = Path(__file__).resolve().parents[1]


def run(args):
    meta = json.loads(base64.b64decode(args.meta))
    sample = meta['sample_id']
    root = Path(sample + '.run').resolve()
    root.mkdir()
    paths = dict(vcf_path=Path(args.vcf).resolve(), cn_path=Path(args.cn).resolve(),
                 purity_path=Path(args.purity).resolve())
    for field, path in paths.items():
        expected = meta.get(field.replace('_path', '_sha256'))
        if expected and hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError('Input checksum mismatch: ' + field)
    manifest = root/'manifest.tsv'
    row = {k: meta.get(k, '') for k in ['sample_id', 'sex', 'par_mode', 'vcf_sample']}
    row.update({k: str(v) for k, v in paths.items()})
    with manifest.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(row), delimiter='\t')
        writer.writeheader()
        writer.writerow(row)
    results = root/'results'
    prepare(manifest, root/'input', str(results))
    output = results/sample
    output.mkdir(parents=True)
    common = ['-i', sample, '--genome_build', 'hg38', '--seed', str(args.seed),
              '--driver_genes_file', str(Path(args.drivers).resolve()),
              '--cn_peaks', str(Path(args.peaks).resolve())]
    for stage, extra in [
        ('Cluster', ['--maf_input_type','calc_ccf','--seg_input_type','alleliccapseg',
                     '--use_indels','--n_iter',str(args.iterations)]),
        ('Timing', ['--comparison_method','original'])]:
        sif = root/'input'/(stage.lower()+'_sif')/(sample+'.sif')
        with (output/(stage.lower()+'.log')).open('w') as log:
            subprocess.run([sys.executable, str(CORE/'PhylogicNDT.py'), stage,
                            '-sif', str(sif), *common, *extra], cwd=output,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
    for suffix in ['mut_ccfs.txt', 'cluster_ccfs.txt', 'comp.tsv', 'timing.tsv']:
        path = output/(sample+'.'+suffix)
        if not path.is_file() or not path.stat().st_size:
            raise ValueError('Missing output: '+str(path))
    (root/'run_metadata.json').write_text(json.dumps(meta, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['meta','vcf','cn','purity','drivers','peaks']:
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--iterations', type=int, required=True)
    parser.add_argument('--seed', type=int, required=True)
    run(parser.parse_args())
