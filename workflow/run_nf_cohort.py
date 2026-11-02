"""Validate the complete specimen denominator and run explicitly assigned groups."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

from prepare_comparisons import aggregate


def run(args):
    roots = sorted(p for p in Path.cwd().glob('*.run') if p.is_dir() and not p.name.startswith('.'))
    records = [json.loads((p/'run_metadata.json').read_text()) for p in roots]
    ids = [r['sample_id'] for r in records]
    if len(ids) != args.expected or len(ids) != len(set(ids)):
        raise ValueError('Incomplete or duplicated cohort; refuse LeagueModel')
    out = Path('cohort')
    out.mkdir()
    results = out/'comparisons'
    results.mkdir()
    groups = {}
    for root, row in zip(roots, records):
        sample = row['sample_id']
        (results/sample).mkdir()
        shutil.copyfile(root/'results'/sample/(sample+'.comp.tsv'), results/sample/(sample+'.comp.tsv'))
        group = row['league_group']
        if group != 'EXCLUDED':
            groups.setdefault(group, []).append(sample)
    all_samples = out/'all_samples.txt'
    all_samples.write_text(''.join(s+'\n' for s in ids))
    aggregate(all_samples, results, out/'all.comparisons.tsv')
    for group, samples in sorted(groups.items()):
        folder = (out/group).resolve()
        folder.mkdir()
        sample_list = folder/'samples.txt'
        sample_list.write_text(''.join(s+'\n' for s in samples))
        comparisons = folder/'comparisons.tsv'
        aggregate(sample_list, results, comparisons)
        with (folder/'league.log').open('w') as log:
            subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1]/'PhylogicNDT.py'),
                'LeagueModel','--cohort',group,'--comparison_fn',str(comparisons),
                '--final_sample_list',str(sample_list),'--n_perms',str(args.perms),
                '--n_seasons',str(args.seasons),'--n_jobs',str(args.jobs),
                '--perm_chunk_size','10','--percent_subset','0.8',
                '--min_event_prevalence','0.05','--seed',str(args.seed)],
                cwd=folder, stdout=log, stderr=subprocess.STDOUT, check=True)
    (out/'cohort_validation.json').write_text(json.dumps(dict(completed_samples=len(ids),
        groups={k:len(v) for k,v in groups.items()},
        excluded_from_league=sum(r['league_group']=='EXCLUDED' for r in records)), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['expected','perms','seasons','jobs','seed']:
        parser.add_argument('--'+name, type=int, required=True)
    run(parser.parse_args())
