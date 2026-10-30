"""Aggregate exactly the listed specimens' Timing outputs for one LeagueModel group."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

FIELDS = ['sample','event1','event2','p_event1_win','p_event2_win','unknown']
LEGACY = ['samp','eve1','eve2','p1->2','p2->1','unknown']


def aggregate(samples_path, results_root, output):
    samples = samples_path.read_text().split()
    if not samples or len(samples) != len(set(samples)):
        raise ValueError('Require a nonempty unique explicit sample list')
    if any(Path(s).name != s or s in {'.','..'} for s in samples):
        raise ValueError('Unsafe sample ID')
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    records = []
    with output.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for sample in samples:
            source = results_root/sample/(sample+'.comp.tsv')
            seen = set()
            with source.open(newline='') as f:
                reader = csv.DictReader(f, delimiter='\t')
                if reader.fieldnames not in [FIELDS, LEGACY]:
                    raise ValueError('Unexpected comparison header: '+str(source))
                count = informative = 0
                for original in reader:
                    r = dict(zip(FIELDS, (original[k] for k in reader.fieldnames)))
                    if r['sample'] != sample:
                        raise ValueError('Comparison specimen mismatch: '+str(source))
                    pair = tuple(sorted((r['event1'], r['event2'])))
                    if pair[0] == pair[1] or pair in seen:
                        raise ValueError('Duplicate/self event pair: '+str(source))
                    seen.add(pair)
                    p = [float(r[k]) for k in FIELDS[3:]]
                    if not all(math.isfinite(v) for v in p) or min(p)<-1e-8 or max(p)>1+1e-8 or abs(sum(p)-1)>1e-5:
                        raise ValueError('Invalid comparison probabilities: '+str(source))
                    writer.writerow(r)
                    count += 1
                    informative += p[0]+p[1] > 1e-12
            if count == 0:
                raise ValueError('Header-only comparison file: '+str(source))
            records.append(dict(sample=sample,comparison_rows=count,informative_rows=informative,
                                sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    report=dict(samples=len(samples), sample_list_sha256=hashlib.sha256(samples_path.read_bytes()).hexdigest(),
                output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(), inputs=records,
                note='Unknown comparisons retained for presence; they do not provide ordering evidence.')
    output.with_suffix('.audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(samples=len(samples),rows=sum(r['comparison_rows'] for r in records))))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples',type=Path,required=True)
    parser.add_argument('--results-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    aggregate(args.samples,args.results_root,args.output)
