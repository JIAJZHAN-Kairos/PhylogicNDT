"""Check a separately obtained licensed core against the accepted code snapshot."""
import argparse
import hashlib
import json
from pathlib import Path


def check(root):
    manifest=Path(__file__).resolve().parents[1]/'provenance/accepted_source_files.json'
    records=json.loads(manifest.read_text())
    failures=[]
    for r in records:
        path=root/r['path']
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=r['sha256']:
            failures.append(r['path'])
    if failures:
        raise ValueError('This is not the accepted repaired core. Missing/different files: '+', '.join(failures))
    print(f'Accepted core verified: {len(records)} files. This does not grant redistribution rights.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phylogic-root',type=Path,required=True)
    check(parser.parse_args().phylogic_root)
