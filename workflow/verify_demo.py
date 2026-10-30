"""Verify this workflow's own synthetic smoke outputs (never load untrusted pickles)."""
import argparse
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
import pickle
import platform
import sys

from check_core import check


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f,delimiter='\t'))


def verify(root,core,report):
    check(core)
    sys.path.insert(0,str(core))
    records=[]
    for sample,expected in [('DEMO_MALE',303),('DEMO_FEMALE',304)]:
        p=root/'results'/sample
        muts=read(p/(sample+'.mut_ccfs.txt'))
        clusters=read(p/(sample+'.cluster_ccfs.txt'))
        comps=read(p/(sample+'.comp.tsv'))
        event_qc=read(p/(sample+'.event_qc.tsv'))
        assert len(muts)==expected and clusters and comps and event_qc
        for r in muts:
            assert r['Sample_ID']==sample
            assert all(math.isfinite(float(r[k])) for k in ['Allelic_CN_minor','Allelic_CN_major'])
            posterior=[float(v) for k,v in r.items() if k.startswith('preDP_ccf_') and k.rsplit('_',1)[-1].replace('.','',1).isdigit()]
            assert len(posterior)==101 and all(math.isfinite(v) and v>=0 for v in posterior)
            assert abs(sum(posterior)-1)<1e-4
        for r in clusters:
            posterior=[float(v) for k,v in r.items() if k.startswith('postDP_ccf_') and k.rsplit('_',1)[-1].replace('.','',1).isdigit()]
            assert len(posterior)==101 and all(math.isfinite(v) and v>=0 for v in posterior)
            assert abs(sum(posterior)-1)<1e-4
        for r in comps:
            assert r['samp']==sample
            pvals=[float(r[k]) for k in ['p1->2','p2->1','unknown']]
            assert all(math.isfinite(v) and -1e-8<=v<=1+1e-8 for v in pvals)
            assert abs(sum(pvals)-1)<1e-5
        par_in_muts=any(r['Chromosome'].removeprefix('chr')=='X' and r['Start_position']=='20000' for r in muts)
        assert par_in_muts==(sample=='DEMO_FEMALE')
        sex_qc=read(p/(sample+'.Cluster.sex_cn_qc.tsv'))
        if sample=='DEMO_MALE':assert any(r['reason'].startswith('PAR_reference_unresolved') for r in sex_qc)
        for name in [sample+'.phylogic_report.html',sample+'.phylogic_timing_report.html']:
            assert (p/name).stat().st_size>100
        records.append(dict(sample=sample,mutation_rows=len(muts),clusters=len(clusters),comparison_rows=len(comps),
                            event_qc_rows=len(event_qc),unresolved_sex_cn_qc_rows=len(sex_qc),par_mutation_retained=par_in_muts))
    with (root/'league/SYNTHETIC_DEMOfig_pkl.pkl').open('rb') as f:
        model=pickle.load(f)
    assert model.full_n_samps==2 and model.n_perms==2 and model.num_seasons==10
    assert all(len(v)==2 and all(math.isfinite(x) for x in v) for v in model.log_odds_full_run.values())
    for suffix in ['positions.pdf','log_odds.pdf']:
        assert (root/'league'/('SYNTHETIC_DEMO.'+suffix)).read_bytes().startswith(b'%PDF')
    import numpy,scipy,pandas
    result=dict(pass_validation=True,validated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                dataset='Fabricated example only; no patient data',python=platform.python_version(),
                platform=platform.system()+' '+platform.machine(),
                runtime=dict(numpy=numpy.__version__,scipy=scipy.__version__,pandas=pandas.__version__),
                core_manifest_sha256=hashlib.sha256((Path(__file__).resolve().parents[1]/'provenance/accepted_source_files.json').read_bytes()).hexdigest(),
                cluster_iterations=250,comparison_method='original',samples=records,
                league=dict(samples=model.full_n_samps,permutations=model.n_perms,seasons=model.num_seasons,selected_events=len(model.final_event_list)),
                limitation='Smoke compatibility only, not a biological benchmark or a full production LeagueModel run.')
    report.parent.mkdir(parents=True,exist_ok=True)
    report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--demo-root',type=Path,required=True)
    parser.add_argument('--phylogic-root',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    a=parser.parse_args()
    verify(a.demo_root,a.phylogic_root.resolve(),a.report)
