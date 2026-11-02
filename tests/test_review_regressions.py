"""Regression cases from the 2026-11-02 SASH/PURPLE integration review."""
import csv
import gzip
import io
import itertools
import os
from pathlib import Path
import pickle
import random
import subprocess
import sys
import tempfile
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from data.Enums import Genome
from data.Patient import Patient
from data.Sample import TumorSample
from SinglePatientTiming.TimingEngine import TimingEngine
from SinglePatientTiming.TimingRepair import eligible
from SinglePatientTiming.SinglePatientTiming import compare_events_original
from utils.calc_ccf import calc_ccf
from output.PhylogicOutput import PhylogicOutput
from workflow.prepare_sash_purple import (MAF_COLUMNS, alleliccapseg_values, parse_vcf,
                                         prepare, read_tsv, write_tsv)


def cn_sample(rows, sex='FEMALE'):
    sample = TumorSample.__new__(TumorSample)
    sample.sample_name, sample.sex = 'S', sex
    sample.sex_cn_unresolved = []
    sample.purity = .6
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'cn.tsv'
        with path.open('w') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['Chromosome','Start.bp','End.bp','mu.minor','sigma.minor','mu.major','sigma.major'])
            for chrom, start, end, minor, major in rows:
                writer.writerow([chrom, start, end,
                    *alleliccapseg_values(minor, .6, 2, 100, 100),
                    *alleliccapseg_values(major, .6, 2, 100, 100)])
        sample.CnProfile = sample._resolve_CnEvents(str(path), 'alleliccapseg', .6, 2, Genome('hg38').CHROMS)
    sample.mutations, sample.low_coverage_mutations = [], {}
    sample._mut_varstring_hashtable = {}
    return sample


def patient(sample):
    return NS(sample_list=[sample], genome=Genome('hg38'), driver_genes={'KRAS'},
              saved_cluster_ccfs={1:np.eye(101)[[100]], 2:np.eye(101)[[40]]})


def focal_engine():
    sample = cn_sample([('1', 1000000, 2000000, 1, 2)])
    for i in range(6):
        pos = 1000100 + i * 100
        sample.mutations.append(NS(chrN='chr1', pos=pos, type='SNP', var_str=f'chr1:{pos}:A:T',
            gene='KRAS', alt='T', ref='A', alt_cnt=55 if i<3 else 27, ref_cnt=45 if i<3 else 73,
            ccf_1d=np.eye(101)[100], prot_change=f'p.Gly{i+1}Asp', cluster_assignment=1,
            local_cn_a1=1., local_cn_a2=2.))
    sample.low_coverage_mutations['focal'] = NS(type='CNV', var_str='focal', cn_category='Focal_gain',
        chrN='chr1', start=1000000, end=2000000, arm='p', a1=False, ccf_hat=1.,
        local_cn_a1=1., local_cn_a2=2.)
    return TimingEngine(patient(sample))


class TimingRegressions(unittest.TestCase):
    def test_arm_subclonal_ccf_survives_discretization(self):
        for chrom, sex, minor, major, name in [('1','FEMALE',1,1.4,'gain_chr1p'),
                ('1','FEMALE',.6,1,'loss_chr1p'), ('X','MALE',0,1.4,'gain_chrXp')]:
            with self.subTest(chrom=chrom, major=major):
                end = Genome('hg38').CENT_DICT['chr'+chrom]
                sample = cn_sample([(chrom, 1, end, minor, major)], sex)
                engine = TimingEngine(patient(sample))
                event = engine.all_cn_events[name][0]
                self.assertEqual(event.cluster_id, 2)
                self.assertFalse(event.is_clonal)
                engine.time_events()
                self.assertEqual(event.timing_status, 'subclonal')

    def test_arm_ccf_is_length_weighted_and_rebased_for_wgd(self):
        sample = cn_sample([('1',1,40000001,1,1.2),('1',40000001,123400000,1,1.8)])
        state = TimingEngine(patient(sample)).sample_list[0].cn_states['chr1p']
        expected = (.2*40000000 + .8*83399999)/123399999
        self.assertAlmostEqual(state.cn_events[0].ccf_hat, expected)
        sample = cn_sample([('1',1,123400000,2,2.4)])
        state = TimingEngine(patient(sample)).sample_list[0].cn_states['chr1p']
        state.call_events(wgd=True)
        self.assertEqual(len(state.cn_events), 1)
        self.assertAlmostEqual(state.cn_events[0].ccf_hat, .4)

    def test_gain_intervals_use_same_scale_as_point_estimate(self):
        for chrom, sex, minor, major in [('1','FEMALE',1,1.4), ('X','MALE',0,1.4),
                                        ('Y','MALE',0,1.4), ('1','FEMALE',1,3.4)]:
            sample = cn_sample([(chrom,3000000,4000000,minor,major)], sex)
            data = next(iter(sample.CnProfile['chr'+chrom])).data[1]
            scale = data['cn_a2'] - 1
            self.assertAlmostEqual(data['ccf_low_a2'], (major-1-.1)/scale)
            self.assertAlmostEqual(data['ccf_high_a2'], min(1.,(major-1+.1)/scale))
            self.assertLessEqual(data['ccf_low_a2'], data['ccf_hat_a2'])
            self.assertGreaterEqual(data['ccf_high_a2'], data['ccf_hat_a2'])

    def test_identical_fragmented_gains_form_one_arm_event(self):
        sample = cn_sample([('1',a,b,1,1.4) for a,b in
                            [(1,40000000),(40000000,80000000),(80000000,123400000)]])
        data = Patient(indiv_name='S',genome_build='hg38')
        data.sample_list = [sample]
        data.get_arm_level_cn_events()
        self.assertEqual(len(sample.low_coverage_mutations), 1)
        self.assertAlmostEqual(next(iter(sample.low_coverage_mutations.values())).ccf_hat, .4)

    def test_focal_mutation_posterior_reaches_output_and_comparisons(self):
        engine = focal_engine()
        random.seed(17); np.random.seed(17)
        engine.time_events()
        event = next(e for evs in engine.all_cn_events.values() for e in evs)
        self.assertTrue(eligible(event))
        for support in event.supporting_muts:
            mutation = engine.mutations[support.var_str]
            self.assertTrue(eligible(mutation))
            np.testing.assert_allclose(mutation.pi_dist, support.pi_dist)
        comps = compare_events_original(engine, drivers=['KRAS'])
        self.assertTrue(any(p[0]+p[1] > 0 for p in comps.values()))
        with tempfile.TemporaryDirectory() as tmp:
            engine.patient.indiv_name = str(Path(tmp)/'S')
            PhylogicOutput().write_timing_tsv(engine)
            rows = read_tsv(Path(tmp)/'S.timing.tsv')
            drivers = [r for r in rows if r['Event Name'].startswith('KRAS_')]
            self.assertEqual(len(drivers),6)
            self.assertTrue(all(r['Timing_status']=='estimated' and r['pi_mean'] for r in drivers))

    def test_conflicting_local_mutation_estimates_are_order_independent(self):
        for reverse in [False, True]:
            engine = focal_engine()
            event = next(e for evs in engine.all_cn_events.values() for e in evs)
            import copy
            other = copy.copy(event)
            other.start, other.end = 1000001, 1999999
            events = [event, other]
            if reverse:
                events.reverse()
            engine.all_cn_events = {'focal':events}
            def gain_posterior(ev, **kwargs):
                ev.pi_dist = np.eye(101)[25 if ev.start == 1000000 else 75]
            with patch('SinglePatientTiming.TimingEngine.TimingCNEvent.get_pi_dist_for_gain', gain_posterior):
                engine.time_events()
            for mutation in engine.mutations.values():
                self.assertFalse(eligible(mutation))
                self.assertEqual(mutation.timing_reason, 'conflicting_local_gain_estimates')


class AdapterRegressions(unittest.TestCase):
    def parse(self, records, annotation_fields='Allele|SYMBOL|HGVSp|PICK|IMPACT|ALLELE_NUM'):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'input.vcf.gz'
            self.write_vcf(path, records, annotation_fields)
            out = io.StringIO()
            writer = csv.DictWriter(out, MAF_COLUMNS, delimiter='\t'); writer.writeheader()
            parse_vcf(path,'S',{'1':[(1,1000,1,1,100,100)]},{'1':[1]},writer)
            return list(csv.DictReader(io.StringIO(out.getvalue()), delimiter='\t'))

    @staticmethod
    def write_vcf(path, records, fields):
        with gzip.open(path,'wt') as handle:
            handle.write('##fileformat=VCFv4.2\n')
            handle.write(f'##INFO=<ID=PCGR_CSQ,Number=.,Type=String,Description="Format: {fields}">\n')
            handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n')
            handle.write('\n'.join(records)+'\n')

    def test_indel_annotation_uses_allele_number(self):
        rows = self.parse(['1\t10\t.\tATC\tA,AT\t.\tPASS\t'
            'PCGR_CSQ=-|TP53|ENSP:p.Arg10ValfsTer2|1|HIGH|1,-|KRAS|ENSP:p.Gly12del|1|HIGH|2\tGT:AD\t1/2:70,20,10'])
        self.assertEqual([r['Hugo_Symbol'] for r in rows], ['TP53','KRAS'])
        self.assertEqual([r['Protein_Change'] for r in rows], ['p.Arg10ValfsTer2','p.Gly12del'])

    def test_indel_annotation_fallback_handles_anchored_and_minimal_alleles(self):
        for ref,alt,annotated in [('AT','A','-'),('A','AT','T'),('ATC','AGC','G'),('A','T','T')]:
            rows = self.parse([f'1\t10\t.\t{ref}\t{alt}\t.\tPASS\t'
                f'PCGR_CSQ={annotated}|TP53|p.Arg10ValfsTer2|1|HIGH\tGT:AD\t0/1:70,30'],
                'Allele|SYMBOL|HGVSp|PICK|IMPACT')
            self.assertEqual(rows[0]['Hugo_Symbol'],'TP53')

    def test_ambiguous_indel_without_allele_number_stays_unknown(self):
        rows = self.parse(['1\t10\t.\tATC\tA,AT\t.\tPASS\t'
            'PCGR_CSQ=-|TP53|p.Arg10ValfsTer2|1|HIGH\tGT:AD\t1/2:70,20,10'],
            'Allele|SYMBOL|HGVSp|PICK|IMPACT')
        self.assertTrue(all(r['Hugo_Symbol']=='unknown' for r in rows))

    def test_multiallelic_depth_and_source_counts(self):
        rows = self.parse(['1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:AD\t1/2:70,20,10'])
        for row,alt,other in zip(rows,[20,10],[10,20]):
            self.assertEqual(int(row['t_ref_count'])+int(row['t_alt_count']),100)
            self.assertEqual(int(row['vcf_ref_count']),70)
            self.assertEqual(int(row['vcf_other_alt_count']),other)
            actual = calc_ccf(1,1,alt,int(row['t_ref_count']),.6)
            np.testing.assert_array_equal(actual,calc_ccf(1,1,alt,100-alt,.6))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'multi.maf'
            write_tsv(path,MAF_COLUMNS,rows)
            sample = TumorSample(str(path),'calc_ccf',sample_name='S',purity=.6,
                genome=Genome('hg38'),sex='FEMALE',min_coverage=8,
                artifact_blacklist=str(ROOT/'data/supplement_data/Blacklist_SNVs.txt'))
            for mutation in sample.mutations:
                self.assertEqual(mutation.alt_cnt+mutation.ref_cnt,100)
                np.testing.assert_allclose(mutation.ccf_1d,
                    calc_ccf(1,1,mutation.alt_cnt,100-mutation.alt_cnt,.6))

    def test_invalid_other_alt_count_does_not_corrupt_depth(self):
        self.assertEqual(self.parse(['1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:AD\t1/2:70,20,-10']), [])

    def test_preparation_converts_closed_segment_ends_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cn = [dict(chromosome=c,start=a,end=b,minorAlleleCopyNumber=n,majorAlleleCopyNumber=m)
                  for c,a,b,n,m in [('1',1,100,1,1),('1',101,101,1,2),
                                   ('X',2781479,2781479,0,1),('X',2781480,2781480,0,1)]]
            write_tsv(root/'cn.tsv',list(cn[0]),cn)
            write_tsv(root/'purity.tsv',['purity','ploidy','gender'],[dict(purity=.6,ploidy=2,gender='MALE')])
            records = [f'{c}\t{p}\t.\tA\tT\t.\tPASS\t.\tGT:AD\t0/1:70,30'
                       for c,p in [('1',100),('1',101),('X',2781479),('X',2781480)]]
            self.write_vcf(root/'input.vcf.gz',records,'Allele|SYMBOL')
            manifest = [dict(sample_id='S',vcf_path='input.vcf.gz',cn_path='cn.tsv',purity_path='purity.tsv')]
            write_tsv(root/'manifest.tsv',list(manifest[0]),manifest)
            prepare(root/'manifest.tsv',root/'prepared','/results')
            maf = read_tsv(root/'prepared/maf_with_cn/S.withCN.maf')
            self.assertEqual(len(maf),4)
            sample = cn_sample([],sex='MALE')
            tree = sample._resolve_CnEvents(str(root/'prepared/segments/S.alleliccapseg.tsv'),
                                           'alleliccapseg',.6,2,Genome('hg38').CHROMS)
            self.assertEqual(next(iter(tree['chr1'][100])).data[1]['cn_a2'],1)
            self.assertEqual(next(iter(tree['chr1'][101])).data[1]['cn_a2'],2)
            self.assertFalse(tree['chr1'][102])
            self.assertFalse(tree['chrX'][2781479])
            self.assertTrue(tree['chrX'][2781480])


class LeagueSeedRegressions(unittest.TestCase):
    def test_serial_and_parallel_repeat_with_fixed_seed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = []
            for sample in ['A','B','C','D']:
                for a,b in itertools.combinations(['KRAS_p.G12D','TP53_p.R175H','gain_chr1p','loss_chr9p'],2):
                    rows.append(dict(sample=sample,event1=a,event2=b,p_event1_win=.6,p_event2_win=.3,unknown=.1))
            write_tsv(root/'comparisons.tsv',list(rows[0]),rows)
            for jobs in [1,2]:
                results = []
                for repeat in [1,2]:
                    out = root/f'{jobs}_{repeat}';out.mkdir()
                    env = dict(os.environ,PYTHONHASHSEED='0',MPLBACKEND='Agg',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
                    subprocess.run([sys.executable,str(ROOT/'PhylogicNDT.py'),'LeagueModel','--cohort','SEED',
                        '--comparison_fn',str(root/'comparisons.tsv'),'--n_perms','4','--n_seasons','10',
                        '--n_jobs',str(jobs),'--perm_chunk_size','1','--seed','0'],
                        cwd=out,env=env,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
                    with (out/'SEEDfig_pkl.pkl').open('rb') as handle:
                        model = pickle.load(handle)
                    results.append(((out/'SEED.log_odds.tsv').read_text(),model.event_positions))
                self.assertEqual(results[0],results[1],f'n_jobs={jobs}')


if __name__ == '__main__':
    unittest.main()
