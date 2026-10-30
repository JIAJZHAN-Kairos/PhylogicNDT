import ast
import itertools
from pathlib import Path
from types import SimpleNamespace as NS
import sys
import os
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from SinglePatientTiming.TimingRepair import (nearest_cn_cluster, protein_change_kind, eligible,
                                              dedup_cn_events, is_driver, read_cluster_ccfs,
                                              aggregate_driver_events, assign_driver_record_ids)


def comparison_function(source):
    tree = ast.parse(source.read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'compare_events_original')
    ns = dict(np=np, itertools=itertools, is_driver=is_driver, eligible=eligible,
              _dedup_cn_events=dedup_cn_events, aggregate_driver_events=aggregate_driver_events)
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), 'exec'), ns)
    return ns['compare_events_original']


def event(name='gain_chr1p', status='estimated', center=20, clonal=True):
    p = np.zeros(101)
    p[center] = 1
    return NS(event_name=name, pi_dist=p, timing_status=status, is_clonal=clonal)


class Repairs(unittest.TestCase):
    def test_overlapping_driver_representations_not_independent_events(self):
        from SinglePatientTiming.TimingEngine import TimingMut
        args=dict(sample_list=[],gene='KRAS',chrN='chr12',pos=25227343,
                  alt_cnt=np.array([85]),ref_cnt=np.array([85]),local_cn_a1=np.array([1.]),
                  local_cn_a2=np.array([2.]),ccf_dist=np.eye(101)[100],
                  prot_change='p.Gln61Lys',cluster_assignment=1)
        a=TimingMut(alt='T',ref='G',**args)
        b=TimingMut(alt='TC',ref='GA',**args)
        for mut,center in [(a,40),(b,50)]:
            mut.pi_dist=np.eye(101)[center];mut.timing_status='estimated'
        assign_driver_record_ids([a,b],['KRAS'])
        self.assertNotEqual(a.event_name,b.event_name)
        grouped=aggregate_driver_events([a,b],['KRAS'])
        self.assertEqual(len(grouped),1)
        self.assertEqual(grouped[0].event_name,'KRAS_p.Gln61Lys')
        self.assertFalse(eligible(grouped[0]))
        self.assertEqual(grouped[0].source_record_count,2)

    def test_timing_network_has_bounded_event_nodes(self):
        from output.PhylogicOutput import PhylogicOutput
        import tempfile
        out=PhylogicOutput()
        comps={(str(i),str(j)):(.01,.01,.98) for i in range(90) for j in range(i+1,90)}
        graph=out._get_timing_graph(comps)
        self.assertEqual(len(graph),90)
        self.assertEqual(graph.number_of_edges(),0)
        large={(str(i),str(i+1)):(1.,0.,0.) for i in range(220)}
        old_cwd=os.getcwd()
        try:
            with tempfile.TemporaryDirectory() as directory:
                os.chdir(directory)
                self.assertEqual(out.draw_timing_graph('large',large),([],[],{}))
                self.assertTrue(Path('large.comp_graph.pdf').is_file())
        finally:
            os.chdir(old_cwd)
    def test_protein(self):
        for text in ['p.Gly12Asp', 'p.G12D', 'p.Arg282Trp', 'p.Trp53Ter', 'p.G12del', 'p.Ala23ValfsTer4']:
            self.assertEqual(protein_change_kind(text), 'nonsynonymous', text)
        for text in ['p.Gly12Gly', 'p.G12G', 'p.Gly12=']:
            self.assertEqual(protein_change_kind(text), 'synonymous', text)
        self.assertEqual(protein_change_kind('garbage'), 'unparsed')

    def test_cn_mode_and_ties(self):
        grids = {c: np.eye(101)[[mode]] for c, mode in [(1,100),(2,67),(3,42),(4,24),(5,7)]}
        self.assertEqual(nearest_cn_cluster([.99], grids), 1)
        self.assertEqual(nearest_cn_cluster([.758126], grids), 2)
        self.assertIsNone(nearest_cn_cluster([.5], {1:np.eye(101)[[100]],2:np.eye(101)[[0]]}))
        self.assertIsNone(nearest_cn_cluster([0], {1:np.eye(101)[[100]]}))

    def test_allele_conflict_and_permutation(self):
        a,b = event(),event(status='subclonal',center=100,clonal=False)
        for evs in [[a,b],[b,a]]:
            self.assertFalse(eligible(dedup_cn_events({'x':evs})[0]))
        a,b = event(),event(status='prior_only',center=70)
        self.assertTrue(eligible(dedup_cn_events({'x':[b,a]})[0]))
        self.assertTrue(eligible(dedup_cn_events({'x':[a,a]})[0]))
        self.assertFalse(eligible(dedup_cn_events({'x':[a,event(center=70)]})[0]))

    def test_uninformative_unknown_and_presence(self):
        compare = comparison_function(ROOT/'SinglePatientTiming/SinglePatientTiming.py')
        a,b = event(),event(name='loss_chr2p', status='prior_only')
        engine = NS(WGD=None, all_cn_events={'a':[a],'b':[b]}, mutations={})
        self.assertEqual(compare(engine)[(a.event_name,b.event_name)],(0,0,1))
        a.pi_dist = np.ones(101)/101
        self.assertFalse(eligible(a))

    def test_original_formula_unchanged_for_eligible(self):
        fixed = comparison_function(ROOT/'SinglePatientTiming/SinglePatientTiming.py')
        original = ROOT/'tests/reference'
        old = comparison_function(original/'comparison_original.py')
        rng = np.random.RandomState(42)
        for _ in range(30):
            a,b = event(),event(name='gain_chr2q')
            a.pi_dist,b.pi_dist = rng.dirichlet(np.ones(101)),rng.dirichlet(np.ones(101))
            engine = NS(WGD=None, all_cn_events={'a':[a],'b':[b]}, mutations={})
            self.assertEqual(fixed(engine),old(engine))

    def test_native_allele_ids_and_homdel(self):
        from SinglePatientTiming.TimingEngine import TimingCNState
        state = TimingCNState([], 'chr3', 'q', (2.,2.), [], supporting_muts=[])
        self.assertEqual({e.allele for e in state.cn_events}, {'a1','a2'})
        self.assertEqual(len({e.detail_event_name for e in state.cn_events}),2)
        homdel = TimingCNState([], 'chr9', 'p', (0.,0.), [], supporting_muts=[], cn_type_prefix='Focal')
        self.assertEqual(len(homdel.cn_events),1)
        self.assertEqual(homdel.cn_events[0].allele,'both')
        homdel.call_events(cn_type_prefix='Focal',wgd=True)
        self.assertEqual([e.Type for e in homdel.cn_events],['Focal_homdel'])

    def test_focal_reconstruction_preserves_presence_without_inventing_timing(self):
        from SinglePatientTiming.TimingEngine import TimingEngine
        from intervaltree import IntervalTree
        import pandas as pd
        def engine(category, pair, wgd=False, caller_cn=None):
            cn=NS(type='CNV',cn_category=category,chrN='chr9',start=100,end=200,arm='p',a1=True,ccf_hat=1.,
                  local_cn_a1=pair[0] if caller_cn is None else caller_cn,local_cn_a2=pair[1])
            profile=IntervalTree()
            profile[100:200]=('sample',{'cn_a1':pair[0],'cn_a2':pair[1]})
            sample=NS(sample=NS(low_coverage_mutations={'cn':cn}),CnProfile={'chr9':profile},
                      mutation_intervaltree={'chr9':IntervalTree()},purity=.8,sample_name='sample')
            e=TimingEngine.__new__(TimingEngine)
            e.sample_list=[sample];e.genome=NS(cytoband_table=pd.DataFrame())
            e.patient=NS(saved_cluster_ccfs={1:np.eye(101)[[100]]})
            e.all_cn_events={};e.truncal_cn_events={};e.timeable_muts={};e.WGD=NS() if wgd else None;e.borrow_arm_prior=True
            e.get_focal_cn_events()
            return next(iter(e.all_cn_events.values()))[0] if e.all_cn_events else None
        mismatch=engine('Focal_loss',(0.,0.))
        self.assertEqual(mismatch.Type,'Focal_loss')
        self.assertIsNone(mismatch.cluster_id)
        self.assertIsNone(mismatch.pi_dist)
        self.assertEqual(mismatch.reconstruction_reason,'no_matching_local_cn_event')
        gain=engine('Focal_gain',(2.,2.))
        self.assertEqual(gain.allele,'a1')
        self.assertEqual(gain.borrowed_prior_key[-1],'a1')
        self.assertIsNone(engine('Focal_gain',(2.,2.),wgd=True))
        qualifying=engine('Focal_gain',(2.,2.),wgd=True,caller_cn=3.)
        self.assertEqual(qualifying.reconstruction_reason,'no_matching_local_cn_event')

    def test_missing_timing_segment_preserves_input_cn_without_event_support(self):
        from SinglePatientTiming.TimingEngine import TimingSample
        from intervaltree import IntervalTree
        mutation = NS(chrN='chr1',pos=20,type='SNP',var_str='chr1:20:A:G',gene='KRAS',alt='G',ref='A',
                      alt_cnt=20,ref_cnt=20,ccf_1d=np.eye(101)[100],prot_change='p.Gly12Asp',
                      cluster_assignment=1,local_cn_a1=1.5143,local_cn_a2=2.0188)
        sample = TimingSample.__new__(TimingSample)
        sample.sample=NS(mutations=[mutation],low_coverage_mutations={})
        sample.CnProfile={'chr1':IntervalTree()}
        sample._chromosomes=['chr1']
        sample.purity=.99
        sample.fill_mutation_intervaltree()
        observed=sample.mut_lookup_table[mutation.var_str]
        self.assertEqual(observed.local_cn_a1[0],1.5143)
        self.assertEqual(observed.local_cn_a2[0],2.0188)
        self.assertEqual(len(sample.mutation_intervaltree['chr1']),0)

    def test_wgd_cumulative_roundoff_and_zero_mass(self):
        from SinglePatientTiming.TimingEngine import TimingCNState
        state=TimingCNState([], 'chr1','p',(0.,4.),[],supporting_muts=[])
        weights=np.arange(1.,102.)
        wgd=NS(pi_dist=weights/weights.sum()*(1+1e-15))
        loss,gain=state.cn_events
        loss.get_pi_dist_for_loss(wgd)
        gain.get_pi_dist_for_higher_gain(wgd)
        for e in (loss,gain):
            self.assertTrue((e.pi_dist>=0).all())
            self.assertAlmostEqual(e.pi_dist.sum(),1)
        wgd.pi_dist=np.eye(101)[0]
        loss.get_pi_dist_for_loss(wgd)
        self.assertIsNone(loss.pi_dist)


if __name__ == '__main__':
    unittest.main()
