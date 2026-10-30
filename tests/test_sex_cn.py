import importlib.util
from pathlib import Path
from types import SimpleNamespace as NS
import sys
import tempfile
import unittest

import numpy as np
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from utils.sex_cn import baseline, normal_cn, region_baseline
from utils.calc_ccf import calc_ccf
from data.Sample import TumorSample
from data.Enums import Genome
from SinglePatientTiming.TimingEngine import TimingCNState, TimingMut, TimingEngine, TimingWGD


def sample(sex):
    return NS(sex=sex, par_mode='unresolved', purity=.6, sample_name='test')


class SexCN(unittest.TestCase):
    def test_neutral_post_wgd_states_remain_clock_candidates(self):
        for sex, chrom, copies in [('FEMALE', 'chr1', (2, 2)), ('MALE', 'chrX', (0, 2))]:
            smp = sample(sex)
            state = TimingCNState([smp], chrom, 'p', copies, [.6], supporting_muts=[])
            state.call_events(wgd=True)
            self.assertEqual(state.cn_events, [])
            smp.WGD, smp.cn_states = NS(), {chrom + 'p': state}
            engine = TimingEngine.__new__(TimingEngine)
            engine.sample_list, engine.arm_regions = [smp], [(chrom, 'p')]
            engine.get_concordant_cn_states()
            self.assertIn(chrom + 'p', engine.concordant_cn_states)

    def test_pooled_wgd_clock_is_counted_once_per_region(self):
        distribution = np.arange(1, 102, dtype=float)
        distribution /= distribution.sum()
        events = [NS(Type='Arm_gain', pi_dist=distribution.copy(), supporting_muts=[1, 2, 3]) for _ in range(2)]
        wgd = TimingWGD([NS(cn_events=events, supporting_muts=[])])
        wgd.get_pi_dist()
        np.testing.assert_allclose(wgd.pi_dist, distribution)
        self.assertEqual(wgd.n_clock_regions, 1)
        self.assertEqual(sum(e.wgd_clock_used for e in events), 1)

    def test_uninformative_wgd_cannot_create_cn_timing(self):
        wgd = TimingWGD([])
        wgd.get_pi_dist()
        self.assertEqual(wgd.timing_status, 'not_estimable')
        state = TimingCNState([sample('FEMALE')], 'chr1', 'p', (0, 2), [.6], supporting_muts=[])
        event = next(e for e in state.cn_events if e.Type.endswith('loss'))
        event.get_pi_dist_for_loss(wgd)
        self.assertIsNone(event.pi_dist)

    def test_subclonal_mutations_cannot_create_a_clonal_gain_clock(self):
        mutations = [NS(is_clonal=False, local_cn_a1=np.array([0.]), local_cn_a2=np.array([2.])) for _ in range(4)]
        state = TimingCNState([sample('MALE')], 'chrX', 'p', (0, 2), [.6], supporting_muts=mutations)
        event = state.cn_events[0]
        self.assertEqual(event.timeable_supporting_muts, [])
        event.get_pi_dist_for_gain()
        self.assertIsNone(event.pi_dist)

    def test_negative_sex_cn_in_maf_is_bounded_and_audited(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'input.maf'
            path.write_text('Hugo_Symbol\tChromosome\tStart_position\tEnd_position\tReference_Allele\tTumor_Seq_Allele2\tVariant_Classification\tVariant_Type\tProtein_Change\tt_ref_count\tt_alt_count\tlocal_cn_a1\tlocal_cn_a2\n'
                            'TEST\tchrX\t3000000\t3000000\tA\tT\tMissense_Mutation\tSNP\tp.Ala1Val\t85\t15\t-0.2\t2\n')
            obj = TumorSample(str(path), 'calc_ccf', sample_name='test', purity=.6, ploidy=2.6,
                              genome=Genome('hg38'), sex='MALE', min_coverage=8,
                              artifact_blacklist=str(ROOT/'data/supplement_data/Blacklist_SNVs.txt'))
            mut = obj.mutations[0]
            self.assertEqual(mut.local_cn_a1, 0)
            np.testing.assert_allclose(mut.ccf_1d, calc_ccf(0, 2, 15, 85, .6, normal_cn=1))
            self.assertTrue(any(reason=='negative_tumor_CN_estimate_clipped_to_zero' for _, reason in obj.sex_cn_unresolved))

    def test_normal_baselines_and_par_boundaries(self):
        male, female = sample('MALE'), sample('FEMALE')
        self.assertEqual(baseline(male, 'chrX'), (0, 1))
        self.assertEqual(baseline(female, 'chrX'), (1, 1))
        self.assertEqual(baseline(female, 'chrY'), (0, 0))
        for chrom, start, end in [('X', 10001, 2781479), ('X', 155701383, 156030895),
                                  ('Y', 10001, 2781479), ('Y', 56887903, 57217415)]:
            self.assertEqual(normal_cn(male, chrom, start - 1), 1)
            self.assertIsNone(normal_cn(male, chrom, start))
            self.assertIsNone(normal_cn(male, chrom, end))
            self.assertEqual(normal_cn(male, chrom, end + 1), 1)
        male.par_mode = 'x_diploid_y_masked'
        self.assertEqual(normal_cn(male, 'X', 10001), 2)
        self.assertIsNone(normal_cn(male, 'Y', 10001))
        self.assertIsNone(region_baseline(male, 'X', 9999, 10002))
        with self.assertRaises(ValueError):
            baseline(NS(), 'X', 3000000)

    def test_neutral_absent_alleles_do_not_generate_losses(self):
        for sex, chrom, copies, wgd in [('MALE', 'X', (0, 1), False), ('MALE', 'Y', (0, 1), False),
                                         ('MALE', 'X', (0, 2), True), ('FEMALE', 'X', (1, 1), False),
                                         ('FEMALE', 'X', (2, 2), True), ('FEMALE', 'Y', (0, 0), False)]:
            state = TimingCNState([sample(sex)], chrom, 'q', copies, [.6], supporting_muts=[])
            state.call_events(wgd=wgd)
            self.assertEqual(state.cn_events, [], (sex, chrom, copies, wgd))

    def test_true_haploid_changes_remain(self):
        state = TimingCNState([sample('MALE')], 'Y', 'q', (0, 0), [.6], supporting_muts=[])
        self.assertEqual([(e.Type, e.allele) for e in state.cn_events], [('Arm_loss', 'a2')])
        event = state.cn_events[0]
        event.get_pi_dist_for_loss(NS(pi_dist=np.ones(101) / 101))
        self.assertIsNone(event.pi_dist)
        state.call_events(cn_type_prefix='Focal')
        self.assertEqual(state.cn_events[0].Type, 'Focal_loss')
        state = TimingCNState([sample('MALE')], 'X', 'p', (0, 1), [.6], supporting_muts=[])
        state.call_events(wgd=True)
        self.assertEqual([(e.Type, e.allele) for e in state.cn_events], [('Arm_loss', 'a2')])
        state = TimingCNState([sample('MALE')], 'X', 'p', (0, 2), [.6], supporting_muts=[])
        self.assertEqual([(e.Type, e.allele) for e in state.cn_events], [('Arm_gain', 'a2')])

    def test_ccf_autosome_exact_and_haploid_changes(self):
        old = ROOT/'tests/reference/calc_ccf_pre_sexaware.py'
        spec = importlib.util.spec_from_file_location('old_ccf', old)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        for a1, a2 in [(0, 1), (1, 1), (.3, 1.8), (2, 3), (3, 54)]:
            for purity in [.3, .6, .99]:
                args = (a1, a2, 15, 85, purity)
                np.testing.assert_array_equal(calc_ccf(*args), mod.calc_ccf(*args))
        diploid = calc_ccf(0, 1, 30, 70, .6)
        haploid = calc_ccf(0, 1, 30, 70, .6, normal_cn=1)
        self.assertLess(np.argmax(haploid), np.argmax(diploid))
        self.assertAlmostEqual(haploid.sum(), 1)

    def test_timing_multiplicity_uses_matching_normal_cn(self):
        for sex, ncn in [('MALE', 1), ('FEMALE', 2)]:
            mut = TimingMut([sample(sex)], 'TEST', 'X', 3000000, 'T', 'A', np.array([30]),
                            np.array([70]), np.array([0.]), np.array([2.]), np.eye(101)[100])
            actual = mut.get_multiplicity_likelihoods()[1][0]
            self.assertAlmostEqual(actual, binom.pmf(30, 100, .6 / (2 * .6 + ncn * .4)))

    def test_legacy_mu_decode_and_par_split(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'cn.tsv'
            for major in [0., 1., 2.]:
                # absolute 0/major encoded with legacy autosomal bridge.
                path.write_text('Chromosome\tStart.bp\tEnd.bp\tmu.minor\tsigma.minor\tmu.major\tsigma.major\n'
                                'chrX\t10000\t3000000\t{}\t.01\t{}\t.01\n'.format(
                                    2 * (.4 + .6 * (-.17 if major == 0 else 0)) / 2.6, 2 * (.4 + .6 * major) / 2.6))
                obj = TumorSample.__new__(TumorSample)
                obj.sample_name, obj.sex, obj.par_mode = 'test', 'MALE', 'unresolved'
                obj.sex_cn_unresolved = []
                tree = obj._resolve_CnEvents(str(path), 'alleliccapseg', purity=.6, ploidy=2.6,
                                            chroms=Genome('hg38').CHROMS)['chrX']
                self.assertFalse(tree[20000])
                row = next(iter(tree[2800000])).data[1]
                self.assertEqual(row['cn_a1'], 0)
                self.assertEqual(row['ccf_hat_a1'], 0)
                self.assertAlmostEqual(row['cn_a2'], major)
                self.assertTrue(obj.sex_cn_unresolved)


if __name__ == '__main__':
    unittest.main()
