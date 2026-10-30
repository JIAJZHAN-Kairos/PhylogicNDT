import csv
import gzip
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'workflow'))
from prepare_sash_purple import prepare, read_tsv, write_tsv, alleliccapseg_values
from prepare_comparisons import aggregate, FIELDS, LEGACY
spec=importlib.util.spec_from_file_location('make_demo',ROOT/'examples/make_demo.py')
demo=importlib.util.module_from_spec(spec);spec.loader.exec_module(demo)


class Inputs(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.raw=self.root/'raw'
        demo.make(self.raw)
        self.manifest=self.raw/'manifest.tsv';self.out=self.root/'prepared'

    def run_prepare(self):
        prepare(self.manifest,self.out,'/results','/input')

    def change_manifest(self,**changes):
        rows=read_tsv(self.manifest);rows[0].update(changes)
        write_tsv(self.manifest,list(rows[0]),rows)

    def test_mapping_counts_and_independent_sifs(self):
        self.run_prepare()
        for sample,sex in [('DEMO_MALE','MALE'),('DEMO_FEMALE','FEMALE')]:
            rows=read_tsv(self.out/'cluster_sif'/(sample+'.sif'))
            self.assertEqual(len(rows),1)
            self.assertEqual(list(rows[0]),['sample_id','maf_fn','seg_fn','purity','timepoint','ploidy','sex','par_mode'])
            self.assertEqual(rows[0]['sex'],sex)
            self.assertEqual(rows[0]['maf_fn'],f'/input/maf_with_cn/{sample}.withCN.maf')
            timing=read_tsv(self.out/'timing_sif'/(sample+'.sif'))[0]
            self.assertEqual(timing['cluster_ccf_fn'],f'/results/{sample}/{sample}.cluster_ccfs.txt')
            maf=read_tsv(self.out/'maf_with_cn'/(sample+'.withCN.maf'))
            kras=next(r for r in maf if r['Hugo_Symbol']=='KRAS')
            self.assertEqual(kras['Protein_Change'],'p.Gly12Asp')
            with gzip.open(self.raw/(sample+'.sash.vcf.gz'),'rt') as f:
                raw=next(x for x in f if not x.startswith('#') and x.split('\t')[0]=='12')
            ad=raw.rstrip().split('\t')[-1].split(':')[1].split(',')
            self.assertEqual([kras['t_ref_count'],kras['t_alt_count']],ad)
            self.assertTrue(any(r['Chromosome']=='X' and r['Start_position']=='20000' for r in maf))
        provenance=json.loads((self.out/'input_provenance.json').read_text())
        self.assertEqual(len(provenance['samples']),2)
        self.assertEqual(provenance['samples'][0]['par_mode'],'unresolved')

    def test_exact_tumour_column_required(self):
        self.change_manifest(vcf_sample='DOES_NOT_EXIST')
        with self.assertRaisesRegex(ValueError,'explicitly named tumour'):
            self.run_prepare()

    def test_runtime_alias_is_not_vcf_sample_guess(self):
        self.change_manifest(sample_id='RENAMED')
        self.run_prepare()
        maf=read_tsv(self.out/'maf_with_cn/RENAMED.withCN.maf')
        self.assertTrue(all(r['Tumor_Sample_Barcode']=='RENAMED' for r in maf))

    def test_sex_mismatch_rejected(self):
        self.change_manifest(sex='FEMALE')
        with self.assertRaisesRegex(ValueError,'concordant'):
            self.run_prepare()

    def test_duplicate_samples_rejected(self):
        self.change_manifest(sample_id='DEMO_FEMALE')
        with self.assertRaisesRegex(ValueError,'unique'):
            self.run_prepare()

    def test_abnormal_germline_cn_requires_review(self):
        p=self.raw/'DEMO_MALE.purple.purity.tsv';rows=read_tsv(p)
        rows[0]['germlineAberrations']='KLINEFELTER'
        write_tsv(p,list(rows[0]),rows)
        with self.assertRaisesRegex(ValueError,'individual review'):
            self.run_prepare()

    def test_overlapping_cn_rejected(self):
        p=self.raw/'DEMO_MALE.purple.cnv.somatic.tsv';rows=read_tsv(p)
        rows.append(dict(rows[0]));write_tsv(p,list(rows[0]),rows)
        with self.assertRaisesRegex(ValueError,'Overlapping'):
            self.run_prepare()

    def test_filtering_multiallelic_and_unknown_annotation(self):
        p=self.raw/'DEMO_MALE.sash.vcf.gz'
        records=[
            '12\t50000001\t.\tA\tC,G\t.\tPASS\t.\tGT:AD\t0/0:100,0,0\t0/1:70,20,10',
            '12\t50000002\t.\tA\tT\t.\tLowQual\t.\tGT:AD\t0/0:100,0\t0/1:70,30',
            '12\t50000003\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1',
            'Y\t5000000\t.\tA\tT\t.\tPASS\t.\tGT:AD\t0/0:100,0\t0/1:70,30',
            '22\t5000000\t.\tA\tT\t.\tPASS\t.\tGT:AD\t0/0:100,0\t0/1:70,30',
            '12\t50000004\t.\tA\t<DEL>\t.\tPASS\t.\tGT:AD\t0/0:100,0\t0/1:70,30']
        with gzip.open(p,'at') as f:f.write('\n'.join(records)+'\n')
        self.run_prepare()
        maf=read_tsv(self.out/'maf_with_cn/DEMO_MALE.withCN.maf')
        multi=[r for r in maf if r['Start_position']=='50000001']
        self.assertEqual([(r['Tumor_Seq_Allele2'],r['t_alt_count']) for r in multi],[('C','20'),('G','10')])
        self.assertTrue(all(r['Hugo_Symbol']=='unknown' for r in multi))
        qc=read_tsv(self.out/'preparation_qc.tsv')[0]
        for key in ['homozygous_deletion','no_cn_segment','missing_ad','unsupported_alleles']:
            self.assertEqual(qc[key],'1')

    def test_mu_bridge_roundtrip_is_not_sex_recalibrated_twice(self):
        for cn in [0.,1.,2.,54.]:
            for purity in [.3,.6,.99]:
                mu,sigma=alleliccapseg_values(cn,purity,2.6,100,400)
                self.assertAlmostEqual(1+(mu*2.6/2-1)/purity,cn)
                self.assertGreater(sigma,0)

    def test_existing_preparation_is_preserved(self):
        self.run_prepare()
        with self.assertRaisesRegex(ValueError,'new or empty'):
            self.run_prepare()


class Comparisons(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.samples=self.root/'samples.txt';self.samples.write_text('SYNTHETIC_A\n')
        self.result=self.root/'SYNTHETIC_A';self.result.mkdir()
        self.file=self.result/'SYNTHETIC_A.comp.tsv'
        self.output=self.root/'group.comparisons.tsv'

    def write(self,prob=(.8,.1,.1),sample='SYNTHETIC_A',legacy=False):
        fields=LEGACY if legacy else FIELDS
        values=[sample,'KRAS','gain_chr1p',*prob]
        write_tsv(self.file,fields,[dict(zip(fields,values))])

    def test_legacy_header_and_unknown_presence_retained(self):
        self.write(prob=(0,0,1),legacy=True)
        aggregate(self.samples,self.root,self.output)
        rows=read_tsv(self.output)
        self.assertEqual(list(rows[0]),FIELDS)
        report=json.loads(self.output.with_suffix('.audit.json').read_text())
        self.assertEqual(report['inputs'][0]['informative_rows'],0)

    def test_wrong_sample_rejected(self):
        self.write(sample='WRONG')
        with self.assertRaisesRegex(ValueError,'mismatch'):
            aggregate(self.samples,self.root,self.output)

    def test_invalid_probabilities_rejected(self):
        self.write(prob=(.9,.9,0))
        with self.assertRaisesRegex(ValueError,'probabilities'):
            aggregate(self.samples,self.root,self.output)

    def test_header_only_rejected(self):
        self.file.write_text('\t'.join(FIELDS)+'\n')
        with self.assertRaisesRegex(ValueError,'Header-only'):
            aggregate(self.samples,self.root,self.output)


if __name__=='__main__':unittest.main()
