import argparse
import base64
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow'))
import run_nf_sample
import run_nf_cohort
from prepare_sash_purple import read_segments, local_cn


class NextflowGuardsTest(unittest.TestCase):
    def test_negative_cn_is_zero_and_original_values_are_audited(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'cn.tsv'
            path.write_text('chromosome\tstart\tend\tminorAlleleCopyNumber\tmajorAlleleCopyNumber\n'
                            'chr1\t1\t100\t0\t-0.2752\n'
                            'chr1\t101\t200\t-0.01\t1.4\n')
            audit=[]
            segments, starts=read_segments(path,audit)
            self.assertEqual(local_cn(100,segments['1'],starts['1']),(0.,0.))
            self.assertEqual(local_cn(101,segments['1'],starts['1']),(0.,1.4))
            self.assertEqual(audit[0]['original_major'],-0.2752)
            self.assertEqual(audit[0]['clipped_major'],0.)

    def test_incomplete_cohort_cannot_run_league(self):
        import os
        old=Path.cwd()
        with tempfile.TemporaryDirectory() as d:
            try:
                os.chdir(d)
                # Nextflow creates this shell script beside staged sample directories.
                Path('.command.run').write_text('#!/bin/bash\n')
                with self.assertRaisesRegex(ValueError,'Incomplete'):
                    run_nf_cohort.run(argparse.Namespace(expected=380))
                self.assertFalse(Path('cohort').exists())
            finally:
                os.chdir(old)

    def test_bad_input_checksum_stops_before_preparation(self):
        import os
        old=Path.cwd()
        with tempfile.TemporaryDirectory() as d:
            try:
                os.chdir(d)
                Path('input').write_text('wrong bytes')
                meta=base64.b64encode(json.dumps(dict(sample_id='S',vcf_sha256='0'*64)).encode()).decode()
                args=argparse.Namespace(meta=meta,vcf='input',cn='input',purity='input')
                with self.assertRaisesRegex(ValueError,'checksum mismatch'):
                    run_nf_sample.run(args)
                self.assertFalse(Path('S.run/input').exists())
            finally:
                os.chdir(old)
