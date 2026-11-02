"""Generate fabricated SASH-like VCF/PURPLE files; no patient data."""
import argparse
import csv
import gzip
import random
from pathlib import Path


def write(path,fields,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def make(root):
    root.mkdir(parents=True,exist_ok=True)
    rows=[]
    rng=random.Random(1234)
    for sample,sex in [('DEMO_MALE','MALE'),('DEMO_FEMALE','FEMALE')]:
        segments=[('1',1,248956422,1,2),('8',1,45200000,1,1),('8',45200001,145138636,1,3),
                  ('9',1,138394717,0,1),('12',1,133275309,1,1),('17',1,83257441,0,1),
                  ('X',1,156040895,0 if sex=='MALE' else 1,1),('Y',1,57227415,0,0)]
        cn=[dict(chromosome=c,start=a,end=b,minorAlleleCopyNumber=n,majorAlleleCopyNumber=m,
                 bafCount=100,depthWindowCount=100) for c,a,b,n,m in segments]
        write(root/(sample+'.purple.cnv.somatic.tsv'),list(cn[0]),cn)
        purity=dict(purity=.6,ploidy=2.2,gender=sex,wholeGenomeDuplication='FALSE',germlineAberrations='NONE')
        write(root/(sample+'.purple.purity.tsv'),list(purity),[purity])
        vcf=['##fileformat=VCFv4.2','##reference=GRCh38',
             '##INFO=<ID=PCGR_CSQ,Number=.,Type=String,Description="Synthetic annotation. Format: Allele|Consequence|IMPACT|SYMBOL|HGVSp|PICK">',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
             '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\t'+sample]
        variants=[]
        for chrom,pos,gene,protein in [('12',25245350,'KRAS','p.Gly12Asp'),('17',7674220,'TP53','p.Arg175His'),
                                       ('X',3000000,'DEMO','p.Ala1Val'),('X',20000,'DEMO','p.Ala2Val')]:
            variants.append((chrom,pos,gene,protein,1.))
        for chrom,start in [('1',10000000),('8',90000000),('9',10000000),('12',30000000),('17',10000000)]:
            for i in range(60):
                variants.append((chrom,start+1000*i,'DEMO','p.Ala3Val',.5 if i%5==0 else 1.))
        for chrom,pos,gene,protein,ccf in variants:
            minor,major=next((a,b) for c,s,e,a,b in segments if c==chrom and s<=pos<=e)
            ncn=1 if sex=='MALE' and chrom=='X' and pos>2781479 else 2
            mult=2 if chrom in {'1','8'} and pos%2000==0 else 1
            probability=.6*mult*ccf/(.6*(minor+major)+.4*ncn)
            alt=max(1,min(99,sum(rng.random()<probability for _ in range(100))))
            info=f'PCGR_CSQ=T|missense_variant|MODERATE|{gene}|{protein}|1'
            vcf.append(f'{chrom}\t{pos}\t.\tA\tT\t.\tPASS\t{info}\tGT:AD\t0/0:100,0\t0/1:{100-alt},{alt}')
        with gzip.open(root/(sample+'.sash.vcf.gz'),'wt') as f:f.write('\n'.join(vcf)+'\n')
        rows.append(dict(sample_id=sample,vcf_path=sample+'.sash.vcf.gz',
                         cn_path=sample+'.purple.cnv.somatic.tsv',purity_path=sample+'.purple.purity.tsv',
                         vcf_sample=sample))
    write(root/'manifest.tsv',list(rows[0]),rows)
    (root/'driver_genes.txt').write_text('KRAS\nTP53\nCDKN2A\nSMAD4\n')
    (root/'focal_peaks.bed').write_text('chr8\t127188796\t127328999\tAMP\nchr9\t21865961\t21995127\tDEL\n')
    for g,sample in [('demo_male','DEMO_MALE'),('demo_female','DEMO_FEMALE')]:
        (root/(g+'.samples.txt')).write_text(sample+'\n')
    (root/'README.txt').write_text('All records and read counts are fabricated. These are format/smoke-test inputs, not biological evidence.\n')
    print('Synthetic files:',root)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,required=True)
    make(parser.parse_args().outdir)
