#!/usr/bin/env python3
"""Prepare explicit per-specimen PhylogicNDT inputs from SASH VCF + PURPLE files.

Homozygous-deletion (total CN=0) sites are omitted because they
cannot yield a valid mutation CCF posterior.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import re
from collections import Counter
import hashlib
import json
from pathlib import PurePosixPath
from pathlib import Path


MAF_COLUMNS = [
    "Hugo_Symbol", "Chromosome", "Start_position", "End_position",
    "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
    "Variant_Type", "Protein_Change", "Tumor_Sample_Barcode", "t_ref_count", "t_alt_count",
    "local_cn_a1", "local_cn_a2", "vcf_ref_count", "vcf_other_alt_count",
]
STANDARD_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}


def canonical_chromosome(value: str) -> str:
    value = value.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def read_segments(path: Path, clipped_segments=None) -> tuple[dict[str, list[tuple[int, int, float, float, int, int]]], dict[str, list[int]]]:
    by_chrom: dict[str, list[tuple[int, int, float, float, int, int]]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chromosome", "start", "end", "minorAlleleCopyNumber", "majorAlleleCopyNumber"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError(f"PURPLE CN file misses required columns: {path}")
        for row in reader:
            chrom = canonical_chromosome(row["chromosome"])
            if chrom not in STANDARD_CHROMS:
                continue
            try:
                start, end = int(row["start"]), int(row["end"])
                minor, major = float(row["minorAlleleCopyNumber"]), float(row["majorAlleleCopyNumber"])
                baf_count = int(float(row.get("bafCount", 0)))
                depth_windows = int(float(row.get("depthWindowCount", 0)))
            except ValueError:
                continue
            if start < 1 or end < start or not all(map(math.isfinite, (minor, major))):
                continue
            if minor < 0 or major < 0:
                if clipped_segments is not None:
                    clipped_segments.append(dict(chromosome=chrom, start=start, end=end,
                        original_minor=minor, original_major=major,
                        clipped_minor=max(0., minor), clipped_major=max(0., major)))
                minor, major = max(0., minor), max(0., major)
            by_chrom.setdefault(chrom, []).append((start, end, min(minor, major), max(minor, major), baf_count, depth_windows))
    starts: dict[str, list[int]] = {}
    for chrom, rows in by_chrom.items():
        rows.sort()
        starts[chrom] = [row[0] for row in rows]
    return by_chrom, starts


def local_cn(position: int, rows: list[tuple[int, int, float, float, int, int]], starts: list[int]) -> tuple[float, float] | None:
    index = bisect.bisect_right(starts, position) - 1
    if index < 0:
        return None
    start, end, minor, major, _, _ = rows[index]
    return (minor, major) if start <= position <= end else None


def variant_type(ref: str, alt: str) -> tuple[str, str]:
    if len(ref) == len(alt) == 1:
        return "SNV", "SNP"
    if len(ref) == len(alt):
        return "MNV", "MNP"
    if len(ref) > len(alt):
        return "Frame_Shift_Del", "DEL"
    return "Frame_Shift_Ins", "INS"


def alleliccapseg_values(allele_cn: float, purity: float, ploidy: float, baf_count: int, depth_windows: int) -> tuple[float, float]:
    """Map PURPLE absolute allele CN to Phylogic's mu/sigma representation.

    PURPLE does not report an allele-CN uncertainty in this TSV.  sigma is a
    conservative support-based proxy and is labelled as such in preparation_qc.
    """
    mu = 2.0 * ((1.0 - purity) + purity * allele_cn) / ploidy
    baf_se = 0.5 / math.sqrt(max(1, baf_count))
    depth_se = 0.5 / math.sqrt(max(1, depth_windows))
    cn_sigma = max(0.05, baf_se + depth_se)
    sigma_mu = max(0.01, 2.0 * purity * cn_sigma / ploidy)
    return mu, sigma_mu


def gene_annotation(info: str, alt: str, csq_fields: list[str] | None,
                    ref: str = '', alt_index: int = 1, all_alts=None) -> tuple[str, str]:
    """Return PCGR/VEP SYMBOL and a concise non-synonymous protein change."""
    if not csq_fields or "SYMBOL" not in csq_fields:
        return "unknown", ""
    csq = next((value for value in info.split(";") if value.startswith("PCGR_CSQ=")), "")
    if not csq:
        return "unknown", ""
    symbol_index = csq_fields.index("SYMBOL")
    protein_index = csq_fields.index("HGVSp") if "HGVSp" in csq_fields else None
    pick_index = csq_fields.index("PICK") if "PICK" in csq_fields else None
    impact_index = csq_fields.index("IMPACT") if "IMPACT" in csq_fields else None
    number_index = csq_fields.index("ALLELE_NUM") if "ALLELE_NUM" in csq_fields else None
    all_alts = all_alts or [alt]
    allele_indices = {}
    trim_anchor = (ref and any(len(a) != len(ref) for a in all_alts)
                   and all(a[0] == ref[0] for a in all_alts))
    for index, alternate in enumerate(all_alts, start=1):
        aliases = {alternate}
        if trim_anchor:
            aliases.add(alternate[1:] or '-')
        # Also accept VEP --minimal representation. Ambiguous aliases (e.g.
        # two deletions both represented by '-') require ALLELE_NUM.
        r, a = ref, alternate
        while r and a and r[0] == a[0]:
            r, a = r[1:], a[1:]
        while r and a and r[-1] == a[-1]:
            r, a = r[:-1], a[:-1]
        aliases.add(a or '-')
        for allele in aliases:
            allele_indices.setdefault(allele, set()).add(index)
    candidates = []
    for annotation in csq.removeprefix("PCGR_CSQ=").split(","):
        values = annotation.split("|")
        if len(values) <= symbol_index or not values[symbol_index]:
            continue
        number = values[number_index] if number_index is not None and len(values) > number_index else ''
        if number not in ('', '.'):
            if number != str(alt_index):
                continue
        elif allele_indices.get(values[0]) != {alt_index}:
            continue
        candidates.append(values)
    if not candidates:
        return "unknown", ""

    # A locus can overlap transcripts from different genes.  PCGR/VEP may mark
    # more than one transcript as PICK=1 (e.g. an upstream/downstream CSDE1
    # transcript and a coding NRAS transcript at NRAS Q61).  Do not therefore
    # take the first PICK entry: prioritise an actual protein consequence and
    # then biological impact, using PICK only as a tie breaker.
    impact_rank = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}

    def annotation_rank(values: list[str]) -> tuple[int, int, int]:
        protein = values[protein_index] if protein_index is not None and len(values) > protein_index else ""
        has_protein_change = bool(re.search(r"(?:^|:)p\.[^:]+$", protein))
        impact = values[impact_index] if impact_index is not None and len(values) > impact_index else ""
        is_pick = pick_index is not None and len(values) > pick_index and values[pick_index] == "1"
        return (
            0 if has_protein_change else 1,
            impact_rank.get(impact, len(impact_rank)),
            0 if is_pick else 1,
        )

    values = min(candidates, key=annotation_rank)
    protein = values[protein_index] if protein_index is not None and len(values) > protein_index else ""
    # PCGR writes HGVSp as e.g. ``ENSP00000286719.6:p.Pro349=``.  The
    # accession prefix is not needed by PhylogicNDT; synonymous ``p.=``
    # annotations do not represent a protein change and are left blank.
    match = re.search(r"(p\.[^:]+)$", protein)
    protein = match.group(1) if match else ""
    if protein.endswith("="):
        protein = ""
    return values[symbol_index], protein


def parse_vcf(vcf_path: Path, sample_id: str, segments, starts, writer: csv.DictWriter, vcf_sample=None) -> Counter:
    stats: Counter = Counter()
    tumour_index: int | None = None
    csq_fields: list[str] | None = None
    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=PCGR_CSQ"):
                match = re.search(r'Format: ([^"]+)', line)
                if match:
                    csq_fields = match.group(1).split("|")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                expected = vcf_sample or sample_id
                tumour_index = samples.index(expected) if expected in samples else None
                if tumour_index is None:
                    raise ValueError(f"Cannot select explicitly named tumour sample {expected} in {vcf_path}")
                continue
            if tumour_index is None:
                raise ValueError(f"VCF has no #CHROM header: {vcf_path}")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or fields[6] != "PASS":
                continue
            stats["pass_records"] += 1
            chrom = canonical_chromosome(fields[0])
            if chrom not in STANDARD_CHROMS:
                stats["nonstandard_chromosome"] += 1
                continue
            try:
                position = int(fields[1])
            except ValueError:
                stats["bad_position"] += 1
                continue
            cn = local_cn(position, segments.get(chrom, []), starts.get(chrom, []))
            if cn is None:
                stats["no_cn_segment"] += 1
                continue
            minor, major = cn
            if minor + major <= 0:
                stats["homozygous_deletion"] += 1
                continue
            format_keys = fields[8].split(":")
            tumour_values = fields[9 + tumour_index].split(":")
            genotype = dict(zip(format_keys, tumour_values))
            try:
                ad = [int(x) for x in genotype["AD"].split(",")]
            except (KeyError, ValueError):
                stats["missing_ad"] += 1
                continue
            if len(ad) != len(fields[4].split(",")) + 1:
                stats["malformed_ad"] += 1
                continue
            if any(count < 0 for count in ad):
                stats["bad_ad"] += 1
                continue
            depth = sum(ad)
            ref_count = ad[0]
            for alt_index, alt in enumerate(fields[4].split(","), start=1):
                alt_count = ad[alt_index]
                if not re.fullmatch("[ACGT]+", fields[3], re.I) or not re.fullmatch("[ACGT]+", alt, re.I):
                    stats["unsupported_alleles"] += 1
                    continue
                if depth == 0:
                    stats["zero_depth"] += 1
                    continue
                classification, kind = variant_type(fields[3], alt)
                gene, protein = gene_annotation(fields[7], alt, csq_fields,
                                                fields[3], alt_index, fields[4].split(','))
                stats["annotated_gene" if gene != "unknown" else "unannotated_gene"] += 1
                writer.writerow({
                    "Hugo_Symbol": gene, "Chromosome": chrom,
                    "Start_position": position, "End_position": position + len(fields[3]) - 1,
                    "Reference_Allele": fields[3], "Tumor_Seq_Allele2": alt,
                    "Variant_Classification": classification, "Variant_Type": kind,
                    "Protein_Change": protein,
                    # The core's binomial model needs all non-target-ALT reads.
                    # Preserve literal VCF REF/other-ALT counts separately.
                    "Tumor_Sample_Barcode": sample_id, "t_ref_count": depth - alt_count,
                    "t_alt_count": alt_count, "local_cn_a1": f"{minor:.6g}",
                    "local_cn_a2": f"{major:.6g}",
                    "vcf_ref_count": ref_count, "vcf_other_alt_count": depth - ref_count - alt_count,
                })
                stats["written_variants"] += 1
    return stats


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, fields, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def prepare(manifest, outdir, results_root, input_root=None):
    rows = read_tsv(manifest)
    required = {'sample_id', 'vcf_path', 'cn_path', 'purity_path'}
    if not rows or not required <= rows[0].keys():
        raise ValueError('Manifest requires sample_id, vcf_path, cn_path, purity_path')
    names = [row['sample_id'] for row in rows]
    if len(names) != len(set(names)) or any(not re.fullmatch(r'[A-Za-z0-9_.-]+', n) or n in {'.', '..'} for n in names):
        raise ValueError('Sample IDs must be unique safe basenames')
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError('Output directory must be new or empty; preserve existing preparations')
    prepared = []
    # Validate the explicit specimen-to-file mapping before writing any outputs.
    for row in rows:
        item = dict(row)
        for field in ['vcf_path', 'cn_path', 'purity_path']:
            p = Path(row[field]).expanduser()
            item[field] = (manifest.parent / p).resolve() if not p.is_absolute() else p.resolve()
            if not item[field].is_file():
                raise FileNotFoundError(item[field])
        purities = read_tsv(item['purity_path'])
        if len(purities) != 1:
            raise ValueError('Expected exactly one PURPLE purity row: '+str(item['purity_path']))
        purity_row = purities[0]
        item['purity'], item['ploidy'] = float(purity_row['purity']), float(purity_row['ploidy'])
        if not 0 < item['purity'] <= 1 or not math.isfinite(item['ploidy']) or item['ploidy'] <= 0:
            raise ValueError('Invalid PURPLE purity/ploidy for '+item['sample_id'])
        purple_sex = purity_row.get('gender', '').upper()
        sex = row.get('sex', '').upper() or purple_sex
        if sex not in {'MALE', 'FEMALE'} or purple_sex in {'MALE', 'FEMALE'} and purple_sex != sex:
            raise ValueError('Require validated MALE/FEMALE sex concordant with PURPLE: '+item['sample_id'])
        aberrations = purity_row.get('germlineAberrations', 'NONE').strip().upper()
        if aberrations not in {'', 'NONE', '[]'}:
            raise ValueError('Non-standard germline CN/sex karyotype requires individual review: '+item['sample_id'])
        item['sex'] = sex
        item['par_mode'] = row.get('par_mode', '') or 'unresolved'
        if item['par_mode'] not in {'unresolved', 'x_diploid_y_masked', 'separate_haploid'}:
            raise ValueError('Unknown par_mode: '+item['par_mode'])
        item['purple_wgd'] = purity_row.get('wholeGenomeDuplication', '')
        item['timepoint'] = float(row.get('timepoint') or 0)
        if not math.isfinite(item['timepoint']):
            raise ValueError('Non-finite timepoint')
        prepared.append(item)
    for folder in ['maf_with_cn', 'segments', 'cluster_sif', 'timing_sif']:
        (outdir/folder).mkdir(parents=True, exist_ok=True)
    input_root = PurePosixPath(input_root) if input_root else outdir.resolve()
    results_root = PurePosixPath(results_root)
    sif_fields = ['sample_id', 'maf_fn', 'seg_fn', 'purity', 'timepoint', 'ploidy', 'sex', 'par_mode']
    qc, provenance = [], []
    for row in prepared:
        sample_id = row['sample_id']
        clipped_segments = []
        segments, starts = read_segments(row['cn_path'], clipped_segments)
        if not segments:
            raise ValueError('No usable PURPLE CN segments for '+sample_id)
        for chrom, values in segments.items():
            if any(a[1] >= b[0] for a, b in zip(values, values[1:])):
                raise ValueError('Overlapping CN segments: '+sample_id+' '+chrom)
        segment_path = outdir/'segments'/f'{sample_id}.alleliccapseg.tsv'
        with segment_path.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['Chromosome','Start.bp','End.bp','mu.minor','sigma.minor','mu.major','sigma.major'])
            for chrom, values in segments.items():
                for start, end, minor, major, baf_count, depth_windows in values:
                    a = alleliccapseg_values(minor, row['purity'], row['ploidy'], baf_count, depth_windows)
                    b = alleliccapseg_values(major, row['purity'], row['ploidy'], baf_count, depth_windows)
                    # PURPLE is one-based closed; core intervals are half-open.
                    writer.writerow([chrom, start, end + 1, *(f'{v:.8g}' for v in (*a, *b))])
        maf_path = outdir/'maf_with_cn'/f'{sample_id}.withCN.maf'
        with maf_path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, MAF_COLUMNS, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            stats = parse_vcf(row['vcf_path'], sample_id, segments, starts, writer, row.get('vcf_sample'))
        if not stats['written_variants']:
            raise ValueError('No usable variants for '+sample_id+'; inspect VCF tumour name, AD and CN mapping')
        sif = {key: row[key] for key in ['sample_id','purity','ploidy','timepoint','sex','par_mode']}
        sif.update(maf_fn=str(input_root/'maf_with_cn'/maf_path.name),
                   seg_fn=str(input_root/'segments'/segment_path.name))
        write_tsv(outdir/'cluster_sif'/f'{sample_id}.sif', sif_fields, [sif])
        timing = dict(sif, maf_fn=str(results_root/sample_id/f'{sample_id}.mut_ccfs.txt'),
                      cluster_ccf_fn=str(results_root/sample_id/f'{sample_id}.cluster_ccfs.txt'))
        write_tsv(outdir/'timing_sif'/f'{sample_id}.sif', sif_fields+['cluster_ccf_fn'], [timing])
        qc.append(dict(sample_id=sample_id, sigma_method='proxy_support',
                       negative_cn_segments_clipped=len(clipped_segments), **stats))
        provenance.append(dict(sample_id=sample_id, vcf_sample=row.get('vcf_sample') or sample_id,
                               purity=row['purity'], ploidy=row['ploidy'], sex=row['sex'], par_mode=row['par_mode'],
                               purple_wgd=row['purple_wgd'], negative_cn_clipped_segments=clipped_segments,
                               inputs={k:dict(path=str(row[k]),sha256=sha256(row[k])) for k in ['vcf_path','cn_path','purity_path']}))
    qc_fields = ['sample_id','sigma_method']+sorted(set().union(*(r.keys() for r in qc))-{'sample_id','sigma_method'})
    write_tsv(outdir/'preparation_qc.tsv', qc_fields, qc)
    (outdir/'samples.txt').write_text(''.join(n+'\n' for n in names))
    (outdir/'input_provenance.json').write_text(json.dumps(dict(genome_build='hg38',
        segment_coordinates='one_based_half_open', count_contract='t_ref_count=sum(AD)-target_ALT; vcf_ref_count=AD[0]',
        samples=provenance), indent=2)+'\n')
    print(f'Prepared {len(names)} specimens; one independent SIF per specimen in {outdir}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--results-root', required=True, help='Absolute runtime path for per-sample result directories')
    parser.add_argument('--input-root', help='Absolute runtime input mount, e.g. /input; default: absolute outdir')
    args = parser.parse_args()
    if not PurePosixPath(args.results_root).is_absolute() or args.input_root and not PurePosixPath(args.input_root).is_absolute():
        parser.error('Runtime roots must be absolute')
    prepare(args.manifest.resolve(), args.outdir.resolve(), args.results_root, args.input_root)


if __name__ == '__main__':
    main()
