"""Normal copy state for the hg38 PURPLE absolute-CN input contract.

Coordinates are one-based, inclusive. PAR mode must describe the alignment
reference, not be inferred from a filename containing 'masked'.
"""

PAR38 = {'X': ((10001, 2781479), (155701383, 156030895)),
         'Y': ((10001, 2781479), (56887903, 57217415))}


def chromosome(chrom):
    value = str(chrom)
    return value[3:] if value.startswith('chr') else value


def context(sample):
    if hasattr(sample, 'sex'):
        return sample
    if hasattr(sample, 'sample'):
        return context(sample.sample)
    if getattr(sample, 'sample_list', None):
        return context(sample.sample_list[0])
    return sample


def baseline(sample, chrom, pos=None):
    chrom = chromosome(chrom)
    if chrom not in ('X', 'Y'):
        return (1., 1.)
    sample = context(sample)
    sex = getattr(sample, 'sex', None)
    if sex not in ('MALE', 'FEMALE'):
        raise ValueError('X/Y requires validated MALE/FEMALE sex metadata')
    if sex == 'FEMALE':
        return (1., 1.) if chrom == 'X' else (0., 0.)
    if pos is not None and any(a <= int(pos) <= b for a, b in PAR38[chrom]):
        mode = getattr(sample, 'par_mode', 'unresolved')
        if mode == 'x_diploid_y_masked':
            return (1., 1.) if chrom == 'X' else None
        if mode == 'separate_haploid':
            return (0., 1.)
        return None
    return (0., 1.)


def normal_cn(sample, chrom, pos):
    state = baseline(sample, chrom, pos)
    return None if state is None else sum(state)


def region_baseline(sample, chrom, start, end):
    """Return a single valid baseline for a half-open interval, or unknown."""
    points = [int(start), int(end) - 1]
    for a, b in PAR38.get(chromosome(chrom), ()):
        if start < b + 1 and end > a:
            points.extend([max(start, a), min(end - 1, b)])
    states = [baseline(sample, chrom, p) for p in points]
    return states[0] if states[0] is not None and all(s == states[0] for s in states) else None


def write_qc(samples, phase):
    import csv
    import json
    for sample in samples:
        with open(sample.sample_name + '.' + phase + '.sex_cn_qc.tsv', 'w') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['sample', 'sex', 'par_mode', 'record', 'reason'])
            for record, reason in sample.sex_cn_unresolved:
                writer.writerow([sample.sample_name, sample.sex, sample.par_mode, json.dumps(record), reason])
