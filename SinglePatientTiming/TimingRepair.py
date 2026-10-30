"""Small compatibility repairs; retain native Timing and comparison likelihoods."""
import copy
import csv
import re
from collections import defaultdict
import numpy as np


def read_cluster_ccfs(path, sample_names):
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        bins = sorted((k for k in reader.fieldnames if k.startswith('postDP_ccf_')
                       and k not in ('postDP_ccf_mean', 'postDP_ccf_CI_low', 'postDP_ccf_CI_high')),
                      key=lambda k: float(k[11:]))
        if len(bins) != 101:
            raise ValueError('Expected 101 postDP cluster bins')
        rows = {(r['Sample_ID'], int(r['Cluster_ID'])): [float(r[k]) for k in bins] for r in reader}
    clusters = sorted({c for _, c in rows})
    grids = {c: np.asarray([rows[(s, c)] for s in sample_names]) for c in clusters}
    if 1 not in grids or any(not np.isfinite(v).all() or (v.sum(axis=1) <= 0).any() for v in grids.values()):
        raise ValueError('Invalid saved cluster distributions')
    return grids


def nearest_cn_cluster(ccfs, grids):
    """Native ClusterEngine CN mode-distance rule; leave exact ties unresolved."""
    target = np.rint(np.clip(np.asarray(ccfs, dtype=float), 0, 1) * 100) / 100
    if not np.isfinite(target).all():
        return None
    distances = {c: np.abs(np.argmax(g, axis=1) / 100. - target) for c, g in grids.items()}
    best = min(v.sum() for v in distances.values())
    winners = [c for c, v in distances.items() if np.isclose(v.sum(), best, atol=1e-12, rtol=0)]
    if len(winners) != 1 or distances[winners[0]].max() > .5:
        return None
    return winners[0]


AA3 = dict(zip(('Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Ter').split(),
               'ARNDCQEGHILKMFPSTWYV*'))


def protein_change_kind(text):
    """Parse the HGVS/simple protein forms present in MAFs; fail closed if unknown."""
    text = str(text or '').strip()
    if text.startswith('p.'):
        text = text[2:]
    text = text.strip('()')
    m = re.fullmatch(r'([A-Z][a-z]{2}|[A-Z*])(\d+)([A-Z][a-z]{2}|[A-Z*=])', text)
    if m:
        ref, _, alt = m.groups()
        ref, alt = AA3.get(ref, ref), AA3.get(alt, alt)
        return 'synonymous' if alt == '=' or ref == alt else 'nonsynonymous'
    if re.fullmatch(r'(?:[A-Z][a-z]{2}|[A-Z*])\d+(?:_(?:[A-Z][a-z]{2}|[A-Z*])\d+)?(?:delins|del|dup|ins|fs|ext).*', text):
        return 'nonsynonymous'
    if re.fullmatch(r'(?:[A-Z][a-z]{2}|[A-Z*])\d+(?:[A-Z][a-z]{2}|[A-Z*])fs.*', text):
        return 'nonsynonymous'
    return 'unparsed'


def is_driver(mut, drivers):
    return mut.gene in drivers and bool(mut.prot_change) and protein_change_kind(mut.prot_change) != 'synonymous'


def assign_driver_record_ids(mutations, drivers):
    groups = defaultdict(list)
    for mut in mutations:
        if is_driver(mut, drivers):
            groups[mut.protein_event_name].append(mut)
    for name, members in groups.items():
        if len(members) > 1:
            for mut in members:
                mut.event_id_override = name + '__locus=' + mut.var_str


def aggregate_driver_events(mutations, drivers, include_synonymous=False):
    groups = defaultdict(list)
    for mut in mutations:
        if mut.gene in drivers and mut.prot_change and (include_synonymous or is_driver(mut, drivers)):
            groups[mut.protein_event_name].append(mut)
    result = []
    for name, members in groups.items():
        chosen = copy.copy(sorted(members, key=lambda m: m.var_str)[0])
        chosen.event_id_override = name
        chosen.source_record_count = len(members)
        if len(members) > 1 and (not all(eligible(m) for m in members)
                or len({m.is_clonal for m in members}) > 1
                or any(not np.allclose(members[0].pi_dist, m.pi_dist, atol=1e-10, rtol=0) for m in members[1:])):
            chosen.pi_dist = None
            chosen.timing_status = 'not_estimable'
            chosen.timing_source = None
            chosen.timing_reason = 'ambiguous_driver_representation'
        result.append(chosen)
    return result


def eligible(event):
    if getattr(event, 'prot_change', None) and protein_change_kind(event.prot_change) == 'unparsed':
        return False
    p = getattr(event, 'pi_dist', None)
    if p is None or getattr(event, 'timing_status', 'estimated') not in ('estimated', 'subclonal'):
        return False
    p = np.asarray(p, dtype=float)
    return (p.shape == (101,) and np.isfinite(p).all() and (p >= 0).all()
            and p.sum() > 0 and not np.allclose(p / p.sum(), 1. / 101, atol=1e-12, rtol=0))


def dedup_cn_events(events):
    """Keep one biological event for prevalence, preserving all allele details elsewhere.

    Prefer local estimates over unsupported distributions. Conflicting allele
    clonality or distinct eligible distributions have no unique event ordering.
    """
    grouped = defaultdict(list)
    for evs in events.values():
        for event in evs:
            grouped[event.event_name].append(event)
    result = []
    for name, evs in grouped.items():
        valid = [e for e in evs if eligible(e)]
        classes = {getattr(e, 'is_clonal', None) for e in evs}
        conflict = len(classes) > 1 or any(not np.allclose(valid[0].pi_dist, e.pi_dist, atol=1e-10, rtol=0)
                                         for e in valid[1:])
        chosen = copy.copy((valid or evs)[0])
        chosen.aggregate_alleles = len(evs)
        if conflict or not valid:
            chosen.pi_dist = None
            chosen.timing_status = 'not_estimable'
            chosen.timing_source = None
            chosen.timing_reason = 'discordant_alleles' if conflict else 'no_eligible_allele'
        result.append(chosen)
    return result


def write_event_qc(engine, drivers):
    patient = engine.patient
    grids = patient.saved_cluster_ccfs
    with open(patient.indiv_name + '.mutation_cn_provenance.tsv', 'w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample','variant','cn_source','input_minor','input_major','timing_minor','timing_major'])
        for sample in engine.sample_list:
            for mut in sample.mutations:
                original = sample.sample.get_mut_by_varstr(mut.var_str)
                source = 'timing_segment' if sample.CnProfile[mut.chrN][mut.pos] else 'cluster_input_no_timing_segment'
                writer.writerow([sample.sample_name,mut.var_str,source,original.local_cn_a1,original.local_cn_a2,
                                 mut.local_cn_a1[0],mut.local_cn_a2[0]])
    with open(patient.indiv_name + '.cn_cluster_assignments.tsv', 'w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['event', 'allele', 'ccf', 'cluster', 'rule'])
        for event in patient.sample_list[0].low_coverage_mutations.values():
            if getattr(event, 'type', None) != 'CNV':
                continue
            ccfs = [float(s.get_mut_by_varstr(event.var_str).ccf_hat) for s in patient.sample_list]
            cluster = nearest_cn_cluster(ccfs, grids)
            writer.writerow([event.var_str, 'a1' if event.a1 else 'a2', ','.join(map(str, ccfs)),
                             cluster, 'native_nearest_mode' if cluster else 'unresolved'])
    write_comparison_qc(engine, drivers)


def write_comparison_qc(engine, drivers):
    events = dedup_cn_events(engine.all_cn_events)
    events += aggregate_driver_events(engine.mutations.values(), drivers, include_synonymous=True)
    if engine.WGD:
        events.append(engine.WGD)
    fields = ['event', 'timing_status', 'timing_source', 'timing_reason', 'eligible', 'protein_change_kind', 'comparison_included', 'source_record_count']
    with open(engine.patient.indiv_name + '.event_qc.tsv', 'w') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for e in events:
            protein = getattr(e, 'prot_change', None)
            writer.writerow(dict(event=e.event_name, timing_status=getattr(e, 'timing_status', 'estimated'),
                                 timing_source=getattr(e, 'timing_source', ''), timing_reason=getattr(e, 'timing_reason', ''),
                                 eligible=int(eligible(e)), protein_change_kind=protein_change_kind(protein) if protein else '',
                                 comparison_included=int(not protein or is_driver(e, drivers)),
                                 source_record_count=getattr(e, 'source_record_count', getattr(e, 'aggregate_alleles', 1))))
    with open(engine.patient.indiv_name + '.driver_representation_qc.tsv', 'w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['canonical_event','record_id','variant','cluster','timing_status'])
        for mut in engine.mutations.values():
            if is_driver(mut, drivers) and mut.event_name != mut.protein_event_name:
                writer.writerow([mut.protein_event_name,mut.event_name,mut.var_str,mut.cluster_assignment,mut.timing_status])
