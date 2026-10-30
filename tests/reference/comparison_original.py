"""Historical comparison reference for regression testing; upstream LICENSE applies."""
def compare_events_original(timing_engine, drivers=()):
    """Original PhylogicNDT pi-distribution difference comparison method."""
    all_events = []
    if timing_engine.WGD:
        all_events.append(timing_engine.WGD)
    all_events.extend(_dedup_cn_events(timing_engine.all_cn_events))
    all_events.extend(
        mut for mut in timing_engine.mutations.values()
        if mut.prot_change
        and mut.gene in drivers
        and (mut.prot_change[0] != mut.prot_change[-1] or not mut.prot_change[-1].isalpha())
    )

    comps = {}
    for eve1, eve2 in itertools.combinations(all_events, 2):
        if eve1.pi_dist is None or eve2.pi_dist is None:
            continue
        pi_diff = eve1.pi_dist - eve2.pi_dist
        p_before = 0.0
        p_after = 0.0
        diff_sign = 0.0
        for d in pi_diff:
            if d > 0:
                if diff_sign < 0:
                    break
                p_before += d
                diff_sign = 1
            if d < 0:
                if diff_sign > 0:
                    break
                p_after -= d
                diff_sign = -1
        for d in reversed(pi_diff):
            if d > 0:
                if diff_sign < 0:
                    break
                p_after += d
                diff_sign = 1
            if d < 0:
                if diff_sign > 0:
                    break
                p_before -= d
                diff_sign = -1
        p_unknown = 1.0 - p_before - p_after
        comps[(eve1.event_name, eve2.event_name)] = (p_before, p_after, p_unknown)
    return comps
