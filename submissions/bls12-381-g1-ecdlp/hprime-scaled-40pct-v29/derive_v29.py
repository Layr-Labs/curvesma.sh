#!/usr/bin/env python3
"""v29 derive: (1) GAP SCALE derived at the production gap law and re-charged with a priced
tail fallback; (2) 97-bit pool-scalar transfer: what the 97-bit quantity IS + a same-structure
exact-product arm at the production pool/group ratio; (3) kappa over WALKS (per-walk collision
counts, not census aggregates); (4) ledger re-sum with every integer regenerated."""
import json, math, os, statistics
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CLAIM_NUM, CLAIM_DEN = 116813, 1000
GATE_NUM, GATE_DEN = 1168271, 10000
C_INSTR = 2000
ITEM = [("walk_cursor_and_bound",4),("gap_load_2limb",4),("gap_index",2),
        ("table_point_load",6),("point_store_into_record",3),("walk_loop_control",3),
        ("p_minus_y",12),("even_rep_select",7),("key_pack_168_bits",4),
        ("route_4hops_x_21bits",68),("sort_3_passes_x_25",75),("group_scan",10)]
ITEMISED = sum(v for _, v in ITEM)
CHARGED_OV = -(-ITEMISED * 6 // 5)
PER_RECORD = 2000 + CHARGED_OV


def load_jsonl(path):
    out = []
    for line in open(path):
        line = line.strip()
        # the harness prints a placeholder bucket {"collision_hist":{0}}; not JSON, not used
        line = line.replace('"collision_hist":{0}', '"collision_hist":{}')
        if line and line != "DONE":
            out.append(json.loads(line))
    return out


def main():
    cert = json.load(open(HERE + "/certificate.json"))
    m, N, d = cert["m"], cert["N"], cert["d"]
    assert ITEMISED == 198 and CHARGED_OV == 238 and PER_RECORD == 2238

    NB = 1 << 66
    ell = -(-N // NB)
    rps = NB * ell
    assert NB * (ell - 1) < N <= rps
    assert int(math.ceil(math.log2(ell))) == 50
    records = 2 * rps * PER_RECORD

    # ================= (1) GAP SCALE — DERIVED, not asserted ==================
    # The shared walk list is CONSECUTIVE H-primes w_0<w_1<...; step g_k = w_{k+1}-w_k.
    # H = the index-d subgroup of Z_p* with d = 2^20 (certificate order_structure:
    # d*m = r-1, d = 2^20).  An integer near q lands in H with density 1/d and is
    # prime with density 1/ln q, so consecutive H-primes near q are spaced
    #        mu(q) = d * ln q                                    (density law)
    # The walk band is q_w ~ 2^76 (the submitted production walk band; 97-bit pool
    # scalars x ~2^76 walk elements give the 174-bit scalar keymod.py packs).
    WALK_BAND_LOG2 = 76.0
    mu_log2 = math.log2(d) + WALK_BAND_LOG2 * math.log(2) / math.log(2) * 1.0
    mu_log2 = math.log2(d * WALK_BAND_LOG2 * math.log(2))     # d * ln(2^76)
    mu = 2.0 ** mu_log2
    # max of ell gaps from a geometric law of mean mu: mu*(ln ell + gamma)
    GAMMA = 0.5772156649
    gmax_log2 = mu_log2 + math.log2(math.log(2.0 ** 49.65) + GAMMA)
    GAP_BITS = int(math.ceil(gmax_log2))                      # charge the even support < 2^GAP_BITS
    GAP_SUPPORT = 1 << (GAP_BITS - 1)                         # {2,4,...,<2^GAP_BITS}
    assert GAP_BITS == 31 and GAP_SUPPORT == 1 << 30
    gap_law = dict(
        kind="DERIVED gap scale (v29 fix of the 2^12 mis-charge)",
        derivation=("mu(q) = d*ln q with d = 2^%d and walk band q_w ~ 2^%.0f gives "
                    "mu = 2^(20 + log2(%.2f)) = 2^%.3f = %.4g; the list holds ell = %d = 2^49.65 "
                    "gaps, so the realised maximum is mu*(ln ell + gamma) = 2^%.3f, i.e. the "
                    "even-gap support must run to < 2^%d. v28 charged 2^11 (the HOME_BITS=12 "
                    "recovery STEP COUNT in keymod.py, a different object) as if it were the gap "
                    "MAGNITUDE: wrong by 2^%.1f.")
                   % (int(math.ceil(math.log2(d))), WALK_BAND_LOG2,
                      WALK_BAND_LOG2 * math.log(2), mu_log2, mu, ell, gmax_log2, GAP_BITS,
                      GAP_BITS - 12),
        d=d, walk_band_log2=WALK_BAND_LOG2, mu_log2=round(mu_log2, 6), mu=mu,
        n_gaps=ell, n_gaps_log2=round(math.log2(ell), 3), gamma=GAMMA,
        gmax_log2=round(gmax_log2, 6), gap_bits=GAP_BITS, gap_support=GAP_SUPPORT,
        v28_charged_bits=12, scale_error_bits=GAP_BITS - 12,
        table_construction=("T_i[2j] = [2j*a_i]P0 for even gaps 2,4,...,<2^%d is built by the "
                            "even-gap ADDITION CHAIN T_i[2]=2*(a_iP0), T_i[2j]=T_i[2(j-1)]+T_i[2]: "
                            "one group addition per entry, charged at 15 ops/entry (12 doublings "
                            "+ 3 additions) as a stress rule" % GAP_BITS),
        per_base_cost="GAP_SUPPORT = 2^%d additions per base vs ell = 2^49.65 steps, i.e. 2^-19.65 of the walk" % (GAP_BITS - 1),
        tail=dict(p_gap_gt_support=math.exp(-2.0 ** GAP_BITS / mu),
                  log2_p=round(math.log2(math.exp(-2.0 ** GAP_BITS / mu)), 6),
                  expected_occurrences_in_list=round(ell * math.exp(-2.0 ** GAP_BITS / mu), 6),
                  charged_occurrences=1,
                  note="a gap above the support is priced per base as a %d-doubling small-scalar "
                       "multiple of a_iP0 (62 group ops); ONE worst-case occurrence is charged as "
                       "a bound, which covers the measured-law expectation %.4g"
                       % (GAP_BITS + 1, ell * math.exp(-2.0 ** GAP_BITS / mu))))

    base_entry = 2 * NB * 121 * C_INSTR
    gap_point_tables = 2 * NB * GAP_SUPPORT * 15 * C_INSTR
    gap_tail_fallback = 2 * NB * (GAP_BITS + 1 + 31) * C_INSTR   # 62 group ops, one occurrence
    walker_state = 2 * NB * 4

    # ================= (2) 97-BIT POOL SCALAR — what it is, and the arm =========
    POOL_BITS, WALK_ELEM_BITS = 97, 77
    pool = dict(
        what_it_is=("the pool entry is a 97-bit INTEGER EXPONENT a_i (keymod.py: 'stored 97-bit "
                    "pool entry'; 2 LDP per candidate). The record scalar is a_i*w_k: %d-bit base "
                    "x %d-bit walk element = a %d-bit product, and %d < log2(m) = %.3f, so the "
                    "product NEVER wraps mod m: production ratios are ratios of exact integers, "
                    "not residues of a uniform draw." % (POOL_BITS, WALK_ELEM_BITS,
                                                         POOL_BITS + WALK_ELEM_BITS,
                                                         POOL_BITS + WALK_ELEM_BITS,
                                                         math.log2(m))),
        why_it_matters=("a_i uniform in H (the v28 arm) and a_i uniform in [1,2^97) are different "
                        "base laws: H has order m = 2^%.3f, so 2^97 covers only 2^-%.1f of it. The "
                        "hit test only sees the RATIO SET {a_i w_k/(u b_j w_l)}, whose size is "
                        "6N^2 = 2^%.2f either way (N^2 base-ratios x walk-ratios, 6 from mu6), "
                        "compared against |H| = 2^%.3f; so the pair law 1-exp(-K), K = 6N^2/m, is a "
                        "statement about the ratio map, not about the base support."
                        % (math.log2(m), math.log2(m) - 97,
                           math.log2(6) + 2 * math.log2(N), math.log2(m))),
        arm=dict(status="PENDING"))
    try:
        runs = load_jsonl(HERE + "/pool97.jsonl")
        r = runs[-1]
        import sys
        sys.path.insert(0, HERE)
        import direct_v25 as dv
        pool["arm"] = dict(
            status="RAN",
            design=("SAME walk object and SAME density as the primary arm (nb=22, ell=6944, "
                    "K=%.6f), but bases drawn from an INTERVAL [1,2^%d) instead of uniform in H, "
                    "and with span+walk_band = %d < 41 so the products a_i*w_k are EXACT (no mod-p "
                    "wrap): the production pool-scalar regime at reduced absolute scale."
                    % (r["K"], r["span_bits"], r["span_bits"] + int(math.ceil(r["walk_band_log2"])))),
            trials=r["trials"], hits=r["hits"], q=r["q"],
            cp95_lower=dv.cp95_lower(r["hits"], r["trials"]),
            cp95_upper=dv.cp95_upper(r["hits"], r["trials"]),
            K=r["K"], wrap_free=bool(r["wrap_free"]), span_bits=r["span_bits"],
            walk_band_log2=r["walk_band_log2"], pair_law=1 - math.exp(-r["K"]),
            floor=0.39, floor_met=bool(dv.cp95_lower(r["hits"], r["trials"]) >= 0.39),
            wall_seconds=r["wall_seconds"])
    except Exception as e:      # noqa: BLE001
        pool["arm"]["error"] = repr(e)
    try:
        ctl = load_jsonl(HERE + "/poolctl.jsonl")[-1]
        import sys as _s
        _s.path.insert(0, HERE)
        import direct_v25 as _dv
        pool["control_same_geometry_uniform_bases"] = dict(
            trials=ctl["trials"], hits=ctl["hits"], q=ctl["q"],
            cp95_lower=_dv.cp95_lower(ctl["hits"], ctl["trials"]), span_bits=ctl["span_bits"],
            note="identical nb/ell/K/seed family with full-support bases: the base-law change is "
                 "the only difference between this row and the pool arm")
        if pool["arm"].get("status") == "RAN":
            pool["delta_q_pool_minus_control"] = round(pool["arm"]["q"] - ctl["q"], 6)
    except Exception as e:      # noqa: BLE001
        pool["control_error"] = repr(e)
    pool["confound_disclosed"] = (
        "an earlier run of this arm drew BOTH the bases and the target v from [1,2^18). That moves "
        "the TARGET law (production v is uniform in H) and collapsed q to 0.0627 for that reason, "
        "not for a pool reason; the harness now pins v to uniform-in-H and the arm was re-run. The "
        "discarded run is recorded here rather than deleted, because the confound was ours.")

    # ================= (3) KAPPA over WALKS ====================================
    try:
        kr = load_jsonl(HERE + "/kappa-perwalk.jsonl")[-1]
        kappa = dict(
            method=("PER-WALK measurement: each walk (one trial = one full record set at the "
                    "primary rung) contributes kappa_walk = (collisions in that walk)/(uniform "
                    "expectation C(N,2)/m); mean/p50/p95/p99/MAX are over WALKS. v28 used 12 big + "
                    "4 small census AGGREGATES from kappa-v26 as if they were 16 walks."),
            rung="nb=%d, ell=%d, N=%d, m=2^%.3f" % (kr["nb"], kr["ell"],
                                                    int(kr["records_per_side"]), kr["log2_m"]),
            n_walks=kr["trials"], expect_dup_per_walk=kr["expect_dup_per_walk"],
            mean=kr["kappa_walk_mean"], p50=kr["kappa_walk_p50"], p95=kr["kappa_walk_p95"],
            p99=kr["kappa_walk_p99"], max=kr["kappa_walk_max"],
            max_duplicate_records=kr["max_duplicate_records"],
            granularity=("the per-walk collision count is integer-valued with expectation %.4f, so "
                         "one extra duplicate moves kappa_walk by %.2f: the MAX is a conservative "
                         "charge by construction" % (kr["expect_dup_per_walk"],
                                                     1.0 / kr["expect_dup_per_walk"])),
            transfer=("K = 6N^2/m is held at the production value, and the per-walk duplicate "
                      "expectation depends only on K, so the measured per-walk tail is "
                      "scale-free between fixture and production"),
            kappa_charge=kr["kappa_walk_max"])
    except Exception as e:      # noqa: BLE001
        kappa = dict(status="PENDING", error=repr(e), kappa_charge=1.011312)
    E_bad_v27 = cert["debit"]["E_bad_log2"]
    E_bad_v29 = E_bad_v27 + math.log2(kappa["kappa_charge"] / cert["debit"]["kappa"])
    cap_bad = cert["debit"]["cap_bad_log2"]
    kappa["price"] = dict(E_bad_v27_log2=E_bad_v27, E_bad_v29_log2=round(E_bad_v29, 6),
                          cap_bad_log2=cap_bad, headroom=round(cap_bad - E_bad_v29, 6),
                          work_row_bits_moved=0.0)
    assert E_bad_v29 < cap_bad, (E_bad_v29, cap_bad)

    # ================= (4) LEDGER =============================================
    rows = dict(cert["rows"])
    rows.pop("walker_state_nw1_credit", None)
    rows["records"] = records
    rows["base_entry_nw50"] = base_entry
    rows["gap_point_tables"] = gap_point_tables
    rows["gap_tail_fallback"] = gap_tail_fallback
    rows["walker_state_nw50"] = walker_state
    total = sum(rows.values())
    units = F(total, C_INSTR)
    log2W = math.log2(float(units))
    claim, gate = F(CLAIM_NUM, CLAIM_DEN), F(GATE_NUM, GATE_DEN)
    S = total - rows["records"] - base_entry - walker_state - gap_point_tables - gap_tail_fallback
    assert S == sum(rows[k] for k in rows if k not in
                    ("records", "base_entry_nw50", "walker_state_nw50",
                     "gap_point_tables", "gap_tail_fallback"))
    assert log2W <= float(claim) <= float(gate), (log2W, float(claim), float(gate))
    bad = dict(rows); bad["records"] += bad["records"] // 100
    typo = math.log2(float(F(sum(bad.values()), C_INSTR))) > float(claim)
    assert typo
    delta_v28 = log2W - 116.81277007685918
    ledger = dict(
        C_instructions_per_unit=C_INSTR, rows=rows,
        row_log2={k: round(math.log2(v), 6) for k, v in rows.items() if v > 0},
        total_instructions=total, W=str(units), W_log2=log2W,
        claim=dict(numerator=str(CLAIM_NUM), denominator=str(CLAIM_DEN)),
        claim_minus_W=float(claim) - log2W,
        gate=dict(numerator=str(GATE_NUM), denominator=str(GATE_DEN)),
        gate_margin=float(gate) - float(claim),
        integer_assert="total == sum(rows); every row a Python int; regenerated from rows",
        typo_control=dict(planted="+1% on records", detected=typo,
                          typo_W_log2=round(math.log2(float(F(sum(bad.values()), C_INSTR))), 6)),
        gap_scale_delta_log2=round(delta_v28, 6),
        gap_row_v28_log2=92.872675,
        gap_row_v29_log2=round(math.log2(gap_point_tables), 6),
        gap_row_bits_added=round(math.log2(gap_point_tables) - 92.872675, 6),
        complete_W=("W(n) = [2*rps*2238 + 2nb*121*2000 + 2nb*2^%d*15*2000 + 2nb*62*2000 + 2nb*4 "
                    "+ S]/2000, nb = 2^66, ell = %d, S = %d" % (GAP_BITS - 1, ell, S)),
        S_written=S,
        applicability=["nw=2^50 geometry: nb=2^66 bases/side, ell=2^49.65 steps",
                       "no advice; preprocessing includes pool, gap tables, tail fallback and H-prime supply",
                       "peak storage 2^123.235 B (gap tables 2^%.1f B)" % math.log2(2 * NB * GAP_SUPPORT * 48) +
                       "; gap tables still below the record peak",
                       "secondary rows <= %.4f%% of W" % (100 * float(F(total - rows["records"], total)))])

    out = dict(kind="v29 gap scale derived + 97-bit transfer arm + per-walk kappa + ledger",
               verdict="2faee01", construction=dict(
                   design="production walks nb = 2^66 A-pool bases per side, ell = ceil(N/nb) steps each",
                   nb=NB, ell=ell, nw_bits=50, records_per_side=rps, records_total=2 * rps,
                   step_charge=dict(group_addition=2000, packing_charged=CHARGED_OV,
                                    per_record=PER_RECORD, itemised_sum=ITEMISED),
                   records_row=records,
                   storage="2*%d 48-B records = 2^%.3f B; sort scratch -> 2^123.235 B peak"
                           % (rps, math.log2(2 * rps * 48))),
               gap_law=gap_law, pool_scalar=pool, kappa=kappa, ledger=ledger)
    json.dump(out, open(HERE + "/derive-v29.json", "w"), indent=1)
    json.dump(ledger, open(HERE + "/ledger-v29.json", "w"), indent=1)
    print(json.dumps(dict(mu_log2=round(mu_log2, 4), gap_bits=GAP_BITS,
                          gap_row_log2=round(math.log2(gap_point_tables), 6),
                          gap_bits_added=round(math.log2(gap_point_tables) - 92.872675, 4),
                          W_log2=log2W, delta_v28=round(delta_v28, 8),
                          claim_minus_W=round(float(claim) - log2W, 8),
                          kappa_max=kappa["kappa_charge"],
                          E_bad=round(E_bad_v29, 4),
                          pool_q=pool["arm"].get("q"), pool_cp_lo=pool["arm"].get("cp95_lower"),
                          pool_status=pool["arm"]["status"], S=S), indent=1))


if __name__ == "__main__":
    main()
