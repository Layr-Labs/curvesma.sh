#!/usr/bin/env python3
"""v28 structured-walk ladder evidence.

The PRIMARY arm is the 120,000-trial run at the largest affordable structured walk
(nb=22, ell=6944 = 108x the fixture's own walk length of 64) at production density.
The ladder rungs are 20,000-trial runs at the same density with only the nb x ell split
changed, so the split dependence is measured rather than assumed.
"""
import json, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import direct_v25 as dv

FIXTURE_ELL = 64.0
PROD_K = 0.5094


def load(path):
    out = []
    for line in open(path):
        line = line.strip()
        if not line or line == "DONE":
            continue
        r = json.loads(line)
        n, k = r["trials"], r["hits"]
        r.update(q=k / n, cp95_lower=dv.cp95_lower(k, n), cp95_upper=dv.cp95_upper(k, n),
                 pair_law=1 - math.exp(-r["K"]), floor=0.39,
                 floor_met=bool(dv.cp95_lower(k, n) >= 0.39),
                 nb_over_ell=round(r["nb"] / r["ell"], 6),
                 ell_over_fixture=round(r["ell"] / FIXTURE_ELL, 3))
        out.append(r)
    return out


def main():
    rungs = load("ladder2.jsonl")
    prim_raw = json.loads(open("primary120k.json").read().strip().splitlines()[0])
    n, k = prim_raw["trials"], prim_raw["hits"]
    prim = dict(prim_raw, q=k / n, cp95_lower=dv.cp95_lower(k, n), cp95_upper=dv.cp95_upper(k, n),
                pair_law=1 - math.exp(-prim_raw["K"]), floor=0.39,
                floor_met=bool(dv.cp95_lower(k, n) >= 0.39),
                nb_over_ell=round(prim_raw["nb"] / prim_raw["ell"], 6),
                ell_over_fixture=round(prim_raw["ell"] / FIXTURE_ELL, 3),
                note="PRIMARY: largest affordable structured walk at production density")
    Ks = {round(r["K"], 9) for r in rungs} | {round(prim["K"], 9)}
    assert len(Ks) == 1, Ks
    K = Ks.pop()
    rungs_sorted = sorted(rungs, key=lambda r: -r["ell"])
    long_walk = [r for r in rungs_sorted if r["ell"] >= FIXTURE_ELL]
    short_walk = [r for r in rungs_sorted if r["ell"] < FIXTURE_ELL]
    pool_n = sum(r["trials"] for r in long_walk) + prim["trials"]
    pool_k = sum(r["hits"] for r in long_walk) + prim["hits"]
    all_n = pool_n + sum(r["trials"] for r in short_walk)
    all_k = pool_k + sum(r["hits"] for r in short_walk)

    out = dict(
        kind="v28 STRUCTURED-WALK arm (production object: shared walk list, base x walk "
             "products, mu6 test) with an nb x ell ladder at fixed density",
        field_p=prim["p"], m=prim["m"], log2_m=prim["log2_m"],
        records_per_side=prim["records_per_side"], K=K, K_production=PROD_K,
        density_match="K = 6N^2/m held at every rung and on the primary arm; only nb x ell varies",
        design=dict(
            walk_list="consecutive H-primes (prime integers in the index-8 subgroup of order "
                      "m=2^38), ONE list shared by both sides and by every base",
            records="records = a_i * w_k mod p for base a_i in H and shared walk element w_k",
            target="v uniform in H; hit iff exists u in mu6, i,j,k,l: v*u*b_j*w_l == a_i*w_k",
            why_this_is_the_right_object="v27 sampled independent D-th power residues and tested a "
                "ratio SET. Production records are base x walk products over a shared W with shared "
                "mu6 keys; this harness keeps exactly that structure and the same density.",
            fixture_ell=FIXTURE_ELL),
        primary=dict(rung="nb=%d, ell=%d (%.1fx the fixture walk length)" %
                          (prim["nb"], prim["ell"], prim["ell"] / FIXTURE_ELL),
                     trials=prim["trials"], hits=prim["hits"], q=prim["q"],
                     cp95_lower=prim["cp95_lower"], cp95_upper=prim["cp95_upper"],
                     floor=0.39, floor_met=prim["floor_met"], K=prim["K"],
                     wall_seconds=prim["wall_seconds"]),
        ladder=[{k2: v for k2, v in r.items() if k2 != "workers"} for r in rungs_sorted],
        pooled_long_walks=dict(rungs="ladder rungs with ell >= %.0f, plus the primary arm" % FIXTURE_ELL,
                               trials=pool_n, hits=pool_k, q=pool_k / pool_n,
                               cp95_lower=dv.cp95_lower(pool_k, pool_n),
                               cp95_upper=dv.cp95_upper(pool_k, pool_n),
                               floor_met=bool(dv.cp95_lower(pool_k, pool_n) >= 0.39)),
        pooled_all=dict(trials=all_n, hits=all_k, q=all_k / all_n,
                        cp95_lower=dv.cp95_lower(all_k, all_n),
                        note="includes the ell=16 short-walk rung, DISCLOSED as sub-floor; the "
                             "walk list must be long enough for the ratio set to densify, and "
                             "production's ell=2^49.65 is 2^45 times the largest rung"),
        short_walk_disclosure=[{kk: r[kk] for kk in ("ell", "nb", "trials", "hits", "q",
                                                    "cp95_lower", "floor_met")}
                               for r in short_walk],
        spread=dict(q_min=min(r["q"] for r in rungs_sorted),
                    q_max=max(r["q"] for r in rungs_sorted),
                    rung_se=round(math.sqrt(0.396 * 0.604 / 20000), 6)),
        ratio_arm_v27=dict(label="SECONDARY: independent-residue ratio statistics, NOT the walk "
                                 "distribution", trials=20000, hits=7997, q=0.39985,
                           cp95_lower=0.394138, cp95_upper=0.405583, m_log2=23.000283),
        scope="fixture-scale field (p=2^41): density K, the shared-walk product structure, the "
              "mu6 hit test and the ell ordering are matched; the ABSOLUTE split nb/ell=2^16.35 is "
              "bracketed, not attained. Absolute-scale transfer is the law's job, disclosed.")
    json.dump(out, open("structured-ladder.json", "w"), indent=1)
    print(json.dumps(dict(K=K, primary_q=prim["q"], primary_cp=dv.cp95_lower(k, n),
                          primary_floor=prim["floor_met"],
                          pooled_q=pool_k / pool_n, pooled_cp=dv.cp95_lower(pool_k, pool_n),
                          pooled_floor=bool(dv.cp95_lower(pool_k, pool_n) >= 0.39),
                          rung_qs=[(r["ell"], r["q"]) for r in rungs_sorted]), indent=1))
    return out


if __name__ == "__main__":
    main()
