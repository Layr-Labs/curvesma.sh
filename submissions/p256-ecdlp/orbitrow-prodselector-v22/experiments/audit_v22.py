"""v22 planted controls -- four, each of which MUST be able to fail.

(a) HAZARD-MATCH CONTROL.  The bridge in the ledger is an assertion, so the audit
    corrupts ONE fixture hazard input (the row-occupancy chi-square z, then the
    per-trail collision mean) and shows the match test FAILS on the corrupted row and
    PASSES on the as-measured rows.  A match test that cannot fail proves nothing.

(b) CROSS-COLLISION CONTROL.  Two walkers are FORCED into the same row sequence on a
    real P-256 twin (R = 64, real curve arithmetic, same secret) by giving them the
    same trail id; the cross-trail detector must fire and the recovery yield must be
    computed from the actual coefficient difference.  A planted collision that the
    detector misses would void the cross-trail accounting.

(c) 128-BIT LOCATOR CONTROL.  Page index 2^64 + 5 must round-trip through the shipped
    locator split (two 64-bit words).  The u64 MUTANT -- the v21 design -- truncates
    and must return a different page, i.e. the v21 entry size really is insufficient.

(d) REGEN 1-ULP CARRY.  The claim must be the ceiling of the exact log2 work at 1e-4
    granularity: claim - 1 ulp must fall BELOW the exact value (so the claim is the
    tightest 4-decimal value that still bounds the work), and the reduced rational must
    equal the claim exactly.
"""
import importlib.util
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")


def load_selector():
    spec = importlib.util.spec_from_file_location(
        "psel", os.path.join(HERE, "production_selector_v22.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def match_test(ps, chi2_tol=4.0, ent_tol=1e-4, lag_tol=4.0, mean_sigmas=3.0):
    """The bridge's arithmetic: do the two functions' measured inputs agree?"""
    P, F = ps["armP_prod"], ps["armP_fix"]
    CP_, CF = ps["armC_prod"], ps["armC_fix"]
    z_diff = abs(P["occupancy"]["chi2_z"] - F["occupancy"]["chi2_z"])
    ent_diff = abs(P["occupancy"]["entropy_over_log_cells"] - F["occupancy"]["entropy_over_log_cells"])
    lag_diff = abs(P["lag1_z"] - F["lag1_z"])
    lam = CP_["lambda_expected"]
    se = math.sqrt(max(lam, 1e-12) / CP_["trails"] + max(lam, 1e-12) / CF["trails"])
    mean_diff = abs(CP_["observed_mean"] - CF["observed_mean"]) / se
    rows = {
        "occupancy_chi2_z_difference": {"value": z_diff, "tolerance": chi2_tol,
                                        "pass": z_diff <= chi2_tol},
        "entropy_difference": {"value": ent_diff, "tolerance": ent_tol,
                               "pass": ent_diff <= ent_tol},
        "lag1_z_difference": {"value": lag_diff, "tolerance": lag_tol,
                              "pass": lag_diff <= lag_tol},
        "per_trail_mean_difference_sigmas": {"value": mean_diff, "tolerance": mean_sigmas,
                                             "pass": mean_diff <= mean_sigmas},
    }
    return all(r["pass"] for r in rows.values()), rows


def control_a(ps):
    ok_real, rows_real = match_test(ps)
    assert ok_real, "as-measured rows must pass the match test"
    out = {"as_measured": {"pass": ok_real, "rows": rows_real}}
    # corrupt ONE fixture input: the occupancy chi-square z moves by +10
    bad = json.loads(json.dumps(ps))
    bad["armP_fix"]["occupancy"]["chi2_z"] += 10.0
    ok1, rows1 = match_test(bad)
    # corrupt a different single input: the per-trail mean moves by 10 sigma
    bad2 = json.loads(json.dumps(ps))
    bad2["armC_fix"]["observed_mean"] = bad2["armC_prod"]["observed_mean"] + 10 * math.sqrt(
        bad2["armC_prod"]["lambda_expected"] / bad2["armC_prod"]["trails"])
    ok2, rows2 = match_test(bad2)
    out["corrupt_occupancy_z"] = {"pass": ok1, "failed_rows": [k for k, v in rows1.items() if not v["pass"]]}
    out["corrupt_per_trail_mean"] = {"pass": ok2, "failed_rows": [k for k, v in rows2.items() if not v["pass"]]}
    assert not ok1 and not ok2, "the match test must FAIL when a hazard input is corrupted"
    out["verdict"] = ("PASS: the as-measured rows satisfy the match test; corrupting either the "
                      "occupancy z or the per-trail mean breaches it, so the bridge is CHECKED, "
                      "not asserted")
    return out


def control_b():
    """Force two real-P-256 walkers into the same row sequence; detector must fire.

    Walker A starts at T = [c_A + d_A*x]P.  Walker B is planted with a DIFFERENT
    record (d_B = d_A + 1, c_B = s_A - d_B*x) that denotes the SAME point, which is
    exactly what an attacker's collision looks like: two distinct trail records with
    equal x.  The row sequence is identical by construction (same x, same trail id),
    so the cross-trail detector must fire, and the coefficient difference must recover
    the secret.  (The identical-coefficient variant is the TRIVIAL class: same record,
    (da,db) = (0,0); the shipped replay rule rejects it, and the audit reports it
    separately so the two classes are not confused.)
    """
    m = load_selector()
    m.R = 64
    m.MASK = 63
    P = m.P
    secret = 0x1234567890abcdef1234567890abcdef1234567890abcdef1234567890abcdef % m.N
    qX, qY, qZ = m.scalar_mul(secret)
    zi = m.inv(qZ)
    qx, qy = qX * zi * zi % P, qY * zi * zi * zi % P
    tab, _ = m.build_row_table(qx, qy)
    S = 256
    trail = 7
    cA = 0xdeadbeefcafe1234 % m.N
    dA = 0x0badf00d11223344 % m.N
    sA = (cA + dA * secret) % m.N
    dB = (dA + 1) % m.N
    cB = (sA - dB * secret) % m.N
    plants = [("A", cA, dA), ("B", cB, dB), ("A_twin", cA, dA)]
    walks = []
    for name, c, d in plants:
        T = m.scalar_mul((c + d * secret) % m.N)
        rows, coeffs = [], []
        cc, cd = c, d
        for t in range(S):
            X, Y, Z = T
            xa = m.affine_x(X, Y, Z)
            pr = ((xa & m.MASK) ^ ((t * m.TMUL) & m.MASK) ^ ((trail * m.TRAILM) & m.MASK)) & m.MASK
            rows.append((pr, xa))      # (selector row, canonical store key = affine x)
            coeffs.append((cc, cd))
            a, b = m.row_coeffs(pr)
            T = m.jac_madd(X, Y, Z, *tab[pr])
            cc, cd = (cc + a) % m.N, (cd + b) % m.N
        walks.append((name, rows, coeffs))
    same = sum(1 for i in range(S)
               if walks[0][1][i][0] == walks[1][1][i][0] == walks[2][1][i][0])
    store = {}
    cross = trivial = 0
    recoveries = 0
    for widx, (name, rows, coeffs) in enumerate(walks):
        for (pr, xkey), (cc, cd) in zip(rows, coeffs):
            got = store.get(xkey)      # the store key is the point, not the selector row
            if got is None:
                store[xkey] = (cc, cd, widx)
            elif got[2] != widx:
                cross += 1
                pc, pd, _ = got
                da, db = (cc - pc) % m.N, (cd - pd) % m.N
                if da == 0 and db == 0:
                    trivial += 1
                else:
                    cand = (-da * m.inv(db, m.N)) % m.N if db else None
                    if cand == secret:
                        recoveries += 1
    nontrivial = cross - trivial
    assert same == S, "the planted walkers are not row-identical"
    assert cross > 0, "the cross-trail detector did not fire on a planted collision"
    assert recoveries == nontrivial, "every nontrivial planted collision must recover x"
    return {"R": m.R, "steps": S, "records_per_walker": S, "trail_ids": [trail, trail, trail],
            "row_sequence_identical_steps": same,
            "cross_trail_collisions_detected": cross,
            "trivial_class_0_0": trivial,
            "nontrivial_cross_collisions": nontrivial,
            "recoveries": recoveries,
            "yield_measured_nontrivial": recoveries / max(1, nontrivial),
            "real_curve": "P-256 affine x selector, planted c_B = s_A - d_B*x",
            "store_key": "affine x (the 32 B canonical key of the shipped record), NOT the selector row: a row match is an address coincidence, a key match is the success-relevant event",
            "verdict": "PASS: planted same-x records with DIFFERENT coefficients fire the "
                       "cross-trail detector and every nontrivial one recovers the secret; the "
                       "identical-coefficient twin is classified TRIVIAL (0,0), so the two "
                       "classes are separated, not conflated"}


def locator_split(page, words=2):
    return [(page >> (64 * i)) & ((1 << 64) - 1) for i in range(words)]


def locator_join(parts):
    return sum(p << (64 * i) for i, p in enumerate(parts))


def control_c():
    pages = [0, 1, 2 ** 40 + 17, 2 ** 64 + 5, 2 ** 121 + 12345, 2 ** 121 - 1]
    round_trip = []
    for p in pages:
        parts = locator_split(p)
        round_trip.append({"page": str(p), "words": [str(w) for w in parts],
                           "round_trip_equal": locator_join(parts) == p})
    assert all(r["round_trip_equal"] for r in round_trip)

    def u64_mutant(p):
        return p & ((1 << 64) - 1)          # the v21 8-byte entry

    mutant = [{"page": str(p), "mutant_value": str(u64_mutant(p)),
               "corrupts": u64_mutant(p) != p} for p in pages]
    big = [m for m in mutant if m["corrupts"]]
    assert any(int(m["page"]) == 2 ** 64 + 5 for m in big)
    assert big, "the u64 mutant must corrupt at least one page index above 2^64"
    return {"pages_tested": len(pages), "round_trip": round_trip,
            "u64_mutant": mutant,
            "required_id_bits": 122, "shipped_id_bits": 128,
            "verdict": "PASS: 128-bit locators round-trip page 2^64+5 and every tested index; the "
                       "u64 mutant (the v21 entry size) returns the wrong page for %d of %d "
                       "indices, so the v21 8-byte entry cannot address the 2^121.26 pages"
                       % (len(big), len(pages))}


def control_d():
    claim = json.load(open(os.path.join(OUT, "claim_v22.json")))
    ledger = json.load(open(os.path.join(OUT, "ledger_v22.json")))
    exact = ledger["work"]["log2_work_final"]
    c = claim["claim"]
    rat = claim["rational"]
    ulp = 1e-4
    assert c >= exact, "claim must bound the exact work"
    assert c - ulp < exact, "claim is not the tightest 4-decimal bound: one ulp lower still bounds"
    assert abs(rat["numerator"] / rat["denominator"] - c) < 1e-12
    assert int(round(c * 10000)) % (10000 // rat["denominator"]) == 0
    return {"claim": c, "exact_log2_work": exact, "ulp": ulp,
            "claim_minus_ulp": c - ulp, "one_ulp_lower_would_understate": c - ulp < exact,
            "rational": rat,
            "verdict": "PASS: claim is the ceiling of the exact log2 work at 1e-4 granularity and "
                       "the reduced rational is exact"}


def main():
    ps = json.load(open(os.path.join(OUT, "production_selector_v22.json")))
    out = {"a_hazard_match": control_a(ps),
           "b_cross_collision": control_b(),
           "c_locator_128": control_c(),
           "d_regen_ulp": control_d()}
    json.dump(out, open(os.path.join(OUT, "audit_v22.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({k: v["verdict"][:80] for k, v in out.items()}, indent=1))
    print("ALL FOUR CONTROLS FIRED")


if __name__ == "__main__":
    main()
