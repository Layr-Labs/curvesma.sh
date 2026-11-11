"""Audit of the capped-descent (6,2) special TNFS package (Python standard library only).

Input (JSON on stdin): {"seed": int, "fresh_cells": int}.  Output: one JSON object.
[E1] tower and field facts       [E2] Dickman spot values          [E3] small-ideal hit sums
[E4] main sample replay + fresh  [E5] special-q replay + scaling   [E6] boot replay
[E7] ECM re-measurement          [E8] ledger = declared = text = manifest
[E9] toy special-q ball counts (fresh sweeps + audit of the stored campaign)
Nothing here runs the production attack; see proof.md section 11 for the scope of each check.
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
for cand in (HERE, os.path.join(HERE, "files"), "/snapshot/files"):
    if os.path.exists(os.path.join(cand, "model.py")):
        BASE = cand
        break
sys.path.insert(0, BASE)

import ecm      # noqa: E402
import lat      # noqa: E402
import model    # noqa: E402

T0 = time.time()
U = 0x44E992B44A6909F1
P = 36 * U**4 + 36 * U**3 + 24 * U**2 + 6 * U + 1
R = 36 * U**4 + 36 * U**3 + 18 * U**2 + 6 * U + 1
P_PROFILE = 21888242871839275222246405745257275088696311157297823662689037894645226208583
R_PROFILE = 21888242871839275222246405745257275088548364400416034343698204186575808495617


def e1_tower():
    out = {}
    out["p_matches_profile"] = P == P_PROFILE and lat.P == P
    out["r_matches_profile"] = R == R_PROFILE
    out["r_divides_p4_p2_1"] = (P**4 - P**2 + 1) % R == 0
    out["gcd_r_p6_1"] = math.gcd(R, P**6 - 1) == 1
    out["u_factorisation"] = 3 * 13 * 29 * 983 * 11003 * 405928799 == U
    # h irreducible mod p (Rabin): y^(p^6) = y, and gcd(y^(p^k) - y, h) = 1 for k = 2, 3
    y = [0, 1, 0, 0, 0, 0]
    frob = [y]
    for _ in range(6):
        frob.append(lat.hpow_mod(frob[-1], P) if hasattr(lat, "hpow_mod") else _hpow(frob[-1], P))
    out["h_frobenius_order_6"] = frob[6] == y and frob[1] != y and frob[2] != y and frob[3] != y
    ok = True
    for k in (1, 2, 3):
        d = [(a - b) % P for a, b in zip(frob[k], y)]
        ok &= _gcd_with_h_is_one(d)
    out["h_irreducible_mod_p"] = bool(ok and out["h_frobenius_order_6"])
    # u*y is a non-square in F_(p^6)
    uy = [0, U % P, 0, 0, 0, 0]
    e = _hpow(uy, (P**6 - 1) // 2)
    out["uy_nonsquare"] = e == [P - 1, 0, 0, 0, 0, 0]
    # common root: with x^2 = u y, f = y^4 P(u) = y^4 p
    x2 = [0, U, 0, 0, 0, 0]
    x4 = lat.hmul(x2, x2)
    x6 = lat.hmul(x4, x2)
    x8 = lat.hmul(x4, x4)
    Y1, Y2, Y3, Y4 = ([0] * i + [1] + [0] * (5 - i) for i in (1, 2, 3, 4))
    f = [36 * a + 36 * b + 24 * c + 6 * d + e_ for a, b, c, d, e_ in zip(
        x8, lat.hmul(Y1, x6), lat.hmul(Y2, x4), lat.hmul(Y3, x2), Y4)]
    out["f_at_root_of_g_is_y4_p"] = f == [P * v for v in Y4]
    out["units_norm"] = [abs(lat.hnorm(x)) for x in (lat.Y, lat.E2, lat.E3)] == [1, 1, 1]
    return out


def _hpow(a, e):
    r = [1, 0, 0, 0, 0, 0]
    a = [x % P for x in a]
    while e:
        if e & 1:
            r = lat.hmul(r, a, P)
        a = lat.hmul(a, a, P)
        e >>= 1
    return r


def _gcd_with_h_is_one(d):
    """gcd(d(y), h(y)) over F_p is constant."""
    a = [x % P for x in lat.H]
    b = [x % P for x in d]

    def trim(v):
        while v and v[-1] == 0:
            v.pop()
        return v
    a, b = trim(a), trim(b)
    while b:
        inv = pow(b[-1], P - 2, P)
        while len(a) >= len(b):
            c = a[-1] * inv % P
            sh = len(a) - len(b)
            for i, bv in enumerate(b):
                a[i + sh] = (a[i + sh] - c * bv) % P
            trim(a)
        a, b = b, a
    return len(a) == 1


def e2_dickman():
    D = model._D()
    return {"rho_2": D(2.0), "rho_5": D(5.0), "rho_8": D(8.0), "rho_10": D(10.0),
            "ok": abs(D(2.0) - 0.3068528) < 1e-4 and abs(D(5.0) / 3.5472e-4 - 1) < 1e-3 and abs(D(10.0) / 2.7701e-11 - 1) < 1e-3}


def e3_small_sums():
    sums = []
    for side in "fg":
        tot = 0.0
        for q in range(5, 4096):
            if q != 197 and all(q % k for k in range(2, int(q**0.5) + 1)):
                tot += len(lat.degree_one_ideals(q, side)) / q
        sums.append(tot)
    return {"sum_f": sums[0], "sum_g": sums[1],
            "ok": all(abs(a - b) < 1e-9 for a, b in zip(sums, model.SMALL_SUMS))}


def _unsieved_bits(n):
    """Upper bound on the mass the sieve does not credit: sum over primes q < 4096 with q^v || n of (v-1) log2 q."""
    x = 0.0
    for q in lat.PRIMES_1000 if False else _P4096:
        if n % q == 0:
            v = 0
            while n % q == 0:
                n //= q
                v += 1
            x += (v - 1) * math.log2(q)
    return x


_P4096 = [q for q in range(2, 4096) if all(q % k for k in range(2, int(q**0.5) + 1))]


def e4_main(data, seed, fresh):
    smp = data["main"]
    A = smp["A"]
    rows = {c[0]: c[1:] for c in lat.unpack_cells(smp["adm"], "33333")}
    rng = lat.Rng(smp["seed"])
    ok, checked, adm_seen = True, 0, 0
    for i in range(60):
        a, b = lat.uniform_cell(rng, A)
        nf, ng = lat.norm_f(a, b), lat.norm_g(a, b)
        adm = lat.admissible(a, b, nf, ng)
        ok &= adm == (i in rows)
        if adm:
            adm_seen += 1
            ok &= rows[i] == [lat.t10(nf), lat.t10(ng), lat.t10(lat.smooth_part(nf)), lat.t10(lat.smooth_part(ng))]
        checked += 1
    out = {"stored_cells_replayed": checked, "stored_admissible_among_them": adm_seen, "replay_ok": bool(ok)}
    rng = lat.Rng("fresh|%d" % seed)
    table = lat.unit_table(4)
    nfs, ngs, n_adm, kept, inv_m, xmax = [], [], 0, 0, [], 0.0
    for _ in range(fresh):
        a, b = lat.uniform_cell(rng, A)
        nf, ng = lat.norm_f(a, b), lat.norm_g(a, b)
        if not lat.admissible(a, b, nf, ng):
            continue
        n_adm += 1
        nfs.append(lat.log2i(nf))
        ngs.append(lat.log2i(ng))
        xf, xg = _unsieved_bits(nf), _unsieved_bits(ng)
        xmax = max(xmax, xf, xg)
        kept += xf <= 47 and xg <= 47
        if len(inv_m) < 80:
            inv_m.append(1.0 / lat.unit_multiplicity(a, b, A, table))
    out.update(fresh_cells=fresh, fresh_admissible=n_adm, fresh_admissible_fraction=n_adm / fresh,
               fresh_mean_log2_nf=sum(nfs) / max(1, n_adm), fresh_mean_log2_ng=sum(ngs) / max(1, n_adm),
               fresh_retention_slack47=kept / max(1, n_adm), fresh_max_unsieved_bits=xmax,
               fresh_mean_inverse_multiplicity=sum(inv_m) / max(1, len(inv_m)), multiplicity_cells=len(inv_m))
    out["fresh_ok"] = bool(n_adm > 0 and out["fresh_retention_slack47"] >= model.CFG["retention"]
                           and (len(inv_m) < 40 or out["fresh_mean_inverse_multiplicity"] >= model.CFG["inv_mult_flat"]))
    return out


def e5_desc(data, seed):
    bins = data["desc"]["bins"]
    rng = lat.Rng("pick|%d" % seed)
    keys = sorted(bins)
    picks = []
    for _ in range(4):
        k = keys[rng.below(len(keys))]
        picks.append((k, rng.below(16)))
    picks.append(("f104", 0))
    ok, cells_checked, scal = True, 0, {}
    for k, i in picks:
        side, n = k[0], int(k[1:])
        sd = "desc|%s|%d|%d" % (side, n, i)
        q, s, t, _ = lat.draw_ideal(lat.Rng(sd), side, n)
        basis = lat.specialq_reduced(q, s, t)
        ok &= abs(lat.det(basis)) == q and all(lat.specialq_member(q, s, t, r) for r in basis)
        stored = {c[0]: c[1:] for c in lat.unpack_cells(bins[k][i], "233")}
        means = {}
        for Rr in ((16, 4, 64) if (k, i) == ("f104", 0) else (16,)):
            cr = lat.Rng(sd + "|R=%d" % Rr)
            lf, lg = [], []
            for j in range(48):
                v = lat.combine(lat.draw_coeffs(cr, Rr), basis)
                a, b = v[:6], v[6:]
                nf, ng = lat.norm_f(a, b), lat.norm_g(a, b)
                ok &= (nf if side == "f" else ng) % q == 0
                lf.append(lat.log2i(nf))
                lg.append(lat.log2i(ng))
                if Rr == 16:
                    adm = lat.admissible(a, b, nf, ng)
                    ok &= adm == (j in stored)
                    if adm:
                        ok &= stored[j] == [lat.t10(nf), lat.t10(ng)]
                    cells_checked += 1
            means[Rr] = (sum(lf) / 48, sum(lg) / 48)
        if len(means) == 3:
            scal = {"R4_shift_f": means[4][0] - means[16][0], "R4_shift_g": means[4][1] - means[16][1],
                    "R64_shift_f": means[64][0] - means[16][0], "R64_shift_g": means[64][1] - means[16][1]}
    # the model's scaling must be on the conservative side of the stored shifts
    cons = True
    for v in data["desc"]["scale_shift_bits"].values():
        # model: R = 64 -> (+96.5, +25); R = 256 -> (+193, +50)   (plans use radii of at least 16)
        cons &= v["R64"][0] <= 96.5 and v["R64"][1] <= 25.0 and v["R256"][0] <= 193.0 and v["R256"][1] <= 50.0
    return {"ideals_regenerated": [list(p_) for p_ in picks], "cells_compared": cells_checked, "replay_ok": bool(ok),
            "scaling_one_ideal": scal, "stored_scaling_shifts": data["desc"]["scale_shift_bits"],
            "model_scaling_conservative": bool(cons)}


def _inv_rows(B):
    from fractions import Fraction
    n = len(B)
    M = [[Fraction(x) for x in r] + [Fraction(int(i == j)) for j in range(n)] for i, r in enumerate(B)]
    for c in range(n):
        p_ = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p_] = M[p_], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [r[n:] for r in M]


def e5c_lambda(seed, led):
    """Per-ideal lambda study (lambda_spread.json): regenerate two of its ideals and audit the stored summary."""
    ls = json.load(open(os.path.join(BASE, "lambda_spread.json")))
    bins = ls["bins"]
    plans = {(r["side"], r["top_bits"]): r for r in led["nodes"]}
    D = model._D()
    af, ag = led["main"]["alpha_f"], led["main"]["alpha_g"]
    rng = lat.Rng("lampick|%d" % seed)
    keys = sorted(bins)
    ok, recs = True, []
    for _ in range(2):
        k = keys[rng.below(len(keys))]
        i = rng.below(ls["ideals_per_bin"])
        side, n = k[0], int(k[1:])
        pl = plans[(side, n)]
        c, Rr = pl["child_bound_bits"], pl["radius"]
        ok &= bins[k]["child_bound"] == c and bins[k]["radius"] == Rr
        df, dg = model.shift_bits(Rr)
        sd = "lam|%s|%d|%d" % (side, n, i)
        q, s, t, _ = lat.draw_ideal(lat.Rng(sd), side, n)
        B = lat.specialq_reduced(q, s, t)
        cr = lat.Rng(sd + "|R=16")
        acc = 0.0
        for _ in range(ls["cells_per_ideal"]):
            v = lat.combine(lat.draw_coeffs(cr, 16), B)
            a, b = v[:6], v[6:]
            nf, ng = lat.norm_f(a, b), lat.norm_g(a, b)
            if not lat.admissible(a, b, nf, ng):
                continue
            lf, lg = lat.log2i(nf) + df, lat.log2i(ng) + dg
            if side == "f":
                lf -= n - 1
            else:
                lg -= n - 1
            acc += D((lf * model.LN2 + af) / (c * model.LN2)) * D((lg * model.LN2 + ag) / (c * model.LN2))
        lam = (2 * Rr + 1) ** 12 / 2.0 * acc / ls["cells_per_ideal"] * (1 - pl["missing"]) * model.CFG["retention"]
        l2 = math.log2(lam) if lam > 0 else -99.0
        ok &= l2 >= bins[k]["min"] - 0.01
        recs.append({"bin": k, "ideal": i, "log2_lambda": round(l2, 2), "bin_min": bins[k]["min"], "bin_median": bins[k]["median"]})
    total = sum(v["ideals"] for v in bins.values())
    lam_min = int(2 ** min(v["min"] for v in bins.values()))
    p_empty = math.exp(-lam_min / (model.CFG["mean_allowance"] * model.CFG["cluster_allowance"]))
    budget = p_empty + 3.0 / total
    unit_cells = sum(v["unit_cells"] for v in bins.values())
    unit_mean = max(v["unit_mult_mean"] for v in bins.values())
    return {"regenerated": recs, "ideals_total": total, "smallest_lambda": lam_min,
            "estimates_below_256": sum(v["below_256"] for v in bins.values()),
            "p_empty_under_H_indep": p_empty, "zero_count_95_upper": 3.0 / total, "sum": budget,
            "eps_charged": model.CFG["eps_empty"], "unit_cells": unit_cells, "largest_bin_mean_unit_multiplicity": unit_mean,
            "largest_unit_multiplicity": max(v["unit_mult_max"] for v in bins.values()),
            "cells_where_exponent_bound_4_and_6_differ": sum(v["unit_bound4_differs"] for v in bins.values()),
            "ok": bool(ok and budget <= model.CFG["eps_empty"] and 2 * unit_mean <= model.CFG["cluster_allowance"])}


def e6_boot(data, seed):
    """One stored boot sample: T uniform in GT from its seed, slope tau = t0/t1, weighted LLL, first vector."""
    cells = lat.unpack_cells(data["boot"]["cells"], "313")
    i = lat.Rng("bootpick|%d" % seed).below(len(cells))
    T = lat.gt_element(lat.Rng("bootgt|%d" % i))
    one = ([1, 0, 0, 0, 0, 0], [0] * 6)
    in_gt = lat.f12_pow(T, R) == one and T != one
    tau = lat.boot_tau(T)
    basis = lat.boot_reduced(tau, data["boot"]["wbits"])
    v = basis[0]
    a, b = v[:6], v[6:]
    # a - b x is an F_(p^6)-multiple of T = t0 + t1 x  iff  a t1 + b t0 = 0
    lin = [(x + y) % P for x, y in zip(lat.hmul([z % P for z in a], T[1], P), lat.hmul([z % P for z in b], T[0], P))]
    ng, nf = lat.norm_g(a, b), lat.norm_f(a, b)
    rec = [lat.t10(ng), 1 if lat.admissible(a, b, nf, ng) else 0, lat.t10(lat.smooth_part(ng))]
    vals = [c[0] / 10 for c in cells]
    return {"boot_index": i, "target_in_GT": bool(in_gt), "member": bool(lat.boot_member(tau, v)) and not any(lin),
            "det_is_p6": abs(lat.det(basis)) == P**6, "log2_ng": lat.log2i(ng), "replay_ok": rec == cells[i],
            "stored_mean_log2_ng": sum(vals) / len(vals), "stored_max_log2_ng": max(vals),
            "stored_admissible_fraction": sum(c[1] for c in cells) / len(cells)}


def e7_ecm(seed):
    tab = json.load(open(os.path.join(BASE, "ecm_table.json")))
    worst, cells_with_captures, below = None, 0, 0
    for bits, u, curves, captures, rho_u in tab["capture"]:
        if captures == 0:
            continue
        cells_with_captures += 1
        ratio = captures / curves / rho_u
        below += ratio < 1.0
        if captures >= 60 and (worst is None or ratio < worst[0]):
            worst = (ratio, bits, u)
    rp = ecm.replay(seed, 1200)
    s = rp["capture_32bit_u4"]
    bits_ok = all(ecm.stage1_ops(b1)["xdbl"] <= 1.8107 * b1 for b1 in (256, 4096))
    return {"stored_cells_with_captures": cells_with_captures, "stored_cells_below_rho": below,
            "stored_worst_ratio_with_60_captures": worst, "stored_planted": tab["planted"],
            "replay_32bit_u4": {k: s[k] for k in ("curves", "captures", "p", "rho_u", "ratio") if k in s},
            "replay_planted_600bit_all_isolated": rp["planted_all_isolated"], "replay_planted_600bit_n": len(rp["planted_600"]),
            "ladder_steps_within_bound": bool(bits_ok),
            "ok": below == 0 and bits_ok and rp["planted_all_isolated"] == len(rp["planted_600"])}


def _close(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(_close(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)):
        return isinstance(b, (list, tuple)) and len(a) == len(b) and all(_close(x, y) for x, y in zip(a, b))
    if isinstance(a, float) or isinstance(b, float):
        return abs(a - b) <= 1e-7 * max(1.0, abs(a), abs(b))
    return a == b


def e8_ledger(led):
    decl = json.load(open(os.path.join(BASE, "declared.json")))
    proof = open(os.path.join(BASE, "proof.md")).read()
    man = json.load(open(os.path.join(BASE, "submission.json")))
    total = led["log2_expected_total_units"]
    num = math.ceil(total * 100)
    sc = man["claims"][0]["score"]
    from fractions import Fraction
    sc_val = Fraction(int(sc["numerator"]), int(sc["denominator"]))     # the host may store a reduced fraction
    score_ok = sc_val == Fraction(num, 100) == Fraction(int(decl["score"]["numerator"]), int(decl["score"]["denominator"]))
    missing = [k for k, v in decl["strings"].items() if v not in proof]
    statement = man["claims"][0]["statement"]
    st_ok = all(decl["strings"][k] in statement for k in ("SCORE", "TOTAL", "RUN_UNITS", "POST", "GATE", "METER"))
    return {"ledger_equals_declared": _close(led, decl["ledger"]),
            "node_table_in_proof": decl["node_table"] in proof,
            "declared_strings": len(decl["strings"]), "declared_strings_missing_from_proof": missing,
            "manifest_score": sc, "manifest_score_matches": bool(score_ok),
            "manifest_statement_numbers_match": bool(st_ok),
            "success": led["success"], "success_at_least_0_39": led["success"] >= 0.39,
            "log2_expected_total_units": total, "score_bounds_total": num / 100 >= total,
            "ledger_log2_instructions": led["ledger_log2_instructions"], "log2_run_units": led["log2_run_units"],
            "main": led["main"], "boot": led["boot"],
            "min_log2_lambda_worst": min(r["log2_lambda_worst"] for r in led["nodes"]),
            "min_log2_lambda_shifted": min(r["log2_lambda_shifted"] for r in led["nodes"]),
            "ok": _close(led, decl["ledger"]) and not missing and decl["node_table"] in proof and st_ok
            and score_ok and led["success"] >= 0.39}


def e9_toy(seed):
    """Fresh complete sweeps of toy special-q balls (true node predicate) and an audit of the stored campaign:
    in every row (population, side, radius) the observed empty fraction is compared with exp(-mean/D), D the
    measured variance/mean of the counts (the independent-cluster law), up to three standard errors."""
    import ballcount
    r = ballcount.run(seed, 3)
    camp = r.get("campaign", {})
    rows, ok, total = [], True, 0
    for p in camp.get("by_radius", []):
        bound = math.exp(-p["mean"] / p["D"])
        expect = p["balls"] * bound
        within = p["empty"] <= expect + 3 * math.sqrt(expect) + 3
        rows.append({"pop": p["pop"], "side": p["side"], "R": p["R"], "balls": p["balls"], "empty": p["empty"],
                     "mean": p["mean"], "D": p["D"], "empty_fraction": p["empty"] / p["balls"],
                     "exp_minus_mean_over_D": bound, "within_law": bool(within),
                     "realized_over_dickman": p["real_over_dickman"]})
        ok &= within and p["real_over_dickman"] >= 1 / model.CFG["mean_allowance"]
        total += p["balls"]
    fresh = [{k: b[k] for k in ("side", "q", "realized", "predicted_dickman", "predicted_table") if k in b}
             for b in r.get("balls", [])]
    return {"fresh_balls": fresh, "campaign_rows": rows, "total_balls": total,
            "ok": bool(ok and rows and fresh)}


def main():
    try:
        inp = json.loads(sys.stdin.read() or "{}")
    except Exception:
        inp = {}
    seed = int(inp.get("seed", 20261111))
    fresh = int(inp.get("fresh_cells", 400))
    data = model.load(BASE)
    out = {"input": {"seed": seed, "fresh_cells": fresh}}
    out["E1_tower"] = e1_tower()
    out["E2_dickman"] = e2_dickman()
    out["E3_small_sums"] = e3_small_sums()
    out["E4_main"] = e4_main(data, seed, fresh)
    out["E5_special_q"] = e5_desc(data, seed)
    led = model.ledger(data)
    out["E5c_lambda_spread"] = e5c_lambda(seed, led)
    out["E6_boot"] = e6_boot(data, seed)
    out["E7_ecm"] = e7_ecm(seed)
    out["E8_ledger"] = e8_ledger(led)
    out["E9_toy_ball_counts"] = e9_toy(seed)
    checks = {
        "E1": all(v is True for v in out["E1_tower"].values()),
        "E2": out["E2_dickman"]["ok"], "E3": out["E3_small_sums"]["ok"],
        "E4": out["E4_main"]["replay_ok"] and out["E4_main"]["fresh_ok"],
        "E5": out["E5_special_q"]["replay_ok"] and out["E5_special_q"]["model_scaling_conservative"]
        and out["E5c_lambda_spread"]["ok"],
        "E6": out["E6_boot"]["replay_ok"] and out["E6_boot"]["member"] and out["E6_boot"]["det_is_p6"] and out["E6_boot"]["target_in_GT"],
        "E7": out["E7_ecm"]["ok"], "E8": out["E8_ledger"]["ok"], "E9": out["E9_toy_ball_counts"]["ok"],
    }
    out["checks"] = checks
    out["all_pass"] = all(checks.values())
    out["seconds"] = round(time.time() - T0, 2)
    json.dump(out, sys.stdout)


if __name__ == "__main__":
    main()
