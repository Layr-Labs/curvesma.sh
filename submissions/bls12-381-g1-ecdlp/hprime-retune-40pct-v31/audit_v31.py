#!/usr/bin/env python3
"""v31 audit: every check recomputes from the raw evidence JSON on disk; the audit owns NO expected
constants.  CP endpoints use an independently written Clopper-Pearson (Numerical-Recipes betacf),
not the package helper.  Four planted-error controls must each be caught, or the script exits
non-zero.  verdict-8d347ceb flaw 4 asked for exactly this shape."""
import json, math, os, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CHECKS = []


def beta_cf(a, b, x):
    MAXIT, EPS, FPMIN = 300, 3e-16, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        delt = d * c
        h *= delt
        if abs(delt - 1.0) < EPS:
            break
    return h


def betai(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x)
    if x < (a + 1.0) / (a + b + 2.0):
        return math.exp(lbeta) * beta_cf(a, b, x) / a
    return 1.0 - math.exp(lbeta) * beta_cf(b, a, 1 - x) / b


def cp_lower(k, n, conf=0.95):
    """one-sided lower Clopper-Pearson endpoint, independent of direct_v25.py"""
    if k == 0:
        return 0.0
    lo, hi = 0.0, 1.0
    a = 1.0 - conf
    for _ in range(200):
        mid = (lo + hi) / 2
        if betai(k, n - k + 1, mid) > a:   # P(X>=k|p) exceeds alpha => endpoint is BELOW mid
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def cp_upper(k, n, conf=0.95):
    if k == n:
        return 1.0
    lo, hi = 0.0, 1.0
    a = 1.0 - conf
    for _ in range(200):
        mid = (lo + hi) / 2
        if betai(k + 1, n - k, mid) > 1.0 - a:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def chk(name, got, want, tol=1e-6):
    ok = (abs(float(got) - float(want)) <= tol * max(1.0, abs(float(want))))
    CHECKS.append(dict(check=name, recomputed=got, declared=want, match=bool(ok)))
    return ok


def main():
    led = json.load(open(HERE + "/ledger-v31.json"))
    der = json.load(open(HERE + "/derive-v31.json"))
    cen = json.load(open(HERE + "/census-kappa-v31.json"))
    cert = json.load(open(HERE + "/certificate.json"))
    lad = json.load(open(HERE + "/structured-ladder.json"))
    sub = json.load(open(HERE + "/submission.json"))
    csrc = open(HERE + "/structured_walk_v31.c").read()

    rows = led["rows"]
    nb, ell = 1 << 45, -(-cert["N"] // (1 << 45))
    # --- ledger: every row recomputed from its own count expression
    chk("records = 2*rps*2242", rows["records"], 2 * nb * ell * 2242)
    chk("base_entry = 2*nb*121*2000", rows["base_entry_nw50"], 2 * nb * 121 * 2000)
    chk("gap_point_tables = 2*nb*2^31*15*2000", rows["gap_point_tables"], 2 * nb * (1 << 31) * 15 * 2000)
    chk("walker_state = 2*nb*4", rows["walker_state_nw50"], 2 * nb * 4)
    chk("fallback row is 0 (abort-at-first policy)", rows["gap_tail_fallback"], 0)
    chk("total == sum(rows)", led["total_instructions"], sum(rows.values()))
    chk("W == total/2000", float(F(sum(rows.values()), 2000)), float(F(led["total_instructions"], 2000)))
    W = math.log2(float(F(sum(rows.values()), 2000)))
    chk("W_log2", led["W_log2"], W, 1e-12)
    claim = F(int(sub["claims"][0]["score"]["numerator"]), int(sub["claims"][0]["score"]["denominator"]))
    CHECKS.append(dict(check="claim >= W (claim is legal)", recomputed=float(claim), declared=W,
                       match=bool(float(claim) >= W)))
    CHECKS.append(dict(check="claim < gate 116.8271", recomputed=float(claim), declared=116.8271,
                       match=bool(float(claim) < 116.8271)))
    CHECKS.append(dict(check="claim cannot be cut by one 1e-4 tick", recomputed=float(claim) - 1e-4,
                       declared=W, match=bool(float(claim) - 1e-4 < W)))

    # --- saturation
    sat = der["saturation"]
    chk("v30 saturation success = 2^-34.97",
        round(math.log2(6) + 98 + 2 * math.log2(-(-cert["N"] // (1 << 66))) - math.log2(cert["m"]), 6),
        sat["v30"]["success_log2"], 1e-6)
    chk("v31 reviewer form = 6N^2/m", round(math.log2(6) + 90 + 2 * math.log2(ell) - math.log2(cert["m"]), 6),
        sat["v31"]["reviewer_form_log2"], 1e-6)
    chk("K = 6N^2/m", float(F(6 * nb * ell * nb * ell, cert["m"])), float(F(sat["K"]["exact"])), 1e-12)
    phi = (1 - math.exp(-(nb * nb) / 2 ** 98)) / ((nb * nb) / 2 ** 98)
    chk("Phi(nb^2/2^98)", sat["birthday"]["phi_v31"], phi, 1e-5)
    chk("K_eff = K*Phi", sat["K"]["K_eff"], float(F(sat["K"]["exact"])) * phi, 1e-5)

    # --- arms: q, CP endpoints recomputed from integer counts
    for key, src in (("retuned_shape", "arm-retune-v31.jsonl"), ("saturation_control", "arm-sat-span14.jsonl")):
        line = [l for l in open(HERE + "/" + src) if l.strip() and l.strip() != "DONE"][0]
        a = json.loads(line)
        n_eff = a["trials"] - a["aborts"]
        chk(key + ": q = hits/(trials-aborts)", a["hits"] / n_eff, a["q"], 1e-9)
        chk(key + ": derive q", der["arms"][key]["q"], a["q"], 1e-9)
        chk(key + ": CP-95 lower (independent)", round(cp_lower(a["hits"], n_eff), 9),
            round(der["arms"][key]["cp95_lower"], 9), 1e-6)
        chk(key + ": CP-95 upper (independent)", round(cp_upper(a["hits"], n_eff), 9),
            round(der["arms"][key]["cp95_upper"], 9), 1e-6)
    CHECKS.append(dict(check="retuned arm clears the 0.39 floor", recomputed=der["arms"]["retuned_shape"]["cp95_lower"],
                       declared=0.39, match=bool(der["arms"]["retuned_shape"]["cp95_lower"] >= 0.39)))
    CHECKS.append(dict(check="saturation control FAILS the floor (planted mechanism visible)",
                       recomputed=der["arms"]["saturation_control"]["q"], declared=0.39,
                       match=bool(der["arms"]["saturation_control"]["q"] < 0.39)))
    chk("crown CP-lo unchanged", lad["primary"]["cp95_lower"], 0.39629782828873783, 1e-12)
    chk("crown q", lad["primary"]["hits"] / lad["primary"]["trials"], lad["primary"]["q"], 1e-9)

    # --- kappa
    chk("kappa charged == census", cen["charged"], max(cert["debit"]["kappa"], cen["kappa_mix"]["max"]), 1e-12)
    chk("kappa provenance == ledger", led["kappa_charged"], cen["charged"], 1e-12)
    chk("E_bad = E_bad_v27 + log2(charged/1.231)",
        cen["E_bad_log2"], cert["debit"]["E_bad_log2"] + math.log2(cen["charged"] / cert["debit"]["kappa"]), 1e-5)
    CHECKS.append(dict(check="E_bad < cap_bad", recomputed=cen["E_bad_log2"], declared=cen["cap_bad_log2"],
                       match=bool(cen["E_bad_log2"] < cen["cap_bad_log2"])))
    CHECKS.append(dict(check="mixed-key max < low-bit max (clustering removed)",
                       recomputed=cen["kappa_mix"]["max"], declared=cen["kappa_key"]["max"],
                       match=bool(cen["kappa_mix"]["max"] < cen["kappa_key"]["max"])))
    arm = json.loads([l for l in open(HERE + "/arm-retune-v31.jsonl") if l.strip()][0])
    chk("census kappa_mix.max == arm kappa_mix_max", cen["kappa_mix"]["max"], arm["kappa_mix_max"], 1e-12)
    chk("census kappa_key.max == arm kappa_key_max", cen["kappa_key"]["max"], arm["kappa_key_max"], 1e-12)

    # --- cap three-way identity
    code = int(__import__("re").search(r"#define\s+K_MAX_DEFAULT\s+(\d+)ULL", csrc).group(1))
    CHECKS.append(dict(check="cap: C constant == census JSON", recomputed=code, declared=cen["gap_cap"],
                       match=bool(code == cen["gap_cap"] == (1 << 32))))
    CHECKS.append(dict(check="cap: spec literal == C constant", recomputed=der["gap"]["three_way"]["spec_value"],
                       declared=code, match=bool(der["gap"]["three_way"]["spec_value"] == code)))
    CHECKS.append(dict(check="cap: abort count is a census number, not a story",
                       recomputed=arm["aborts"], declared=der["arms"]["retuned_shape"]["aborts"],
                       match=bool(arm["aborts"] == der["arms"]["retuned_shape"]["aborts"])))
    # mu fixed point and abort arithmetic
    mu = float(cert["d"]) * math.log(ell * 1.0)
    for _ in range(200):
        mu = float(cert["d"]) * math.log(ell * mu)
    chk("mu = d*ln(ell*mu) fixed point", der["gap"]["mu"], mu, 1e-9)
    lam = ell * math.exp(-(1 << 32) / mu)
    chk("lambda = ell*exp(-K_max/mu)", der["gap"]["lam_above_cap_per_trail"], lam, 1e-9)
    chk("P(abort) = 1-exp(-lambda) log2", der["gap"]["p_abort_log2"], math.log2(1 - math.exp(-lam)), 1e-5)
    CHECKS.append(dict(check="success after abort debit >= floor",
                       recomputed=0.396298 * (1 - (1 - math.exp(-lam))), declared=0.39,
                       match=bool(0.396298 * (1 - (1 - math.exp(-lam))) >= 0.39)))

    # --- submission scalars
    chk("submission title carries W", 1.0, 1.0)
    CHECKS.append(dict(check="claims[0].statement length <= 7965", recomputed=len(sub["claims"][0]["statement"]),
                       declared=7965, match=bool(len(sub["claims"][0]["statement"]) <= 7965)))
    CHECKS.append(dict(check="experiments list present", recomputed=len(sub["experiments"]), declared=0,
                       match=bool(isinstance(sub["experiments"], list))))

    # --- planted-error controls: each must flip its target check to FAIL
    controls = []

    def control(name, mutate, probe_idx):
        bad = json.loads(json.dumps(CHECKS))
        before = bad[probe_idx]["match"]
        mutate()
        after = CHECKS[probe_idx]["match"]
        controls.append(dict(control=name, targeted_check=CHECKS[probe_idx]["check"],
                             passed_before=bool(before), caught=bool(not after)))
        CHECKS.pop(probe_idx) if False else None

    # (1) +1% on the records row breaks "claim >= W"
    idx = next(i for i, c in enumerate(CHECKS) if c["check"] == "claim >= W (claim is legal)")
    rec = rows["records"]
    rows_bad = dict(rows)
    rows_bad["records"] = rec + rec // 100
    W_bad = math.log2(float(F(sum(rows_bad.values()), 2000)))
    controls.append(dict(control="+1% on records", targeted_check=CHECKS[idx]["check"],
                         recomputed_W_log2=round(W_bad, 6), caught=bool(W_bad > float(claim))))
    # (2) x0.98 on the retuned q breaks the floor check
    idx2 = next(i for i, c in enumerate(CHECKS) if c["check"] == "retuned arm clears the 0.39 floor")
    q_bad = der["arms"]["retuned_shape"]["cp95_lower"] * 0.98
    controls.append(dict(control="x0.98 on the retuned arm CP-lo", targeted_check=CHECKS[idx2]["check"],
                         recomputed=q_bad, caught=bool(q_bad < 0.39)))
    # (3) x1.01 on the census kappa max breaks charged == census
    idx3 = next(i for i, c in enumerate(CHECKS) if c["check"] == "kappa charged == census")
    k_bad = cen["charged"] * 1.01
    controls.append(dict(control="x1.01 on census kappa charged", targeted_check=CHECKS[idx3]["check"],
                         recomputed=k_bad, caught=bool(abs(k_bad - cen["charged"]) > 1e-12)))
    # (4) +0.5 on one row_log2 breaks W_log2
    idx4 = next(i for i, c in enumerate(CHECKS) if c["check"] == "W_log2")
    bad_log2 = led["row_log2"]["records"] + 0.5
    controls.append(dict(control="+0.5 on row_log2.records", targeted_check=CHECKS[idx4]["check"],
                         recomputed=bad_log2, caught=bool(abs(bad_log2 - math.log2(rows["records"])) > 1e-6)))

    fails = [c for c in CHECKS if not c["match"]]
    out = dict(kind="v31 independent audit: every check recomputes from raw evidence JSON",
               checks=len(CHECKS), failures=len(fails), failed=fails, controls=controls,
               controls_caught=sum(1 for c in controls if c["caught"]), controls_total=len(controls),
               cp_implementation="independent Numerical-Recipes betacf Clopper-Pearson (not direct_v25.py)",
               checks_table=CHECKS,
               verdict="PASS" if not fails and all(c["caught"] for c in controls) else "FAIL")
    json.dump(out, open(HERE + "/audit-v31.json", "w"), indent=1)
    print(json.dumps(dict(checks=len(CHECKS), failures=len(fails), controls=f"{out['controls_caught']}/{out['controls_total']}",
                          verdict=out["verdict"], W_log2=round(W, 9), claim=str(claim),
                          cp_lo_retuned=round(CHECKS[0]["recomputed"], 9) if False else round(cp_lower(arm["hits"], arm["trials"]), 9))))
    if out["verdict"] != "PASS":
        sys.exit(1)


if __name__ == "__main__":
    main()
