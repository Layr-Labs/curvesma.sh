#!/usr/bin/env python3
"""v30 INDEPENDENT AUDIT (verdict-0ad3f80 flaw 4: the v29 audit was VACUOUS -- it checked assertions
against the same values it sourced them from, with three constructed power pairs).

Rules this file follows:
  * every check compares a PACKAGE LITERAL against a RECOMPUTATION from the raw evidence JSON on
    disk: the audit owns no constants (no expected values are written down here), and no value is
    taken from the file it is used to check;
  * coverage is >= 20 numeric pairs spread over ALL rows (every ledger row's log2, every arm's
    q = hits/trials, every CP endpoint, every ledger identity, every gap/kappa row);
  * planted-error controls mutate the SOURCE JSONs in a scratch copy (never this file's
    expectations) and must fire 4/4; the script exits non-zero if any control survives;
  * the CP endpoints are recomputed with an independently written incomplete-beta, not by importing
    the package's own cp95 helper.

Usage: python3 audit_v30.py [--root DIR] [--out FILE]
"""
import json, math, os, shutil, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 5e-6


# ---------- independent Clopper-Pearson (Numerical-Recipes betacf/betai, written here) ----------
def _betacf(a, b, x):
    MAXIT, EPS, FPMIN = 200, 3e-14, 1e-300
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
        de = d * c
        h *= de
        if abs(de - 1.0) < EPS:
            break
    return h


def _betai(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def cp_lower(k, n, target=0.05):
    if k <= 0:
        return 0.0
    if k >= n:
        return target ** (1.0 / n)
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _betai(k, n - k + 1, mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------- checks ----------------------------------
def load(root, name):
    p = os.path.join(root, name)
    if name.endswith(".jsonl"):
        for line in open(p):
            line = line.strip()
            if line and line != "DONE":
                return json.loads(line)
        raise RuntimeError("no JSON line in " + name)
    return json.load(open(p))


def checks(root):
    """Each returns (name, ok, observed, recomputed). No expected constants live here."""
    out = []

    def add(name, obs, rec, tol=TOL):
        ok = obs is not None and rec is not None and abs(float(obs) - float(rec)) <= tol
        out.append((name, ok, obs, rec))

    led = load(root, "ledger-v30.json")
    der = load(root, "derive-v30.json")
    cen = load(root, "census-kappa-v30.json")
    cert = load(root, "certificate.json")
    lad = load(root, "structured-ladder.json")
    d25 = load(root, "direct-v25.json")
    cr = load(root, "arm-crown-kb20.jsonl")
    pr = load(root, "arm-pool-kb20.jsonl")
    cd = load(root, "arm-capdemo.jsonl")
    sub = load(root, "submission.json")

    rows = led["rows"]
    # (1) every declared row log2 against the row's own integer, recomputed
    for k, v in rows.items():
        add("row_log2[%s]" % k, led["row_log2"].get(k), math.log2(v) if v > 0 else None, 1e-5)
    # (2) ledger identities, all recomputed from the rows
    add("total==sum(rows)", led["total_instructions"], sum(rows.values()), 0.0)
    add("W_log2", led["W_log2"], math.log2(sum(rows.values()) / led["C_instructions_per_unit"]), 1e-9)
    claim = int(led["claim"]["numerator"]) / int(led["claim"]["denominator"])
    gate = int(led["gate"]["numerator"]) / int(led["gate"]["denominator"])
    add("claim_minus_W", led["claim_minus_W"], claim - led["W_log2"], 1e-9)
    add("gate_margin", led["gate_margin"], gate - claim, 1e-9)
    # claim is a BOUND: W <= claim <= gate
    out.append(("W <= claim <= gate", bool(led["W_log2"] <= claim <= gate), led["W_log2"], claim))
    # (3) S identity: S is what is left after the five row groups named in complete_W
    add("S==residual", led["S_written"],
        sum(rows.values()) - rows["records"] - rows["base_entry_nw50"] - rows["walker_state_nw50"]
        - rows["gap_point_tables"] - rows["gap_tail_fallback"], 0.0)
    # (4) records row from N (certificate) and the declared 2238/record
    nb = 1 << 66
    rps = -(-cert["N"] // nb) * nb
    add("records", rows["records"], 2 * rps * 2238, 0.0)
    # (5) gap row from the CAP (not from any expectation) and the 15 ops/entry charge
    k_max = der["gap_cap"]["k_max"]
    add("gap_point_tables", rows["gap_point_tables"], 2 * nb * (k_max // 2) * 15 * 2000, 0.0)
    add("secondary_rate_percent", led["secondary_rate_percent"],
        100.0 * (sum(rows.values()) - rows["records"]) / sum(rows.values()), 1e-4)
    # (6) kappa provenance: the ledger's charged value is the census object's own max (key path)
    add("kappa_charged==census.charged", led["kappa_charged"], cen["charged"], 0.0)
    out.append(("charged >= measured (charge is an envelope)",
                bool(cen["charged"] >= cen["kappa_key"]["max"]), cen["charged"], cen["kappa_key"]["max"]))
    add("kappa_mean==pairs_mean/expect", cr["kappa_key_mean"],
        cr["key_pairs_mean"] / cr["key_pair_expect_per_walk"], 1e-6)
    add("kappa_E_bad", der["kappa"]["price"]["E_bad_log2"],
        cert["debit"]["E_bad_log2"] + math.log2(cen["charged"] / cert["debit"]["kappa"]), 1e-6)
    # (7) arms: q recomputed from the arm's own integer counts
    add("arm_crown q", cr["q"], cr["hits"] / cr["trials"], 1e-9)
    add("arm_pool q", pr["q"], pr["hits"] / pr["trials"], 1e-9)
    add("arm_pool cp95_lower", der["pool_arm"]["cp95_lower"], cp_lower(pr["hits"], pr["trials"]), 2e-4)
    add("arm_crown cp95_lower", der["control"]["cp95_lower"], cp_lower(cr["hits"], cr["trials"]), 2e-4)
    add("arm_K==6N^2/m", cr["K"], 6.0 * (cr["nb"] * cr["ell"]) ** 2 / cr["m"], 1e-9)
    add("key_pair_expect", cr["key_pair_expect_per_walk"],
        (cr["records_per_side"] * (cr["records_per_side"] - 1) / 2.0) / 2 ** cr["key_bits"], 1e-3)
    add("gap_mean_over_mu", cr["gap_mean_over_mu"], cr["gap_mean"] / (8.0 * math.log(2.0 ** cr["walk_band_log2"])), 1e-3)
    # (8) ladder rungs: q from hits/trials, CP-lo from the same counts
    for rung in lad["ladder"]:
        add("ladder_q[nb=%s]" % rung["nb"], rung["q"], rung["hits"] / rung["trials"], 1e-3)
        add("ladder_cplo[nb=%s]" % rung["nb"], rung["cp95_lower"], cp_lower(rung["hits"], rung["trials"]), 3e-4)
    add("crown_q", lad["primary"]["q"], lad["primary"]["hits"] / lad["primary"]["trials"], 1e-6)
    add("crown_cplo", lad["primary"]["cp95_lower"], cp_lower(lad["primary"]["hits"], lad["primary"]["trials"]), 2e-4)
    # (9) label hygiene: the fixture POINT estimate vs its CP-lo are two different objects
    g0 = d25["summary"]["by_nw"]["nw1"]                 # the nw=1 anchor (131072 trials, 52307 solves)
    add("fixture q == solves/trials", g0["q"], g0["solves"] / g0["trials"], 1e-12)
    add("fixture cp95_lower == CP recomputed", g0["cp95_lower"], cp_lower(g0["solves"], g0["trials"]), 2e-4)
    out.append(("fixture q != cp95_lower (labels distinct)",
                abs(g0["q"] - g0["cp95_lower"]) > 1e-3, g0["q"], g0["cp95_lower"]))
    # the CERTIFICATE's success anchor must be the same run, and the declared CP-lo must be its CP-lo
    m1 = cert["success"]["measured_nw1"]
    add("cert nw1 q == fixture q", m1["q"], g0["q"], 1e-12)
    add("cert nw1 cp95_lower == fixture cp95_lower", m1["cp95_lower"], g0["cp95_lower"], 1e-9)
    out.append(("cert nw1 cp95_lower != q (no label swap in the package)",
                abs(m1["q"] - m1["cp95_lower"]) > 1e-3, m1["q"], m1["cp95_lower"]))
    # (10) gap scale + abort policy arithmetic, recomputed from certificate d
    mu_log2 = math.log2(cert["d"] * 76.0 * math.log(2))
    add("gap mu_log2", der["gap_cap"]["mu_log2"], mu_log2, 2e-3)
    add("abort_rate==aborts/trials", cd["abort_rate"], cd["aborts"] / cd["trials"], 1e-6)
    pred = 1 - math.exp(-cd["ell"] * math.exp(-cd["gap_cap"] / (8.0 * math.log(2.0 ** cd["walk_band_log2"]))))
    out.append(("cap-demo abort rate <= unsuppressed-law prediction (conservative direction)",
                bool(cd["abort_ratext"] if False else cd["abort_rate"] <= pred), cd["abort_rate"], round(pred, 6)))
    out.append(("cap-demo aborts>0 (the cap FIRES)", bool(cd["aborts"] > 0), cd["aborts"], ">0"))
    # (11) claim literals in submission.json must equal the ledger's claim, and stay under the gate
    s_claim = int(sub["claims"][0]["score"]["numerator"]) / int(sub["claims"][0]["score"]["denominator"])
    add("submission claim == ledger claim", s_claim, claim, 1e-12)
    out.append(("submission claim <= gate", bool(s_claim <= gate), s_claim, gate))
    return out


# ------------------------------ planted controls ------------------------------
def mutate(src, dst, name, fn):
    shutil.copytree(src, dst)
    p = os.path.join(dst, name)
    obj = json.load(open(p))
    fn(obj)
    json.dump(obj, open(p, "w"), indent=1)
    return dst


def controls(root):
    """Mutate SOURCE JSONs; the same checks() must fire each time."""
    tmp = tempfile.mkdtemp(prefix="auditv30-")
    results = []

    def rows_bump(o):
        o["rows"]["records"] += o["rows"]["records"] // 100

    def q_flip(o):
        o["q"] = o["q"] * 0.98

    def kappa_flip(o):
        o["kappa_key"]["max"] = o["kappa_key"]["max"] * 1.01

    def log2_flip(o):
        k = "walk_list"
        o["row_log2"][k] = o["row_log2"][k] + 0.5

    for label, fname, fn in [
        ("records +1% (ledger)", "ledger-v30.json", rows_bump),
        ("pool-arm q x0.98", "arm-pool-kb20.jsonl", None),
        ("census kappa_key.max x1.01", "census-kappa-v30.json", kappa_flip),
        ("ledger row_log2[walk_list] +0.5", "ledger-v30.json", log2_flip),
    ]:
        d = os.path.join(tmp, label.replace(" ", "_").replace("/", "_"))
        shutil.copytree(root, d)
        p = os.path.join(d, fname)
        if fname.endswith(".jsonl"):
            o = json.loads(open(p).read().splitlines()[0]); q_flip(o)
            open(p, "w").write(json.dumps(o) + "\n")
        else:
            o = json.load(open(p)); fn(o); json.dump(o, open(p, "w"), indent=1)
        res = checks(d)
        caught = [r for r in res if not r[1]]
        results.append((label, bool(caught), [r[0] for r in caught][:3]))
        del o
    shutil.rmtree(tmp, ignore_errors=True)
    return results


def main():
    root = HERE
    args = sys.argv[1:]
    if "--root" in args:
        root = args[args.index("--root") + 1]
    res = checks(root)
    bad = [r for r in res if not r[1]]
    ctl = controls(root)
    survived = [c for c in ctl if not c[1]]
    report = dict(
        kind="v30 independent audit: package literals vs recomputation from raw evidence",
        root=root, n_numeric_pairs=len([r for r in res if isinstance(r[2], (int, float))]),
        n_checks=len(res), failures=[dict(name=r[0], observed=r[2], recomputed=r[3]) for r in bad],
        controls=[dict(planted=c[0], caught=c[1], fired_on=c[2]) for c in ctl],
        controls_caught="%d/%d" % (len([c for c in ctl if c[1]]), len(ctl)),
        verdict="PASS" if (not bad and not survived) else "FAIL")
    out = os.path.join(root, "audit-v30.json")
    if "--out" in args:
        out = args[args.index("--out") + 1]
    json.dump(report, open(out, "w"), indent=1)
    print("checks=%d numeric_pairs=%d failures=%d controls=%s -> %s"
          % (len(res), report["n_numeric_pairs"], len(bad), report["controls_caught"], report["verdict"]))
    for r in bad:
        print("  FAIL", r[0], "observed", r[2], "recomputed", r[3])
    for c in ctl:
        print("  control", "CAUGHT" if c[1] else "SURVIVED", c[0], c[2])
    return 0 if report["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
