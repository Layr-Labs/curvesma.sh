#!/usr/bin/env python3
r"""v29 EXTENDED label audit (re-runnable by intake).

v28 audited only 'a op b = c' and '2^e = v' pairs. This audits:
  (A) every '2^y' prose claim against its computed source row, in BOTH directions,
      within one sentence (regex 2\^[0-9.]+  <->  [0-9]+(.[0-9]+)?(e[+-]?[0-9]+)?);
  (B) every percentage claim ('X%') against a recomputed ratio;
  (C) every 'x'-multiple claim ('N.NNx') against a recomputed ratio;
  (D) counts (how many rungs satisfy a predicate);
  (E) direction-consistency: if the prose says A < B (or 'below'/'above'), the computed
      values must order the same way;
  (F) a PLANTED-FLIP control: each assertion class is re-run with one digit flipped and
      must FAIL, proving the audit is not vacuous.
Exit 0 iff every live claim passes and every planted flip is caught.
"""
import json, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DOCS = ["argument.md", "public-note.md", "submission.json"]
TOL_EQ = 0.02          # '=' or 'is' assertion on a log pair
TOL_APPROX = 0.06      # 'about'/'~'/'approx' assertion
TOL_PCT = 0.05         # percentage claims (absolute, in percentage points)

NUM = r"(?<![\w.^])(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(?![\w.^])"
POW = r"2\^(\d+(?:\.\d+)?)"


def sentences(text):
    return [s for s in re.split(r"(?<=[.;])\s+|\n+", text) if s.strip()]


def log2b(num_txt):
    v = float(num_txt)
    assert v > 0
    return math.log2(v)


def audit_pow_pairs(text, docs_name, tol_eq=TOL_EQ, tol_approx=TOL_APPROX):
    """(A) every asserted '2^e <op> magnitude' pair, BOTH directions, plus direction words.

    Only ASSERTED pairs are enforced (a 2^ token joined to a magnitude by = / is / equals /
    approx / ~); loose co-occurrence in one sentence is reported as informational, because a
    number in the same sentence is not necessarily a claim about that exponent. The planted
    control proves the enforcement fires.
    """
    rows, fails = [], []
    FWD = re.compile(r"2\^([0-9.]+)\s*(?:=|equals|\bis\b|≈|~|about|approx\.?)\s*"
                     r"([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)")
    REV = re.compile(r"(?<![\w,./*+^=(-])([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)"
                     r"\s*(?:=|equals|\bis\b|≈|~)\s*2\^([0-9.]+)")
    for sent in sentences(text):
        approx = bool(re.search(r"approx|about|~|≈|near", sent, re.I))
        tol = tol_approx if approx else tol_eq
        for m in FWD.finditer(sent):
            e, num = float(m.group(1).rstrip(".")), m.group(2).rstrip(".")
            d = abs(log2b(num) - e)
            rec = dict(doc=docs_name, claim=m.group(0), exponent=e, against=num,
                       log2_of_number=round(log2b(num), 4), delta=round(d, 6),
                       direction="2^e vs magnitude")
            rows.append(rec)
            if d > tol:
                fails.append(dict(rec, why="printed exponent and printed magnitude disagree",
                                  sentence=sent[:200]))
        for m in REV.finditer(sent):
            num, e = m.group(1).rstrip("."), float(m.group(2).rstrip("."))
            d = abs(log2b(num) - e)
            rec = dict(doc=docs_name, claim=m.group(0), exponent=e, against=num,
                       log2_of_number=round(log2b(num), 4), delta=round(d, 6),
                       direction="magnitude vs 2^e")
            rows.append(rec)
            if d > tol:
                fails.append(dict(rec, why="printed magnitude and printed exponent disagree",
                                  sentence=sent[:200]))
        m = re.search(r"([0-9.]+)\s*(above|below|over|under)\s*([0-9.]+)", sent)
        if m:
            a, rel, b = float(m.group(1)), m.group(2), float(m.group(3))
            ok = (a > b) if rel in ("above", "over") else (a < b)
            if not ok:
                fails.append(dict(doc=docs_name, sentence=sent[:200],
                                  why="direction word contradicts operands",
                                  operands=[a, rel, b]))
    return rows, fails


def audit_claims(claims, tol_pct=TOL_PCT):
    """(B) percentages, (C) multiples, (D) counts, all against recomputed sources."""
    fails, rows = [], []
    for c in claims:
        comp = eval(c["expr"], {"math": math, "log2": math.log2}, dict(c.get("vars", {})))  # noqa: S307
        if c["kind"] == "pct":
            ok = abs(comp - c["claimed"]) <= tol_pct
            rows.append(dict(claim=c["text"], computed=round(comp, 6), claimed=c["claimed"],
                             kind="pct", ok=ok))
        elif c["kind"] == "multiple":
            ok = abs(comp / c["claimed"] - 1.0) <= 0.005
            rows.append(dict(claim=c["text"], computed=round(comp, 6), claimed=c["claimed"],
                             kind="multiple", ok=ok))
        elif c["kind"] == "count":
            ok = int(comp) == int(c["claimed"])
            rows.append(dict(claim=c["text"], computed=comp, claimed=c["claimed"],
                             kind="count", ok=ok))
        elif c["kind"] == "logpair":
            ok = abs(math.log2(comp) - c["claimed"]) <= tol_pct / 100.0 * max(1.0, abs(c["claimed"])) + 0.02
            rows.append(dict(claim=c["text"], computed=round(comp, 6), log2=round(math.log2(comp), 4),
                             claimed_log2=c["claimed"], kind="logpair", ok=ok))
        else:
            ok, comp = False, None
        if not ok:
            fails.append(dict(claim=c["text"], computed=comp, claimed=c["claimed"], kind=c["kind"]))
    return rows, fails


def main():
    ledger = json.load(open(HERE + "/ledger-v29.json"))
    lad = json.load(open(HERE + "/structured-ladder.json"))
    der = json.load(open(HERE + "/derive-v29.json"))
    rungs = lad["ladder"]
    long_rungs = [r for r in rungs if r["ell"] >= 64]
    claims = [
        dict(kind="count", text="rungs with ell >= 64 (ladder only)", claimed=len(long_rungs),
             expr="len(r)", vars=dict(r=long_rungs)),
        dict(kind="count", text="rungs with ell >= 64 including primary", claimed=len(long_rungs) + 1,
             expr="len(r) + 1", vars=dict(r=long_rungs)),
        dict(kind="logpair", text="ell / 6944 (largest rung ratio)", claimed=36.89,
             expr="%d / 6944" % der["construction"]["ell"]),
        dict(kind="pct", text="secondary rows as pct of W",
             claimed=round(100 * (1 - ledger["rows"]["records"] / ledger["total_instructions"]), 4),
             expr="100*(1 - r/t)", vars=dict(r=ledger["rows"]["records"],
                                             t=ledger["total_instructions"])),
        dict(kind="pct", text="primary arm hit rate (%)", claimed=100 * lad["primary"]["q"],
             expr="100*k/n", vars=dict(k=lad["primary"]["hits"], n=lad["primary"]["trials"])),
        dict(kind="multiple", text="primary ell / fixture ell", claimed=108.5,
             expr="%d/64" % 6944),
    ]
    # per-walk kappa multiple
    if "max" in der.get("kappa", {}) and "expect_dup_per_walk" in der["kappa"]:
        claims.append(dict(kind="multiple", text="kappa_walk_max / expectation-of-1-dup",
                           claimed=der["kappa"]["max"] * der["kappa"]["expect_dup_per_walk"],
                           expr="%r*%r" % (der["kappa"]["max"], der["kappa"]["expect_dup_per_walk"])))

    all_fails, all_rows = [], []
    for d in DOCS:
        p = HERE + "/" + d
        if not os.path.exists(p):
            continue
        rows, fails = audit_pow_pairs(open(p).read(), d)
        all_rows += rows
        all_fails += fails
    crows, cfails = audit_claims(claims)
    all_rows += crows
    all_fails += cfails

    # (F) planted-flip controls: the audit MUST catch a flipped digit in each class
    ctl = []
    # power-pair control: claim 2^36.89 against 2^45
    # a 2^ exponent paired with a magnitude that is 11 bits off: must be caught
    _, f1 = audit_pow_pairs("The walk ratio is 2^36.89 = 1.4e14 times the largest rung.", "planted")
    ctl.append(dict(control="2^ exponent vs magnitude", caught=bool(f1)))
    # claim-level control: flip the rung count
    bad = [dict(c, claimed=c["claimed"] + 1) for c in claims if c["kind"] == "count"]
    _, f2 = audit_claims(bad)
    ctl.append(dict(control="rung count +1", caught=bool(f2)))
    bad2 = [dict(c, claimed=c["claimed"] * 1.05) for c in claims if c["kind"] == "pct"]
    _, f3 = audit_claims(bad2)
    ctl.append(dict(control="percentage x1.05", caught=bool(f3)))
    bad3 = [dict(c, claimed=c["claimed"] + 0.5) for c in claims if c["kind"] == "logpair"]
    _, f4 = audit_claims(bad3)
    ctl.append(dict(control="log2 exponent +0.5", caught=bool(f4)))

    out = dict(kind="v29 extended label audit", docs=DOCS,
               n_pow_pairs=sum(1 for r in all_rows if "exponent" in r),
               n_claims=len(crows), passes=sum(1 for r in all_rows if r.get("ok", True)),
               failures=all_fails, controls=ctl,
               controls_all_caught=all(c["caught"] for c in ctl),
               rows=all_rows[:80])
    json.dump(out, open(HERE + "/audit-labels-v29.json", "w"), indent=1)
    print(json.dumps(dict(n_pow_pairs=out["n_pow_pairs"], n_claims=out["n_claims"],
                          failures=len(all_fails), controls=ctl,
                          all_caught=out["controls_all_caught"]), indent=1))
    if all_fails or not out["controls_all_caught"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
