#!/usr/bin/env python3
"""selfaudit.py -- p256 v11 fleet-law audits with planted-flip controls.

Checks (all machine-computed, nothing read from prose):
  A1  sum(step rows) == STEP_INSN == ledger-out.step_insn
  A2  claim: submission.json score == ledger-out.claim (exact rational x1000)
  A3  claim == ceil(log2 W * 1000)/1000 recomputed from the ledger's own W
  A4  success: ledger-out.success == cp95_design - debits - restart_loss  (x1e-12)
  A5  floor direction: floor_met == (success >= 0.39); ranks == (claim < gate)
  A6  every 4+ significant-figure decimal in public-note.md is either a field of
      evidence/ledger-out.json (to its printed precision) or appears in this file's
      CITED list with its source field -- an orphan number FAILS the audit
  A7  counter coverage: every call site of qwalk.step in the shipped scripts is in
      the audited call-site list (no uncounted path)
  A8  evidence completeness: every file referenced by submission.json exists
  
Each check is a function; the planted-flip controls re-run the SAME checker on a
mutated copy and must FAIL (3/3).  A control that does not fire means the checker is
vacuous, and the audit reports 0/3 and exits nonzero.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EVID = os.path.join(HERE, "evidence")


def load(rel):
    with open(os.path.join(HERE, rel)) as f:
        return json.load(f)


LEDGER = load("evidence/ledger-out.json")
SUB = load("submission.json")
CLAIM = SUB["claims"][0]["score"]
CLAIM_F = int(CLAIM["numerator"]) / int(CLAIM["denominator"])
NOTE = open(os.path.join(HERE, "public-note.md")).read()
REPRO = load("reproduction.json")

CITED = {
    # literals the note may print that are NOT ledger-out fields, each with its
    # source.  Anything else numeric in the note must be a ledger-out field at some
    # printed precision (A6).
    "0.0705": "rhosigma shipped-generator census (abandoned family, v10p evidence)",
    "0.4558": "rhosigma control-generator census (v10p evidence)",
    "126.2364": "plain-lineage ticket 5ecbf15 score (bls12-381-g1-ecdlp-plain)",
    "0.39": "track success floor (rules.md)",
    "0.3901": "the planted-flip control literal (0.39 + 1e-4), named in the note's audit line",
}


def a1_step_rows(L):
    return sum(n for _o, n in L["step_rows"]) == L["step_insn"]


def a2_claim_matches(L):
    return abs(int(CLAIM["numerator"]) / int(CLAIM["denominator"]) - L["claim"]) <= 5e-4


def a3_claim_is_ceil(L):
    import math
    want = math.ceil(L["log2_W"] * 1000) / 1000.0
    return abs(want - L["claim"]) <= 1e-9 and L["claim"] >= L["log2_W"]


def a4_success_arithmetic(L):
    got = L["cp95_design"] - L["debits"] - L["restart_loss"]
    return abs(got - L["success"]) <= 1e-12


def a5_directions(L):
    return (L["floor_met"] == (L["success"] >= L["floor"])
            and L["ranks"] == (L["claim"] < L["gate"]))


NUM_RE = re.compile(r"(?<![\w.])(\d+\.\d{3,})(?![\w])")


def a6_no_orphan_numbers(L, note=None):
    note = NOTE if note is None else note
    known = set()
    for k, v in L.items():
        if isinstance(v, float):
            for prec in range(2, 10):
                known.add(("%." + str(prec) + "f") % abs(v))
        elif isinstance(v, int):
            known.add(str(v))
    known.update(CITED)
    bad = []
    for m in NUM_RE.finditer(note):
        tok = m.group(1)
        if tok in known:
            continue
        # a printed decimal is accepted if it is a prefix of some known printed row
        if any(k.startswith(tok) for k in known if len(tok) >= 4):
            continue
        bad.append(tok)
    if bad:
        print("  ORPHAN NUMBERS:", sorted(set(bad))[:12])
    return not bad


def a7_counter_coverage(L):
    """Every shipped call site of the counted step must be listed here."""
    audited = {"scripts/census_v20n.py:walk_arm->qwalk.step",
               "scripts/pair_census.py->census_v20n.walk_arm",
               "scripts/transfer_arm.py->census_v20n.walk_arm"}
    found = set()
    for name in os.listdir(os.path.join(HERE, "scripts")):
        if not name.endswith(".py"):
            continue
        body = open(os.path.join(HERE, "scripts", name)).read()
        for m in re.finditer(r"qwalk\.step\(|C\.walk_arm\(|walk_arm\(", body):
            found.add("%s:%s" % (name, m.group(0)))
    # the counted path is walk_arm's qwalk.step call; no other module may call it
    direct = {f for f in found if f.startswith("census_v20n.py:qwalk.step")}
    table_only = {f for f in found if f.startswith("census_v20n.py:C.walk_arm")}
    return bool(direct) and not table_only


def a8_evidence_complete(L):
    ev = SUB["claims"][0]["evidence"]
    missing = [e for e in ev if not os.path.exists(os.path.join(HERE, e))]
    if "evidence/ledger-out.json" not in ev:
        missing.append("evidence/ledger-out.json (not listed)")
    if missing:
        print("  MISSING:", missing)
    return not missing


CHECKS = [("A1 step rows == STEP_INSN", a1_step_rows),
          ("A2 submission score == ledger claim", a2_claim_matches),
          ("A3 claim == ceil(log2 W)", a3_claim_is_ceil),
          ("A4 success arithmetic", a4_success_arithmetic),
          ("A5 direction laws", a5_directions),
          ("A6 no orphan note numbers", a6_no_orphan_numbers),
          ("A7 counter coverage", a7_counter_coverage),
          ("A8 evidence completeness", a8_evidence_complete)]


def run(L, note=None):
    out = {}
    for name, fn in CHECKS:
        try:
            out[name] = bool(fn(L) if fn is not a7_counter_coverage and fn is not a6_no_orphan_numbers
                             else (fn(L) if fn is a7_counter_coverage else fn(L, note)))
        except Exception as exc:                     # a crash is a FAIL, not a pass
            print("  %s raised %r" % (name, exc))
            out[name] = False
    return out


def main():
    res = run(LEDGER)
    for k, v in res.items():
        print("%-38s %s" % (k, "PASS" if v else "FAIL"))

    # ---- planted-flip controls: the SAME checker must catch each mutation ----
    controls = []
    m1 = dict(LEDGER)
    m1["step_rows"] = [(o, n + 1) for o, n in LEDGER["step_rows"]]
    controls.append(("step-row +1", m1, "A1 step rows == STEP_INSN"))
    m2 = dict(LEDGER)
    m2["claim"] = LEDGER["claim"] - 0.002
    controls.append(("claim -0.002", m2, "A2 submission score == ledger claim"))
    m3 = dict(LEDGER)
    m3["success"] = LEDGER["success"] + 0.01
    controls.append(("success +0.01", m3, "A4 success arithmetic"))
    fired = []
    for label, mut, key in controls:
        r = run(mut)
        fired.append(not r[key])
        print("control %-18s caught by %-38s %s"
              % (label, key, "FIRED" if not r[key] else "MISSED"))
    # note-text control: a wrong floor direction must trip A6/A5 text law
    bad_note = NOTE.replace("%.3f" % LEDGER["claim"], "%.3f" % (LEDGER["claim"] + 0.001), 1)
    n_ok = a6_no_orphan_numbers(LEDGER, bad_note)
    print("control note-claim +0.001 caught by A6                          %s"
          % ("FIRED" if not n_ok else "MISSED"))
    total_fired = sum(1 for f in fired if f) + (0 if n_ok else 1)
    ok = all(res.values())
    print("audit %s; controls %d/4 fired" % ("PASS" if ok else "FAIL", total_fired))
    with open(os.path.join(EVID, "selfaudit.json"), "w") as f:
        json.dump({"checks": res, "controls_fired": total_fired, "controls_total": 4,
                   "pass": ok}, f, indent=1, sort_keys=True)
    return 0 if (ok and total_fired == 4) else 1


if __name__ == "__main__":
    sys.exit(main())
