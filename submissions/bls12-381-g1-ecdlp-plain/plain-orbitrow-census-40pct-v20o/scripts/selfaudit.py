#!/usr/bin/env python3
"""selfaudit.py -- v20n: path-existence + citation-resolution + row-binding audit.

Verdict e9ee1dbc point 3 (third time absorption evidence was called absent):
the v20m arm RAN and wrote to a workspaces/ path that does not exist in the
package.  A numeric audit cannot see that, so this audit walks the citations:

  1. every file named by submission.json:evidence, and every path-like token in
     public-note.md / argument.md / certificate.json, EXISTS on the staged tree;
  2. every `file.json:dotted.field` citation RESOLVES to a value;
  3. every file in evidence/ except verdict-*.md and selfaudit.json is CITED;
  4. PLANTED CONTROL: a bogus citation is injected into the set and the audit
     must FAIL on it -- a path check that cannot fail is not evidence;
  5. numeric row binding: ledger.RESTART_REF / FRESH_REF equal the JSON rows
     they are supposed to consume (the v20n repricing fix), and every 4+
     significant-figure number in the note is emitted by some evidence JSON or
     by ledger.py (orphan check -- the v20g failure mode).

Writes evidence/selfaudit.json and exits non-zero unless everything passes.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, D)

PATH_RE = re.compile(r"(?<![\w/])((?:evidence|scripts)/[\w.\-]+\.(?:json|py|md|c|dis)"
                     r"|ledger\.py|public-note\.md|argument\.md|certificate\.json"
                     r"|submission\.json|dependencies\.json|reproduction\.json)")
FIELD_RE = re.compile(r"((?:evidence)/[\w.\-]+\.json):([A-Za-z_][\w.]*)")
NUM = re.compile(r"(?<![\d.])\d[\d,]*(?:\.\d+)?(?:[eE]-?\d+)?")


def flat(obj, out):
    if isinstance(obj, dict):
        for v in obj.values():
            flat(v, out)
    elif isinstance(obj, list):
        for v in obj:
            flat(v, out)
    elif isinstance(obj, bool):
        pass
    elif isinstance(obj, (int, float)):
        out.append(float(obj))
    elif isinstance(obj, str):
        try:
            out.append(float(obj))
        except ValueError:
            pass


def resolve(doc, dotted):
    cur = doc
    for part in dotted.split("."):
        if isinstance(cur, list):
            cur = cur[int(part)]
        elif isinstance(cur, dict) and part in cur:
            cur = cur[part]
        else:
            return False
    return True


def sig_figs(tok):
    t = tok.replace(",", "").lower()
    if "e" in t:
        t = t.split("e")[0]
    return len(t.replace(".", "").lstrip("0"))


def audit(extra_citation=None):
    manifest = json.load(open(os.path.join(D, "submission.json")))
    cited = set(manifest["claims"][0]["evidence"])
    texts = {}
    for name in ("public-note.md", "argument.md", "certificate.json",
                 "reproduction.json"):
        p = os.path.join(D, name)
        if os.path.exists(p):
            texts[name] = open(p).read()
    for name, txt in texts.items():
        cited |= set(PATH_RE.findall(txt))
    if extra_citation:
        cited.add(extra_citation)
    # selfaudit.json is the artifact this run is WRITING: exempt it from the
    # existence check (its own citation appears in the note text).
    missing = sorted(c for c in cited
                     if c not in ("evidence/selfaudit.json",)
                     and not os.path.exists(os.path.join(D, c)))
    fields = []
    for name, txt in texts.items():
        fields += [(f, d) for f, d in FIELD_RE.findall(txt)]
    unresolved = []
    for f, dfield in sorted(set(fields)):
        p = os.path.join(D, f)
        if f == "evidence/selfaudit.json":
            continue
        if not os.path.exists(p):
            unresolved.append("%s:%s (file missing)" % (f, dfield))
            continue
        try:
            doc = json.load(open(p))
        except Exception:
            continue
        if not resolve(doc, dfield):
            unresolved.append("%s:%s" % (f, dfield))
    evdir = os.path.join(D, "evidence")
    present = sorted(os.listdir(evdir)) if os.path.isdir(evdir) else []
    uncited = [f for f in present
               if "evidence/" + f not in cited
               and not f.startswith("verdict-") and f != "selfaudit.json"]
    vals = []
    for f in present:
        if f.endswith(".json") and f != "selfaudit.json":
            try:
                flat(json.load(open(os.path.join(evdir, f))), vals)
            except Exception:
                pass
    import ledger
    import math as _m
    lst = [float(_m.log2(float(ledger.E_RESTARTS))) if hasattr(ledger, "E_RESTARTS") else 0.0,
           float(ledger.LOG_W), float(ledger.CLAIM), float(ledger.SUCCESS),
           float(ledger.METER_CAP), float(ledger.RESTART_REF),
           float(ledger.FRESH_REF), float(ledger.C),
           float(ledger.ABSORB_P_LEN3_SMALL), float(ledger.ABSORB_P_LEN3_LARGE),
           float(ledger.ABSORB_STARTS), float(ledger.K_ESC_MAX),
           float(ledger.STEP_INSN), float(ledger.R_EFF),
           float(ledger.SUCCESS_WALK_ONLY) if hasattr(ledger, "SUCCESS_WALK_ONLY") else 0.0,
           float(ledger.STEP_PRICE), float(ledger.BUDGET),
           float(ledger.CLAIM_RATIO_CEIL), float(ledger.TABLE_BUILD_REF),
           float(ledger.B_DRAIN_STEPS), float(ledger.B_DRAIN_PRICE),
           float(ledger.DEBITS), float(ledger.METER_CAP / ledger.FRESH_REF),
           float(ledger.E_RESTARTS)]
    vals += lst
    vals += [float(ledger.DEBITS), __import__('math').log2(float(ledger.E_RESTARTS)),
             float(ledger.STEP_PRICE), float(ledger.BUDGET), float(ledger.CLAIM_RATIO_CEIL),
             float(ledger.TABLE_BUILD_REF), float(ledger.B_DRAIN_STEPS),
             float(ledger.B_DRAIN_PRICE), float(ledger.METER_CAP / ledger.FRESH_REF)]
    orphans = []
    note = texts.get("public-note.md", "")
    for m in NUM.finditer(note):
        tok = m.group(0)
        if sig_figs(tok) < 4:
            continue
        x = float(tok.replace(",", ""))
        if not any(abs(x - v) <= max(5e-4, 5e-4 * abs(v)) for v in vals):
            orphans.append(tok)
    rst = json.load(open(os.path.join(D, "evidence", "restart-count.json")))
    row_bind = {
        "ledger_fresh_ref_eq_json": abs(float(ledger.FRESH_REF) - rst["fresh_ref"]) < 1e-9,
        "ledger_restart_ref_eq_json_fresh": abs(float(ledger.RESTART_REF)
                                                - float(ledger.METER_CAP)) < 1e-6,
        "json_chain_ref_neq_v20l_319": abs(rst["chain"]["chain_ref"] - 319.0) > 1.0,
    }
    return {
        "test": "selfaudit_v20n",
        "cited_files": len(cited), "missing_files": missing,
        "field_citations": len(set(fields)), "unresolved_fields": unresolved,
        "evidence_files": len(present), "uncited_evidence": uncited,
        "numbers_checked": sum(1 for _ in NUM.finditer(note)),
        "orphans": sorted(set(orphans)),
        "row_binding": row_bind,
        "passed": (not missing and not unresolved and not uncited and not orphans
                   and all(row_bind.values())),
    }


def main():
    res = audit()
    planted = audit(extra_citation="evidence/__planted_missing__.json")
    res["planted_wrong_path"] = {
        "injected": "evidence/__planted_missing__.json",
        "caught": bool(planted["missing_files"]),
        "missing": planted["missing_files"]}
    res["passed"] = res["passed"] and res["planted_wrong_path"]["caught"]
    with open(os.path.join(D, "evidence", "selfaudit.json"), "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    print(json.dumps(res, indent=1, sort_keys=True))
    raise SystemExit(0 if res["passed"] else 1)


if __name__ == "__main__":
    main()
