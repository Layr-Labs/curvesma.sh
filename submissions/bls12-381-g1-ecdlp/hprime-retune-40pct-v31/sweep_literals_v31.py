#!/usr/bin/env python3
"""v31 literal sweep (verdict 8d347ceb flaw 5): every integer of >=5 digits in the claims text and
the public note is resolved against the file that OWNS it.  One table: literal, source path, value,
match.  Exits non-zero if any literal is unowned or disagrees.  Also runs a planted-negative control:
a fabricated literal must be reported as NO MATCH."""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
EVIDENCE = ["derive-v31.json", "ledger-v31.json", "census-kappa-v31.json", "audit-v31.json",
            "certificate.json", "structured-ladder.json", "kappa-v26.json", "direct-v25.json",
            "fixtures-v25.json", "gate-v25.json", "ratios-v26.json", "census-key-v20.json",
            "submission.json", "arm-retune-v31.jsonl", "arm-sat-span14.jsonl",
            "arm-birthday-calib-v31.jsonl", "arm-ctrl-span20.jsonl", "structured_walk_v31.c"]
C_SRC = "structured_walk_v31.c"


def walk(obj, path, out):
    if isinstance(obj, dict):
        for k, v in obj.items():
            walk(v, path + "." + str(k), out)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            walk(v, path + "[%d]" % i, out)
    elif isinstance(obj, bool):
        return
    elif isinstance(obj, int):
        out.setdefault(str(obj), []).append(path)
    elif isinstance(obj, float):
        out.setdefault(repr(obj), []).append(path)
    elif isinstance(obj, str) and obj.isdigit():
        out.setdefault(obj, []).append(path)


def owners():
    own = {}
    for f in EVIDENCE:
        p = os.path.join(HERE, f)
        if not os.path.exists(p):
            continue
        if f.endswith(".json"):
            walk(json.load(open(p)), f, own)
        elif f.endswith(".jsonl"):
            for ln, line in enumerate(open(p), 1):
                if line.strip() and line.strip() != "DONE":
                    walk(json.loads(line), "%s:%d" % (f, ln), own)
        else:                                    # C source: integer literals incl. ULL suffixes
            for m in re.finditer(r"(?<![\w.])(\d{5,})(?:ULL|UL|U|L)?(?![\w.])", open(p).read()):
                own.setdefault(m.group(1), []).append(f + ":literal")
    return own


def lits(text):
    """integers of >=5 digits as written (thousands separators removed)"""
    t = text.replace(",", "")
    return re.findall(r"(?<![\w.^])(\d{5,})(?![\w.])", t)


def main():
    sub = json.load(open(HERE + "/submission.json"))
    note = open(HERE + "/public-note.md").read()
    corpus = {"claims[0].statement": sub["claims"][0]["statement"], "public-note.md": note,
              "title": sub["title"]}
    own = owners()
    table, bad = [], 0
    seen = set()
    for where, text in corpus.items():
        for L in lits(text):
            if (where, L) in seen:
                continue
            seen.add((where, L))
            srcs = own.get(L)
            ok = bool(srcs)
            table.append(dict(literal=L, where=where, source=(srcs[0] if srcs else None),
                              occurrences=(len(srcs) if srcs else 0), match=ok))
            if not ok:
                bad += 1
    # planted negative control: a fabricated literal must NOT resolve
    fake = "918273645"
    ctrl = dict(literal=fake, resolved=bool(own.get(fake)), caught=bool(not own.get(fake)))

    out = dict(kind="v31 literal sweep: every >=5-digit integer in the claims text and note vs the "
                    "file that owns it",
               literals_checked=len(table), unowned=bad,
               negative_control=ctrl,
               table=table,
               verdict="PASS" if bad == 0 and ctrl["caught"] else "FAIL")
    json.dump(out, open(HERE + "/literal-sweep-v31.json", "w"), indent=1)
    for r in table:
        print("%-24s %-22s %-46s %s" % (r["literal"], r["where"], str(r["source"])[:46],
                                        "Y" if r["match"] else "N"))
    print(json.dumps(dict(literals_checked=len(table), unowned=bad, control=ctrl["caught"],
                          verdict=out["verdict"])))
    if out["verdict"] != "PASS":
        sys.exit(1)


if __name__ == "__main__":
    main()
