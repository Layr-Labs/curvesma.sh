#!/usr/bin/env python3
"""ARITHMETIC AUDIT (mandatory fleet row): extract every 'a op b = c' and '2^e = v'
pattern from the public note and argument, recompute each with exact Decimal, and
fail if any printed result disagrees beyond its own printed precision."""
import json, re, sys
from decimal import Decimal, getcontext

getcontext().prec = 60
FILES = sys.argv[1:] or ["public-note.md", "argument.md"]

NUM = r"-?\d[\d,]*(?:\.\d+)?"
SUB = re.compile(r"(?<![\w.^])(?P<a>" + NUM + r")\s*(?P<op>[+\-\u2212\u2013])\s*"
                 r"(?P<b>" + NUM + r")\s*=\s*(?P<c>" + NUM + r")")
POW = re.compile(r"2\^(?P<e>-?\d+(?:\.\d+)?)\s*=\s*(?P<v>\d[\d,]*)")


def dec_places(s):
    return len(s.split(".")[1]) if "." in s else 0


def audit(paths):
    rows, bad = [], []
    for p in paths:
        try:
            txt = open(p).read()
        except FileNotFoundError:
            continue
        for m in SUB.finditer(txt):
            a, b, c = (Decimal(m.group(k).replace(",", "")) for k in ("a", "b", "c"))
            got = a - b if m.group("op") in "-\u2212\u2013" else a + b
            tol = Decimal(5) * Decimal(10) ** (-(dec_places(m.group("c")) + 1))
            ok = abs(got - c) <= tol
            rows.append(dict(file=p, expr=m.group(0).strip(), a=str(a), op=m.group("op"),
                             b=str(b), printed=str(c), computed=str(got), ok=ok))
            if not ok:
                bad.append(rows[-1])
        for m in POW.finditer(txt):
            e = Decimal(m.group("e"))
            v = Decimal(m.group("v").replace(",", ""))
            got = Decimal(2) ** e
            ratio = float(got / v) if v else 0.0
            ok = 0.9995 <= ratio <= 1.0005
            rows.append(dict(file=p, expr=m.group(0).strip(), kind="pow2",
                             computed="%.4e" % float(got), printed=m.group("v"), ok=ok))
            if not ok:
                bad.append(rows[-1])
    return rows, bad


if __name__ == "__main__":
    rows, bad = audit(FILES)
    out = dict(kind="v28 arithmetic audit (regex extraction, exact Decimal recomputation)",
               files=FILES, patterns_checked=len(rows), failures=len(bad),
               rule="every 'a op b = c' in the public text must recompute to its printed "
                    "precision; every '2^e = v' must agree within 5e-4 relative",
               rows=rows)
    json.dump(out, open("arithmetic-audit.json", "w"), indent=1)
    for r in rows:
        print(("OK  " if r["ok"] else "FAIL"), r["expr"], "->", r.get("computed"))
    print("\nchecked %d, failures %d" % (len(rows), len(bad)))
    raise SystemExit(1 if bad else 0)
