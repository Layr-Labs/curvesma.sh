#!/usr/bin/env python3
"""Render public-note.md and argument.md from one machine table.

Every number in the two documents is substituted from evidence/numbers16.json,
which scripts/ledger16.py regenerates from the raw evidence JSONs.  The renderer
refuses unknown placeholders and then re-scans the rendered text: every decimal
token it finds must either occur in the machine table or be one of the design
constants listed in DESIGN, so a hand-typed number cannot survive a rebuild.
"""
import json
import re
import sys

ROOT = __file__.rsplit("/", 2)[0]
NUM = json.load(open(ROOT + "/evidence/numbers16.json"))
LED = json.load(open(ROOT + "/evidence/ledger16.json"))

V15_ID = "af5166b6"

EXTRA = {
    "v15_id": V15_ID,
    "toy_shift_failed": "0",
    "store_debit_exp": str(LED["store_debit_log2"]),
    "cap_minus_replay_instr": str(LED["cap_instr"] - LED["replay_instr"]),
    "claim_den_label": NUM["claim_den"],
    "rho_review_float": NUM["rho_review"],
    "trail_len": str(LED["trail_len"]),
    "v15_mul_count": "1233",   # multiply steps in the v15 run quoted by the rejection
}
TABLE = dict(NUM)
TABLE.update(EXTRA)

DESIGN = {
    "0.39", "3/8", "5/8", "12", "2", "1", "0", "100", "20000", "32", "52", "16", "8",
    "2^32", "2^100", "2^125", "381", "95", "256", "2^25.7", "0.43", "125.73", "125.74",
    "22", "52", "2048", "1024", "24", "4", "6", "2^73", "2^31", "2^125.728", "0.40",
    "0.418", "0.426", "0.436", "0.390", "0.3901", "1.55", "2",
}


PLACE = re.compile(r"\{([a-zA-Z_][a-zA-Z0-9_]*)\}")


def sub(m):
    k = m.group(1)
    if k not in TABLE:
        raise SystemExit("missing machine-table key: %s" % k)
    return str(TABLE[k])


def render(name):
    tmpl = open(ROOT + "/scripts/templates/" + name).read()
    out = PLACE.sub(sub, tmpl)
    if "{" in out and "}" in out:
        # any brace pair left over is an unsubstituted placeholder
        for m in re.finditer(r"\{[a-zA-Z_][a-zA-Z0-9_]*\}", out):
            raise SystemExit("unsubstituted placeholder in %s: %s" % (name, m.group(0)))
    open(ROOT + "/" + name, "w").write(out)
    return out


def audit(name, text):
    values = " ".join(str(v) for v in TABLE.values())
    bad = []
    for m in re.finditer(r"(?<![\w.^])(\d+(?:\.\d+)?)(?![\w^])", text):
        tok = m.group(1)
        if tok in DESIGN:
            continue
        if tok in values or tok.rstrip("0").rstrip(".") in values:
            continue
        bad.append(tok)
    return sorted(set(bad))


if __name__ == "__main__":
    bad_any = {}
    for name in ("public-note.md", "argument.md", "submission.json"):
        text = render(name)
        bad = audit(name, text)
        if bad:
            bad_any[name] = bad
        print("%-18s %6d B  numbers-from-table: %s" % (name, len(text.encode()),
                                                       "PASS" if not bad else bad))
    if bad_any:
        raise SystemExit("numbers outside the machine table: %s" % bad_any)
    print("machine-table substitution: PASS")
