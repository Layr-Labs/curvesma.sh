#!/usr/bin/env python3
"""sensitivity_assert.py -- v19: MACHINE-VERIFY every numeric sensitivity assertion.

Verdict 46b1b85b found a hand-checkable assertion error: the package declared that
+-5 on M_insn "moves W by 2^-113 bits", while the reviewer computed ~6e-6 bits from
batch_inversion_total ~ (265*M_insn/1000)*S_OPS/2^18.  The assertion was ANALYTIC
and WRONG.  The repair is not a better analytic claim: it is to PERTURB the ledger,
recompile, and report what actually moves.

This script executes ledger.py's module body three times -- at M_insn = 97 (shipped),
92 and 102 -- in a private namespace, and reports the exact change in the submitted
work bound (log2 of W2 = 2*W), in bits, together with the changed EXTRAS rows.

It also SCANS the shipped text for sensitivity-style assertions ("insensitive",
"negligible", "moves W by", "does not depend on", "2^-1xx") so that each one is
either backed by a measured row in evidence/sensitivity-assert.json or DELETED.
Nothing in the note may quote a sensitivity that is not in that file.
"""
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
LEDGER = os.path.join(D, "ledger.py")
OUT = os.path.join(D, "evidence", "sensitivity-assert.json")
SRC = open(LEDGER).read()
assert re.search(r"^M_insn = 97$", SRC, re.M)


def run(m_insn):
    s = re.sub(r"^M_insn = 97$", "M_insn = %d" % m_insn, SRC, count=1, flags=re.M)
    ns = {"__file__": LEDGER, "__name__": "ledger_variant"}
    exec(compile(s, LEDGER, "exec"), ns)          # noqa: S102 - deliberate
    return ns


base = run(97)
LOGW0 = base["LOG_W"]
rows = []
for delta in (-5, 0, 5):
    v = run(97 + delta)
    dlogw = v["LOG_W"] - LOGW0
    rows.append({
        "M_insn": 97 + delta,
        "delta_insn": delta,
        "log2_W2": v["LOG_W"],
        "delta_log2_W2_bits": dlogw,
        "delta_log2_W2_bits_str": "%.6e" % dlogw,
        "claim_log2_at_1e-3": v["CLAIM_LOG2"],
        "EXTRAS_ratio_change": float(v["EXTRAS"] / base["EXTRAS"] - 1),
        "table_build_ref_change": float(v["TABLE_BUILD_REF"] / base["TABLE_BUILD_REF"] - 1),
        "inv_insn_change": float(v["INV_INSN"] / base["INV_INSN"] - 1),
        "batch_inversion_change": float(v["batch_inversion_total"]
                                        / base["batch_inversion_total"] - 1),
        "fallback_ref_change": float(v["FALLBACK_REF"] / base["FALLBACK_REF"] - 1),
    })

# the analytic form the reviewer used, for the side-by-side row
M0 = 97
inv0 = base["batch_inversion_total"]
claim_old = "2^-113"
max_abs = max(abs(r["delta_log2_W2_bits"]) for r in rows)
# per-instruction slope of log2 W2, exact enough by central difference
slope = (rows[2]["delta_log2_W2_bits"] - rows[0]["delta_log2_W2_bits"]) / 10.0

# ---- scan the shipped text for sensitivity-style assertions -------------------
PAT = re.compile(r"(insensitive|negligible|negligibly|moves? W by|does not depend on|"
                 r"depend\w* on it|2\^-1[0-9]{2}|within rounding|to within)", re.I)
TEXT_EXT = (".md", ".py", ".json", ".txt")
hits = []
for root, _dirs, files in os.walk(D):
    if "/evidence" in root or "/.git" in root:
        continue
    for fn in files:
        if not fn.endswith(TEXT_EXT):
            continue
        rel = os.path.relpath(os.path.join(root, fn), D)
        if rel.startswith("evidence/") or rel == "scripts/sensitivity_assert.py":
            continue
        try:
            txt = open(os.path.join(root, fn)).read()
        except (UnicodeDecodeError, OSError):
            continue
        for ln, line in enumerate(txt.splitlines(), 1):
            if PAT.search(line):
                hits.append({"file": rel, "line": ln, "text": line.strip()[:240]})

out = {
    "test": "sensitivity_assert_machine_verified",
    "what": ("M_insn +-5 perturbation recomputes the ledger with no analytic shortcut; "
             "the declared sensitivity is replaced by the measured delta"),
    "baseline": {"M_insn": 97, "log2_W2": LOGW0, "claim_log2": base["CLAIM_LOG2"]},
    "rows": rows,
    "max_abs_delta_bits": max_abs,
    "slope_bits_per_insn": slope,
    "declared_in_v18": claim_old,
    "declared_is_false": True,
    "declared_value_is_smaller_than_measured": bool(max_abs > 2.0 ** -113),
    "measured_value_2_pow": ("2^%.2f" % math.log2(max_abs)) if max_abs > 0 else "0",
    "claim_log2_moves": bool(any(r["claim_log2_at_1e-3"] != base["CLAIM_LOG2"] for r in rows)),
    "assertion_scan": {"pattern": PAT.pattern, "hits": hits, "n_hits": len(hits)},
}
with open(OUT, "w") as f:
    json.dump(out, f, indent=1, sort_keys=True)
print("M_insn perturbation (exact, recompiled):")
for r in rows:
    print("  M_insn %3d (%+d): log2 W2 = %.8f  delta %+.6e bits  claim %.3f"
          % (r["M_insn"], r["delta_insn"], r["log2_W2"], r["delta_log2_W2_bits"],
             r["claim_log2_at_1e-3"]))
print("max |delta| = %.4e bits = %s (v18 declared %s -> FALSE)"
      % (max_abs, out["measured_value_2_pow"], claim_old))
print("claim (1e-3) moves: %s" % out["claim_log2_moves"])
print("assertion scan: %d hits" % len(hits))
for h in hits:
    print("  %s:%d  %s" % (h["file"], h["line"], h["text"]))
