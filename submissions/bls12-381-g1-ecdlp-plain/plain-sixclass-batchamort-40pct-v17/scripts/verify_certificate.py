#!/usr/bin/env python3
"""verify_certificate.py -- v17 self-check.  No human arithmetic is trusted.

Every check recomputes from a raw source: ledger.py's exact Fractions, the printed
lines of argument.md, evidence/toy.json's counters, or the exact integer test
2^(c-1) < W^1000 <= 2^c for the claim.  Rounded display strings are never compared
against rounded display strings; where a document prints a rounded value, the
checker re-renders it exactly from the Fraction and compares digit for digit.

Writes certificate.json.  Exit status 0 iff every check passes.
"""
import ast
import json
import os
import random
import re
import sys
from decimal import Decimal
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, PKG)
sys.path.insert(0, HERE)
import ledger      # noqa: E402
import selector    # noqa: E402
import gen_tables  # noqa: E402

CHECKS = []
GATE = F(126288, 1000)          # board gate pulled at submit time, see public-note
EVIDENCE_BUDGET = 512 * 1024


def add(name, ok, detail):
    CHECKS.append({"check": name, "pass": bool(ok), "detail": detail})
    return ok


def frac(text):
    """Parse an exact printed value: '1.0665', '101495/262144', '133'."""
    return F(text.strip())


def eval_exact(expr):
    """Evaluate a printed formula in EXACT rational arithmetic.

    Python's `/` on ints returns a float, which would silently turn `133/2000`
    into a binary approximation; this walks the AST and performs every division
    as a Fraction division, so the comparison against the printed value is exact.
    """
    node = ast.parse(expr, mode="eval").body

    def ev(n):
        if isinstance(n, ast.BinOp):
            a, b = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add):
                return a + b
            if isinstance(n.op, ast.Sub):
                return a - b
            if isinstance(n.op, ast.Mult):
                return a * b
            if isinstance(n.op, ast.Div):
                return F(a) / F(b)
            if isinstance(n.op, ast.Pow):
                return F(a) ** int(b)
            raise ValueError("operator not allowed: %r" % type(n.op))
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            return -ev(n.operand)
        if isinstance(n, ast.Constant) and isinstance(n.value, int):
            return F(n.value)
        raise ValueError("token not allowed: %r" % (n,))

    return F(ev(node))


# ---- 1: the shared selector, against exact integer division -----------------
FX_P, FX_Q = 600217, 50047
random.seed(20261023)
keys = [random.randrange(0, FX_P // 2 + 1) for _ in range(100)]   # class keys live in [0, P//2]
def ref_cell(K, R, P):
    """The declared map, re-derived here: NO import, NO shared alias."""
    return (K * R) >> (P.bit_length() - 1)


bad = [K for K in keys if selector.cell(K, 4096, FX_P) != ref_cell(K, 4096, FX_P)]
add("selector_100_random_keys", not bad,
    {"keys": 100, "agree": 100 - len(bad), "R": 4096, "P": FX_P,
     "map": "scripts/selector.py::cell == (K*R) >> (bitlen P - 1), reference re-derived in the checker"})

shift_v = ref_cell(2045, 4096, FX_P)
quot_v = (2045 * 4096) // FX_P
add("selector_is_the_declared_shift_not_floor_KR_over_P",
    shift_v == 15 and quot_v == 13 and selector.cell(2045, 4096, FX_P) == shift_v,
    {"K": 2045, "R": 4096, "declared_shift": shift_v, "floor_KR_over_P": quot_v,
     "shared": selector.cell(2045, 4096, FX_P),
     "note": "the two maps differ by design; the package declares and ships the shift, "
             "and R_eff is the reachable-column count that makes its column law uniform",
     "v15_value_the_old_toy_computed": 3572})

edge_bad, edge_n = [], 0
for Rv in (256, 512, 1024, 2048, 4096):
    for base in (0, 1, FX_P // Rv, FX_P // 2 - 1, FX_P - 1):
        for dk in range(-4, 5):
            K = base + dk
            if 0 <= K <= FX_P // 2:
                edge_n += 1
                if selector.cell(K, Rv, FX_P) != selector.cell_exact(K, Rv, FX_P):
                    edge_bad.append((K, Rv))
add("selector_boundary_neighbourhoods", not edge_bad,
    {"keys": edge_n, "R_swept": [256, 512, 1024, 2048, 4096], "mismatches": len(edge_bad)})

add("selector_cell_within_R", all(selector.cell(K, 2048, FX_P) < 2048 for K in keys),
    {"rule": "K <= P//2 => cell < R, no mask needed"})
add("reachable_columns_matches_definition",
    ledger.R_EFF == (((ledger.p // 2) * ledger.R_TAB) >> (ledger.p.bit_length() - 1)) + 1,
    {"R_eff": ledger.R_EFF, "R": ledger.R_TAB, "P_bits": ledger.p.bit_length()})
rem = ((ledger.p // 2) * ledger.R_TAB) % (1 << (ledger.p.bit_length() - 1))
add("last_column_is_the_only_partial_one", 0 < rem < (1 << (ledger.p.bit_length() - 1)),
    {"remainder_bits": rem.bit_length(),
     "note": "reachable K-range is a prefix, so every column but the last is equally weighted"})

# ---- 2: the document is generated, not typed --------------------------------
with open(os.path.join(PKG, "argument.md")) as f:
    doc = f.read()
block = doc[doc.index(gen_tables.BEGIN):doc.index(gen_tables.END) + len(gen_tables.END)]
add("argument_block_matches_render", block.rstrip("\n") == gen_tables.render().rstrip("\n"),
    {"renderer": "scripts/gen_tables.py", "source": "ledger.py", "bytes": len(gen_tables.render())})

rows = re.findall(r"^\| ([a-z0-9_]+) \| .* \| (\d+) \|$", block, re.M)
row_map = {n: int(i) for n, i in rows}
add("printed_rows_readd_to_step_insn",
    sum(int(i) for _n, i in rows) == ledger.STEP_INSN
    and row_map == {n: v for n, _o, v in ledger.ROWS},
    {"rows_printed": len(rows), "sum": sum(int(i) for _n, i in rows),
     "ledger_step_insn": ledger.STEP_INSN})

# ---- 3: every F: line is evaluated, not believed ----------------------------
flines = re.findall(r"^F: (.+?) = (\S+)\s+<!-- ([a-z0-9_]+) -->$", block, re.M)
f_bad, f_seen = [], []
for lhs, rhs, tag in flines:
    try:
        got = eval_exact(lhs)                            # exact rational arithmetic only
    except Exception as exc:                              # noqa: BLE001
        f_bad.append({"tag": tag, "error": repr(exc)})
        continue
    want = frac(rhs)
    f_seen.append(tag)
    if got != want:
        f_bad.append({"tag": tag, "lhs": lhs, "got": str(got), "printed": rhs})
add("every_printed_formula_evaluates", not f_bad and len(flines) >= 8,
    {"formulas": len(flines), "tags": f_seen, "failures": f_bad})

# ---- 4: the ledger's own arithmetic, recomputed from its integers ----------
add("inversion_instructions", ledger.INV_INSN == (254 + 11) * ledger.M_insn,
    {"INV_INSN": ledger.INV_INSN, "mults": 254 + 11, "M": ledger.M_insn})
add("inversion_share_per_step", ledger.INV_PER_STEP == F(ledger.INV_INSN, ledger.BATCH),
    {"INV_PER_STEP": gen_tables.fmt_exact(ledger.INV_PER_STEP), "B": ledger.BATCH})

want_batches = -((-(ledger.S_ops + 2 ** 60)) // ledger.BATCH)
add("batch_count_is_a_ceiling",
    ledger.n_batches == want_batches and ledger.n_batches + 1 > want_batches,
    {"n_batches": ledger.n_batches, "ceil": want_batches, "priced": ledger.n_batches + 1,
     "note": "one whole tail batch is charged"})
add("batch_inversion_total_exact",
    ledger.batch_inversion_total == F(ledger.INV_INSN * (ledger.n_batches + 1), ledger.C),
    {"value": gen_tables.fmt_exact(ledger.batch_inversion_total), "units": "reference"})
add("phase1_exact", ledger.PHASE1 == ledger.STEP_PRICE * (ledger.S_ops + 2 ** 60),
    {"step_price": gen_tables.fmt_exact(ledger.STEP_PRICE)})
add("W_is_the_sum_of_its_three_terms",
    ledger.W == ledger.PHASE1 + ledger.EXTRAS + ledger.batch_inversion_total,
    {"W": gen_tables.log2_8(ledger.W),
     "phase1": gen_tables.log2_8(ledger.PHASE1),
     "extras": gen_tables.log2_8(ledger.EXTRAS),
     "batch_inversion": gen_tables.log2_8(ledger.batch_inversion_total)})
add("step_price_above_doubling_price", ledger.DOUBLING_PRICE < ledger.STEP_PRICE,
    {"step": gen_tables.fmt_exact(ledger.STEP_PRICE),
     "doubling": gen_tables.fmt_exact(ledger.DOUBLING_PRICE)})
add("S_ops_is_ceil_A_H1_S0",
    ledger.S_ops == -((-ledger.A_H1 * ledger.S0).numerator // (ledger.A_H1 * ledger.S0).denominator)
    and F(ledger.S_ops) >= ledger.A_H1 * ledger.S0,
    {"S_ops": ledger.S_ops, "A_H1*S0 <= S_ops": True})

# ---- 5: claim tightness in both directions, exact integer powers -----------
claim = ledger.CLAIM
c_num = claim.numerator * 1000 // claim.denominator      # claim in thousandths
Wp = ledger.W ** 1000
add("claim_is_the_ceil_of_log2_W_both_directions",
    2 ** (c_num - 1) < Wp <= 2 ** c_num,
    {"claim": gen_tables.fmt_exact(claim), "W^1000 <= 2^claim": Wp <= 2 ** c_num,
     "2^(claim-0.001) < W^1000": 2 ** (c_num - 1) < Wp, "log2 W": gen_tables.log2_8(ledger.W)})
add("claim_is_tightest_in_its_granularity",
    not (2 ** (c_num - 2) >= Wp),
    {"note": "two steps coarser would break the upper bound"})

# ---- 6: window row ---------------------------------------------------------
log2w8 = Decimal(gen_tables.log2_8(ledger.W))
claim_d = Decimal(c_num) / 1000
gate_d = Decimal(GATE.numerator) / Decimal(GATE.denominator)
add("claim_window_row",
    log2w8 + Decimal("0.00002") <= claim_d and claim_d <= gate_d - Decimal("0.008"),
    {"log2_W_8dp": str(log2w8), "claim": str(claim_d), "gate": str(gate_d),
     "claim_minus_log2W": str(claim_d - log2w8), "gate_minus_claim": str(gate_d - claim_d),
     "window": "log2 W + 2e-5 <= claim <= gate - 0.008"})

# ---- 7: success, recomputed and re-rendered exactly -----------------------
rho = ledger.rho_ideal(ledger.A_H1)
recomputed = rho - sum(v for _k, v in ledger.LOSSES)
add("success_recomputed_from_terms", recomputed == ledger.SUCCESS,
    {"rho": gen_tables.dec10(rho), "losses": len(ledger.LOSSES),
     "success": gen_tables.dec10(recomputed)})
add("success_above_floor", ledger.SUCCESS >= ledger.SUCCESS_FLOOR,
    {"success": gen_tables.dec10(ledger.SUCCESS), "floor": gen_tables.fmt_exact(ledger.SUCCESS_FLOOR)})
add("success_printed_exactly", gen_tables.dec10(ledger.SUCCESS) in block,
    {"printed": gen_tables.dec10(ledger.SUCCESS)})
add("each_loss_smaller_than_margin",
    all(v < ledger.SUCCESS - ledger.SUCCESS_FLOOR for _k, v in ledger.LOSSES),
    {"margin": gen_tables.dec10(ledger.SUCCESS - ledger.SUCCESS_FLOOR)})

# ---- 8: the toy record, read raw ------------------------------------------
with open(os.path.join(PKG, "evidence", "toy.json")) as f:
    toy = json.load(f)
cr = toy["collision_report"]
per_e_ok = all(v["e_ok"] == v["trials"] and v["x_ok"] == v["trials"]
               for v in cr["per_e"].values()) and len(cr["per_e"]) == 6
add("toy_verdict_pass", toy["verdict"] == "PASS", {"verdict": toy["verdict"]})
add("toy_e_resolution_zero_failures",
    cr["failed"] == 0 and cr["e_resolved"] == cr["trials"] and cr["x_recovered"] == cr["trials"]
    and per_e_ok,
    {"failed": cr["failed"], "e_resolved": cr["e_resolved"], "trials": cr["trials"],
     "x_recovered": cr["x_recovered"], "elements": sorted(cr["per_e"]),
     "per_element_all_ok": per_e_ok})
add("toy_invariant_every_transition",
    toy["invariant_checks"] == toy["invariant_pass_steps"] and toy["invariant_checks"] > 1000,
    {"checks": toy["invariant_checks"], "steps": toy["invariant_pass_steps"]})
sw = toy["sweep"]["rows"]
sel_total = sum(r["selector_checks"] for r in sw)
add("toy_selector_checks", sel_total > 0 and toy["selector"] == "scripts/selector.py",
    {"selector_checks": sel_total, "rows": len(sw)})
dp0 = [r for r in sw if r["dp_bits"] == 0][0]
rec = dp0["recovery_rate"]
add("measured_recovery_at_charged_budget",
    Decimal(str(rec["kappa_1.005"])) >= Decimal("0.39"),
    {"kappa_1.005": rec["kappa_1.005"], "kappa_1.000": rec["kappa_1.000"],
     "kappa_1.025": rec["kappa_1.025"], "trials": dp0["trials"],
     "charged_budget": "A_H1 = " + gen_tables.fmt_exact(ledger.A_H1)})
add("A_measured_consistent_with_ledger",
    Decimal(str(toy["A_measured"])) <= Decimal(gen_tables.fmt_exact(ledger.A_H1)),
    {"A_measured": toy["A_measured"], "A_H1": gen_tables.fmt_exact(ledger.A_H1)})

# ---- 9: packaging ----------------------------------------------------------
with open(os.path.join(PKG, "dependencies.json")) as f:
    deps = json.load(f)
with open(os.path.join(PKG, "submission.json")) as f:
    sub = json.load(f)
add("dependencies_empty", deps == [], {"dependencies.json": deps})
add("experiments_empty", sub.get("experiments") == [], {"experiments": sub.get("experiments")})
add("claim_literal_matches_ledger",
    frac(sub["claims"][0]["score"]["numerator"]) / frac(sub["claims"][0]["score"]["denominator"])
    == claim,
    {"submission": sub["claims"][0]["score"], "ledger": gen_tables.fmt_exact(claim)})
ev_bytes = 0
for root, _dirs, files in os.walk(PKG):
    if "__pycache__" in root:
        continue
    for fn in files:
        ev_bytes += os.path.getsize(os.path.join(root, fn))
add("evidence_budget", ev_bytes <= EVIDENCE_BUDGET,
    {"bytes": ev_bytes, "limit": EVIDENCE_BUDGET})

# ---- verdict ---------------------------------------------------------------
passed = sum(1 for c in CHECKS if c["pass"])
cert = {"package": "plain-sixclass-batchamort-40pct-v17",
        "checks": CHECKS, "passed": passed, "total": len(CHECKS),
        "verdict": "PASS" if passed == len(CHECKS) else "FAIL",
        "claim": gen_tables.fmt_exact(claim), "log2_W": gen_tables.log2_8(ledger.W),
        "success": gen_tables.dec10(ledger.SUCCESS)}
with open(os.path.join(PKG, "certificate.json"), "w") as f:
    json.dump(cert, f, indent=1)
bad_names = [c["check"] for c in CHECKS if not c["pass"]]
print(json.dumps({"verdict": cert["verdict"], "passed": passed, "total": len(CHECKS),
                  "failed": bad_names, "claim": cert["claim"], "log2_W": cert["log2_W"]}))
sys.exit(0 if cert["verdict"] == "PASS" else 1)
