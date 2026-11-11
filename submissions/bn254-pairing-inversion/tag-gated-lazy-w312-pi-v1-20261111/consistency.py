"""Part C: cross-document parameter consistency (standard library only).

Extracts every numeric parameter the claim depends on from submission.json, proof.md, machine.md, the
code (lazy.py, experiment.py, reference.py) and the ledger output of the same run, and asserts that they
agree.  Rounded figures stated in the text are checked in the safe direction against the exact ledger
value: stated upper bounds must be >= the ledger value, stated lower bounds <= it.  The result lists,
for each parameter, every value found in every source, so a disagreement is visible, not just flagged.
"""
import json
import re
from fractions import Fraction


def _read(base, name):
    return (base / name).read_text()


def _ints(pattern, text):
    return sorted({int(m) for m in re.findall(pattern, text)})


def run(base, ledger, part_l, part_p, base_out):
    sub = json.loads(_read(base, "submission.json"))
    stmt = sub["claims"][0]["statement"]
    score = Fraction(int(sub["claims"][0]["score"]["numerator"]), int(sub["claims"][0]["score"]["denominator"]))
    proof, machine = _read(base, "proof.md"), _read(base, "machine.md")
    lazy, exper, refpy = _read(base, "lazy.py"), _read(base, "experiment.py"), _read(base, "reference.py")
    docs = {"statement": stmt, "proof.md": proof, "machine.md": machine}
    checks = {}

    def same(name, found):
        vals = sorted({v for vs in found.values() for v in vs})
        checks[name] = {"values_by_source": found, "ok": len(vals) == 1 and all(found.values())}

    # 1. distinguished-point predicate: one bit count everywhere
    dp = {}
    for k, t in docs.items():
        dp[k] = _ints(r"t\(Y\)\s*\*\s*2\^(\d+)\s*<\s*p", t) + _ints(r"ceil\(p/2\^(\d+)\)", t) \
            + _ints(r"2 E\[T\]/2\^(\d+) \+ 1", t)
    dp["lazy.py ledger E_walks"] = _ints(r"E_walks = E_steps / Decimal\(2\) \*\* (\d+)", lazy)
    dp["lazy.py DP_BITS_PRODUCTION"] = _ints(r"DP_BITS_PRODUCTION = (\d+)", lazy)
    dp["experiment.py DP_BITS"] = _ints(r"\nDP_BITS = (\d+)", exper)
    dp["reference.py production dbits"] = _ints(r"production: cbits=\d+, dbits=(\d+)", refpy)
    if part_p.get("dbits") is not None:
        dp["part P executed dbits"] = [part_p.get("dbits")]
    dp = {k: sorted(set(v)) for k, v in dp.items() if k != "statement" or v}
    same("distinguished_point_bits", dp)
    checks["no_other_dp_definition"] = {
        "ok": not re.search(r"2\^80\b|<< 80\)|dbits=80|2\^-80\b", proof + machine + stmt + exper + refpy)}

    # 2. caps and counter widths
    same("step_cap_log2", {"proof.md": _ints(r"S_cap = 2\^(\d+)", proof),
                           "statement": _ints(r"2\^(\d+)-step global", stmt),
                           "ledger": [ledger["log2_step_cap"]]})
    same("walk_cap_log2", {"proof.md": _ints(r"at most 2\^(\d+) walks", proof),
                           "statement": _ints(r"at most 2\^(\d+) walks", stmt),
                           "ledger": _ints(r"walk_cap_2\^(\d+)_markov", json.dumps(ledger))})
    same("walk_length_log2", {"proof.md": _ints(r"each at most 2\^(\d+) steps", proof),
                              "statement": _ints(r"walks of at most 2\^(\d+) steps", stmt),
                              "machine.md": _ints(r"per-walk length counter \(initially 2\^(\d+)\)", machine)})
    same("exact_meter_log2", {"proof.md": _ints(r"at most 2\^(\d+) exact evaluations", proof),
                              "statement": _ints(r"at most 2\^(\d+) exact evaluations", stmt),
                              "ledger": _ints(r"exact_meter_2\^(\d+)_markov", json.dumps(ledger))})
    same("probe_meter_log2", {"proof.md": _ints(r"2\^(\d+) store probes", proof),
                              "statement": _ints(r"at most 2\^(\d+) store probes", stmt),
                              "ledger": _ints(r"store_probe_meter_2\^(\d+)_markov", json.dumps(ledger))})
    same("draws_per_scalar", {"proof.md": [int(v.replace(",", "")) for v in re.findall(r"at most ([\d,]+) draws per scalar", proof)],
                              "statement": [int(v.replace(",", "")) for v in re.findall(r"at most ([\d,]+) draws per random scalar", stmt)]})
    widths = _ints(r"three-word \((\d+)-bit\)", proof)
    widths_table = {m.group(1).strip(): (m.group(2), m.group(3)) for m in
                    re.finditer(r"\| ([a-z -]+?) \| ([0-9^,]+) \| (one|two|three) words?", machine)}
    wmap = {"one": 1, "two": 2, "three": 3}
    fits = {}
    for name, (cap, w) in widths_table.items():
        capv = 2 ** int(cap[2:]) if cap.startswith("2^") else int(cap.replace(",", ""))
        fits[name] = capv < 2 ** (64 * wmap[w])
    checks["every_counter_width_fits_its_cap"] = {"counters": fits,
        "ok": len(fits) >= 6 and all(fits.values()) and "exact-evaluation meter" in fits and "store-probe meter" in fits}
    checks["global_counter_fits"] = {"counter_bits": widths, "step_cap_log2": ledger["log2_step_cap"],
                                     "ok": bool(widths) and min(widths) > ledger["log2_step_cap"] and "three-word" in machine}

    # 3. instruction costs
    li = part_l.get("executed_instructions") or [None]
    same("window_product_instr", {"proof.md": [int(v.replace(",", "")) for v in re.findall(r"product\**\s+of ([\d,]+) executed", proof)],
                                  "statement": [int(v.replace(",", "")) for v in re.findall(r"Fp12 product of ([\d,]+) scalar", stmt)],
                                  "part L executed": [li[0]]})
    same("window_end_charged", {"proof.md": [int(v.replace(",", "")) for v in re.findall(r"at most ([\d,]+) instead of", proof)],
                                "statement": [int(v.replace(",", "")) for v in re.findall(r"window end costs at most ([\d,]+)", stmt)],
                                "ledger": [ledger["window_end_charged"]]})
    per_step = Fraction(ledger["per_step_instr"])
    same("per_step_instr", {"proof.md": [Fraction(v) for v in re.findall(r"= \*\*([\d,]+/[\d,]+)\*\*", proof)],
                            "statement": [Fraction(v) for v in re.findall(r"a step 273 \+ 9216/312 = ([\d,]+/[\d,]+)", stmt)],
                            "ledger": [per_step]})
    same("gate_K", {"proof.md": _ints(r"K = (\d{10,})", proof), "statement": _ints(r"K = (\d{10,})", stmt),
                    "ledger": [int(ledger["gate_K"])]})

    # 4. claim
    same("claim", {"submission.json score": [score], "ledger CLAIM": [Fraction(ledger["claim"])],
                   "proof.md": [Fraction(v) for v in re.findall(r"< 2\^(\d+\.\d+)\n", proof)][:1],
                   "statement": [Fraction(v) for v in re.findall(r"< 2\^(\d+\.\d+) reference units", stmt)]})
    checks["claim_covers_ledger"] = {"ok": bool(ledger["claim_ok"])}

    # 5. rounded figures in the safe direction
    def stated(rx, texts):
        out = []
        for t in texts:
            for m in re.findall(rx, t):
                for v in (m if isinstance(m, tuple) else (m,)):
                    if v:
                        out.append(float(v))
        return out
    bounds = [
        ("step_cap_markov <=", stated(r"step cap at most (0\.\d+)|<= E\[T\]/2\^133 <= (0\.\d+)", [stmt, proof]), ledger["losses"]["step_cap_markov"], "upper"),
        ("walk_cap_markov <=", stated(r"walk cap at most (0\.\d+)|<= E\[walks\]/2\^52 <= (0\.\d+)", [stmt, proof]), ledger["losses"]["walk_cap_2^52_markov"], "upper"),
        ("failure budget below", stated(r"failure budget below (0\.\d+)|sum is below (0\.\d+)", [stmt, proof]), ledger["failure_budget_past_gate"], "upper"),
        ("success at least", stated(r"probability at least (0\.9\d+)|success past the gate is at least (0\.9\d+)", [stmt, proof]), ledger["success_past_gate_lower"], "lower"),
        ("every-run log2 <=", stated(r"capped at 2\^(\d+\.\d+) units", [stmt, proof]), ledger["log2_every_run_bound_units"], "upper"),
        ("charged mean steps log2 <=", stated(r"<= 2\^(\d+\.\d+) steps", [stmt, proof]), ledger["log2_charged_mean_steps"], "upper"),
    ]
    import math
    bounds += [
        ("expected exact evaluations log2 <=", stated(r"number is E\[T\] x 97/2\^61 below 2\^(\d+\.\d+)", [proof]), ledger["log2_expected_exact_evaluations"], "upper"),
        ("exact meter Markov log2 <=", stated(r"P\(> 2\^90\) < 2\^(-\d+\.\d+)", [proof]), math.log2(ledger["losses"]["exact_meter_2^90_markov"]), "upper"),
        ("walk-length/self-cycle union log2 <=", stated(r"union over at most 2\^52 walks:\s+below 2\^(-\d+)", [proof]), ledger["log2_walk_cap_union_bound"], "upper"),
        ("ledger walk-length loss covers union bound", [math.log2(ledger["losses"]["walk_length_cap_and_self_cycles"])], ledger["log2_walk_cap_union_bound"], "upper"),
    ]
    p_val = int(json.loads(_read(base, "target-profile.json"))["field_modulus"], 16)
    excess = stated(r"below 2\^-84 \+ 2\^(-\d+);", [proof])
    checks["density_excess_bound"] = {"stated_log2": excess, "exact_log2_one_over_p": -math.log2(p_val),
                                      "ok": bool(excess) and all(v >= -math.log2(p_val) for v in excess)}
    for name, vals, exact, kind in bounds:
        vals = [v for v in vals if v]
        ok = bool(vals) and all((v >= exact) if kind == "upper" else (v <= exact) for v in vals)
        checks["rounded:" + name] = {"stated": vals, "exact": exact, "direction": kind, "ok": ok}

    # 6. predicate actually executed at full width
    checks["full_width_predicate_walks"] = {
        "part_P": part_p, "base_dp_bits": dp.get("experiment.py DP_BITS"),
        "ok": base_out is not None and dp.get("experiment.py DP_BITS") == [84]
              and (not part_p or (part_p.get("dbits") == 84 and part_p.get("dp_flag_checks", 0) > 0 and "error" not in part_p))}
    return {"checks": checks, "all_ok": all(c["ok"] for c in checks.values()),
            "failed": [k for k, c in checks.items() if not c["ok"]]}
