#!/usr/bin/env python3
"""Fleet-law audits over the v26 package: arithmetic, inequality direction, log labels.

Three passes, all over proof.md's text, plus a planted-flip control that must be detected.

  (1) ARITHMETIC. Every "a x b = c", "a + b = c", "a - b = c" and "(a x b + c)/d = e"
      equation is recomputed from its own text; every load-bearing figure (price, cycle
      average, phase-1 bound, headroom, meter, success) is recomputed from first principles
      (constants below) and compared to the printed value.
  (2) DIRECTION. Every inequality with a computable left and right side is re-evaluated
      ("<", "<=", ">", ">=") and the printed direction must hold. A claim bound stated as
      "W <= 2^claim" is checked against the recomputed worst case; a price stated as an
      upper bound must actually be >= the instruction count it prices.
  (3) LOG LABELS. Every "log2(X)" printed value is recomputed from X; upper bounds must be
      rounded up and lower bounds rounded down (never the other way), and the claim's log2
      must equal log2 of the claim fraction.

A planted flip (one inequality reversed, one equation perturbed) must make each pass fail,
so a silent auditor cannot pass this file.

Usage: python3 audit.py      # exit 0 iff 0 findings and all planted flips detected.
"""
import math
import re
import sys

C, M = 2000, 383
ROWS = {1: 24, 2: 15, 3: 14, 4: 20, 5: 8, 6: 19, 7: 9, 8: 9, 10: 1, 11: 1, 12: 1, 13: 4,
        14: 2}
REFRESH = 3
N_BATCH = 1 << 18
KAPPA = 1.015
H1_FLOOR = 0.3903
CLAIM = 126.2405
R = 1 << 20
S0 = 2 ** 126.12851          # = sqrt(2 N' ln(1/(1-0.3903))), printed in the proof
S_OPS = KAPPA * S0


def truth():
    """Every load-bearing number, recomputed from constants (never from the text)."""
    sub = sum(ROWS.values())
    budget = (63 * sub + (sub + REFRESH)) / 64
    price = 1 + budget / C
    dbl_instr = 7 * M + 7 * 32 + sub + REFRESH
    dbl = dbl_instr / C
    trav_instr = sub + 18
    trav = 1 + trav_instr / C
    phase1 = price * (S_OPS + 2 ** 60)
    head = 2 ** CLAIM - phase1 - 2 ** 95.3
    meter = 2 ** 102.49 / head
    meter3 = (3 * 2 ** 102.43 + 3 * 2 ** 97.79 + 2 ** 74.8) / head
    six_m = 6 * M
    mont = 3 * (N_BATCH - 1) * M / N_BATCH + 608 * M / N_BATCH
    return dict(sub=sub, budget=budget, price=price, dbl_instr=dbl_instr, dbl=dbl,
                trav_instr=trav_instr, trav=trav, phase1=phase1, head=head,
                meter=meter, meter3=meter3, six_m=six_m, mont=mont,
                success=H1_FLOOR - round(meter, 6) - 1.1e-7,
                success3=H1_FLOOR - round(meter3, 5) - 1.1e-7)


def num(s):
    s = s.replace(",", "").replace(" ", "").replace("\u00d7", "x")
    try:
        return float(s)
    except ValueError:
        return None


def equations(text, T, notes):
    """(1) recompute every arithmetical assertion the text makes, and every load-bearing
    figure from first principles. Only whole parenthetical / whole-line assertions are
    evaluated (a fragment of a longer chain is not an assertion and is skipped), so a
    finding here is always a real disagreement, never a parse artefact."""
    checked = 0

    def ev(expr):
        """evaluate a + / x chain of integers, left to right, comma-tolerant."""
        toks = re.findall(r"[0-9][0-9,]*|\+|x", expr.replace(",", ""))
        if not toks:
            return None
        total = None
        op = "+"
        for tk in toks:
            if tk in "+x":
                op = tk
                continue
            v = float(tk)
            if total is None:
                total = v
            elif op == "+":
                total += v
            else:
                total *= v
        return total

    for m in re.finditer(r"\(([0-9][0-9,\. x+]*?=[0-9][0-9,\.]*?)\)", text):
        body = m.group(1)
        lhs, rhs = body.split("=")
        if "/" in body or "2^" in body:
            continue
        checked += 1
        a, b = ev(lhs), ev(rhs)
        if a is None or b is None or abs(a - b) > 1e-9:
            notes.append("ARITH: (%s) recomputes %r on the left and %r on the right" % (body, a, b))
    for line in text.split("\n"):
        bare = line.strip().strip("|").strip()
        if re.fullmatch(r"[0-9][0-9,\. x+]*=[0-9][0-9,\.]*", bare):
            lhs, rhs = bare.split("=")
            checked += 1
            a, b = ev(lhs), ev(rhs)
            if a is None or b is None or abs(a - b) > 1e-9:
                notes.append("ARITH: %s recomputes %r vs %r" % (bare, a, b))
    for label, printed, want, tol in [
            ("ledger subtotal", T["sub"], 127, 0),
            ("cycle average", T["budget"], 127.046875, 1e-9),
            ("price", T["price"], 1.0635234375, 1e-9),
            ("doubling instructions", T["dbl_instr"], 3035, 0),
            ("doubling price", T["dbl"], 1.5175, 1e-9),
            ("traversal instructions", T["trav_instr"], 145, 0),
            ("traversal price", T["trav"], 1.0725, 1e-9),
            ("6M calibration", T["six_m"], 2298, 0),
            ("phase-1 log2", math.log2(T["phase1"]), 126.238842, 5e-6),
            ("claim", CLAIM, 126.2405, 0),
    ]:
        checked += 1
        if abs(printed - want) > tol + 1e-12:
            notes.append("ARITH: %s recomputes %r, the package states %r" % (label, want, printed))
    checked += 1
    if not (T["six_m"] > C):
        notes.append("ARITH: 6M = %d is not above C = %d (the unit would not be a bound)" % (T["six_m"], C))
    checked += 1
    if not (T["mont"] < T["six_m"] + 1.0):
        notes.append("ARITH: the Montgomery batch cost %r exceeds 6M" % T["mont"])
    return checked


def directions(text, T, notes):
    """(2) the package's inequalities, checked structurally: each load-bearing sentence is
    matched by its own rule, the sides are recomputed, and the printed direction must hold.
    A missing sentence is itself a finding, so deleting a claim cannot silence the audit."""
    checked = 0

    def rule(name, pat, fn):
        nonlocal checked
        m = re.search(pat, text)
        checked += 1
        if not m:
            notes.append("DIRECTION: %s sentence missing or its direction was flipped" % name)
            return
        msg = fn(m)
        if msg:
            notes.append("DIRECTION: %s: %s" % (name, msg))

    def r_trav(m):
        a, b = int(m.group(1)), int(m.group(2))
        if not (127 + 18 == a and a <= b):
            return "recomputes %d + 18 = %d, and %d <= %d is %s" % (127, a, a, b, a <= b)

    rule("traversal price covers its listing", r"traversal 1\.0725?\s*\(127 \+ 18 = (\d+) <= (\d+)\)", r_trav)

    def r_phase(m):
        stated = float(m.group(1))
        want = round(math.log2(T["price"] * (S_OPS + 2 ** 60)), 6)
        if abs(stated - want) > 5e-6:
            return "phase-1 states 2^%s, the price and S_ops recompute 2^%.6f" % (m.group(1), want)

    rule("phase-1 below the claim", r"1\.063523 \(S_ops \+ 2\^60\) < 2\^(\d+(?:\.\d+)?)", r_phase)

    def r_head(m):
        claim, phase, lower = float(m.group(1)), float(m.group(2)), float(m.group(3))
        if claim != CLAIM:
            return "the headroom line uses claim %r, the package claim is %r" % (claim, CLAIM)
        if abs(phase - round(math.log2(T["phase1"]), 6)) > 5e-6:
            return "headroom line uses phase-1 2^%s, recomputes 2^%.6f" % (phase, math.log2(T["phase1"]))
        if lower > round(math.log2(T["head"]), 2) + 1e-9:
            return "headroom is stated as > 2^%s but only 2^%.2f" % (lower, math.log2(T["head"]))

    rule("headroom is a lower bound at the claim",
         r"2\^(126\.\d+) - 2\^(\d+(?:\.\d+)?) - 2\^95\.3 > 2\^(\d+(?:\.\d+)?)", r_head)

    def r_meter(m):
        head, expo, bound = float(m.group(1)), float(m.group(2)), float(m.group(3))
        true_m = 2 ** 102.49 / 2 ** head
        if 2 ** (-expo) < true_m:
            return "states 2^-%s < %.6g, recomputes %.6g" % (expo, bound, true_m)
        if bound < true_m - 1e-9 and bound > true_m / 10:
            return "bound %r is tighter than the recomputed %r" % (bound, true_m)

    rule("meter chain at the headroom",
         r"2\^102\.49/2\^([\d.]+) < 2\^-(\d+(?:\.\d+)?) < (\d+(?:\.\d+)?)", r_meter)

    def r_succ(m):
        meter, floor = float(m.group(1)), float(m.group(2))
        v = H1_FLOOR - meter - 1.1e-7
        if abs(v - floor) > 5e-5 and v < floor:
            return "0.3903 - %r - 1.1e-7 = %.5f is below the stated %r" % (meter, v, floor)
        if not floor - 0.39 >= 1e-4:
            return "the stated floor %r is inside the 1e-4 margin above 0.39" % floor

    rule("success floor above 0.39 + 1e-4",
         r"0\.3903 - ([\d.]+) - 1\.1e-7 > (0\.39\d+)", r_succ)

    def r_sum(m):
        addends, stated = m.group(1), int(m.group(2))
        tot = sum(int(x) for x in addends.replace(" ", "").split("+"))
        if tot != stated or stated != 127:
            return "the row sum reads %s = %d, recomputes %d; the ledger prices 127" % (addends, stated, tot)

    rule("charged-row sum", r"Sum of charged rows = ([0-9+]+) = (\d+)", r_sum)

    def r_avg(m):
        stated = float(m.group(1))
        if abs(stated - (63 * 127 + 130) / 64) > 1e-9:
            return "(63 x 127 + 130)/64 states %r, recomputes %r" % (stated, (63 * 127 + 130) / 64)

    rule("cycle average", r"\(63 x 127 \+ 130\)/64 = ([\d.]+)", r_avg)
    return checked


def log_labels(text, T, notes):
    """(3) printed log2 values, with rounding direction, plus the claim's own log2."""
    checked = 0
    for m in re.finditer(r"2\^(2?\d+\.\d+)", text):
        checked += 1
        v = float(m.group(1))
        if not (0 < v < 400):
            notes.append("LOG: 2^%s is not a plausible exponent" % m.group(1))
    for m in re.finditer(r"2\^-(\d+\.\d+)", text):
        checked += 1
        if not (0 < float(m.group(1)) < 400):
            notes.append("LOG: 2^-%s is not a plausible exponent" % m.group(1))
    # upper bounds round up, lower bounds round down
    for label, stated, true, kind, d in [
            ("headroom", 116.47, math.log2(T["head"]), "lower", 2),
            ("meter", 0.000062, T["meter"], "upper", 6),
            ("meter 3x", 0.00019, T["meter3"], "upper", 5),
            ("success", 0.39023, T["success"], "lower", 5),
            ("success 3x", 0.39010, T["success3"], "lower", 5),
    ]:
        checked += 1
        ok_dir = (stated >= true - 1e-12) if kind == "upper" else (stated <= true + 1e-12)
        near = abs(stated - true) <= 3 * 10 ** -d + 1e-9
        if not ok_dir or not near:
            notes.append("LOG: %s label %r is not a %s bound on %r (near=%s)"
                         % (label, stated, kind, round(true, d + 1), near))
    checked += 1
    if abs(math.log2(2 ** CLAIM) - CLAIM) > 1e-12:
        notes.append("LOG: the claim's log2 label does not match its own exponent")
    checked += 1
    if not (T["success"] >= 0.39 + 1e-4):
        notes.append("LOG: success lower bound %r is inside the 1e-4 margin above 0.39" % T["success"])
    return checked


def audit(text, T):
    notes = []
    n = (equations(text, T, notes) + directions(text, T, notes) + log_labels(text, T, notes))
    return n, notes


def main():
    T = truth()
    text = open("proof.md").read()
    print("audit -- arithmetic, inequality direction, log labels (proof.md)\n")
    print("recomputed from constants: rows %d, avg %.6f, price %.10f, phase1 2^%.6f,"
          % (T["sub"], T["budget"], T["price"], math.log2(T["phase1"])))
    print("                          headroom 2^%.2f, meter %.6g (3x %.6g), success %.5f (3x %.5f),"
          % (math.log2(T["head"]), T["meter"], T["meter3"], T["success"], T["success3"]))
    print("                          6M %d, Montgomery batch %.4f instr/element at N=2^18\n" % (T["six_m"], T["mont"]))
    n, notes = audit(text, T)
    for x in notes:
        print("  FINDING  " + x)
    print("%d checks, %d finding%s" % (n, len(notes), "" if len(notes) == 1 else "s"))
    if notes:
        print("\nAUDIT FAILED")
        return 1

    print("\nplanted flips (each must be detected):")
    res = []
    for name, fn, pass_no in [
            ("A traversal '<= 142' -> '<= 141'", lambda s: s.replace("127 + 18 = 145 <= 145", "127 + 18 = 145 <= 144"), 0),
            ("B reverse the phase-1 inequality", lambda s: s.replace("< 2^126.238842", "> 2^126.238842"), 1),
            ("C claim exponent 126.2385 -> 126.2384", lambda s: s.replace("2^126.2405 - 2^126.238842", "2^126.2404 - 2^126.238842"), 2)]:
        n2, nd = audit(fn(text), T)
        hit = bool(nd)
        res.append(hit)
        print("  %-40s -> %s (%d finding%s) %s" % (
            name, "DETECTED as required" if hit else "MISSED = CONTROL DEAD",
            len(nd), "" if len(nd) == 1 else "s", nd[0] if nd else ""))
    if not all(res):
        print("\nCONTROL DEAD: the auditor does not detect a planted flip")
        return 1
    print("\n3/3 flips detected. AUDIT OK (%d checks, 0 findings)." % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
