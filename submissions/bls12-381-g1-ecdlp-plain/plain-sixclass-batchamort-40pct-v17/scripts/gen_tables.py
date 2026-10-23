"""gen_tables.py -- renders the generated block of argument.md from ledger.py.

Nothing here is typed by hand: every cell is an attribute of ledger.py or a
function of one.  verify_certificate.py re-renders this module and byte-compares
it against the block in argument.md, then re-evaluates every `F:` line.
"""
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import ledger  # noqa: E402

BEGIN = "<!-- BEGIN GENERATED: rendered by scripts/gen_tables.py from ledger.py; do not edit -->"
END = "<!-- END GENERATED -->"


def fmt_exact(fr):
    """Exact decimal if the denominator is 2^a*5^b, else the exact rational."""
    fr = F(fr)
    n, d = fr.numerator, fr.denominator
    d2, a, b = d, 0, 0
    while d2 % 2 == 0:
        d2 //= 2
        a += 1
    while d2 % 5 == 0:
        d2 //= 5
        b += 1
    if d2 != 1:
        return "%d/%d" % (n, d)
    k = max(a, b)
    scaled = n * 10 ** k // d
    s = str(scaled)
    if k == 0:
        return s
    s = s.rjust(k + 1, "0")
    return s[:-k] + "." + s[-k:]


def dec10(fr):
    """Exact rounding of a Fraction to 10 decimal places, as a string."""
    fr = F(fr)
    q = fr.numerator * 10 ** 10
    neg = q < 0
    q = abs(q)
    scaled = (2 * q + fr.denominator) // (2 * fr.denominator)   # round half up
    s = str(scaled).rjust(11, "0")
    out = s[:-10] + "." + s[-10:]
    return ("-" if neg else "") + out


def log2_8(fr):
    """8-decimal rendering of log2, from integers only (no float compare)."""
    import math
    fr = F(fr)
    v = math.log2(fr.numerator) - math.log2(fr.denominator)
    return "%.8f" % v


def main(argv):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)
    with open(os.path.join(os.path.dirname(HERE), "argument.md")) as f:
        doc = f.read()
    block = doc[doc.index(BEGIN):doc.index(END) + len(END)]
    ok = block.rstrip("\n") == render().rstrip("\n")
    print("generated block matches ledger.py" if ok else "BLOCK MISMATCH")
    return 0 if ok else 1


def render():
    L = ledger
    out = []
    w = out.append
    w(BEGIN)
    w("")
    w("## 7. Step ledger (every row is charged)")
    w("")
    w("One step is one batched incomplete affine G1 addition (1.0 reference unit =")
    w("%d instructions) plus the control work below. `insn` counts scalar AArch64" % L.C)
    w("instructions; 6-limb values are 381-bit `Fp` elements.")
    w("")
    w("| row | operations (count each) | insn |")
    w("|---|---|---:|")
    for name, ops, n in L.ROWS:
        w("| %s | %s | %d |" % (name, ops, n))
    w("| **sum** | ten rows | **%d** |" % L.STEP_INSN)
    w("")
    w("## 8. Budget, constants, losses")
    w("")
    w("| quantity | rule | value |")
    w("|---|---|---:|")
    w("| S_ops | ceil(A_H1 * S0) | %s |" % fmt_exact(L.S_ops))
    w("| phase1 | step_price * (S_ops + 2^60) | %s |" % log2_8(L.PHASE1))
    w("| batch_inversion_total | inv_insn * (n_batches + 1) / %d | %s |"
      % (L.C, log2_8(L.batch_inversion_total)))
    w("| extras | table+setup+meter+routing+verify+starts+drain+report | %s |" % log2_8(L.EXTRAS))
    w("| W | phase1 + extras + batch_inversion_total | %s |" % log2_8(L.W))
    w("| claim | ceil_3dp(log2 W) | %s |" % fmt_exact(L.CLAIM))
    w("")
    w("| constant | value |")
    w("|---|---:|")
    w("| log2 N | %s |" % log2_8(L.N))
    w("| log2 S0 | %s |" % log2_8(L.S0))
    w("| reachable columns R_eff | %d |" % L.R_EFF)
    w("| A_MIN (0.001 grid) | %s |" % fmt_exact(L.A_MIN))
    w("| A_H1 = A_MIN + 0.005 | %s |" % fmt_exact(L.A_H1))
    w("| log2 S_ops | %s |" % log2_8(L.S_ops))
    w("| step instructions | %d |" % L.STEP_INSN)
    w("| step price | %s |" % fmt_exact(L.STEP_PRICE))
    w("| inversion multiplies (254 sq + 11 mul) | %d |" % L.INV_MULTS)
    w("| inversion instructions | %d |" % L.INV_INSN)
    w("| inversion instructions per step at B = 2^18 | %s |" % fmt_exact(L.INV_PER_STEP))
    w("| full batches over S_ops + 2^60 | %d |" % L.n_batches)
    w("| doubling instructions | %d |" % L.dbl_instructions)
    w("| doubling price | %s |" % fmt_exact(L.DOUBLING_PRICE))
    w("| rho_ideal(A_H1) | %s |" % dec10(L.rho_ideal(L.A_H1)))
    w("| success | %s |" % dec10(L.SUCCESS))
    w("| log2 peak memory | %s |" % log2_8(L.MEMORY_BYTES))
    w("| records cap (mean + 4 sigma + workers) | %d |" % L.records)
    w("")
    w("| loss term | log2 | charged value |")
    w("|---|---:|---:|")
    for name, v in L.LOSSES:
        w("| %s | %s | %s |" % (name, log2_8(v), fmt_exact(v)))
    w("")
    w("### 8.1 Formula lines (every line is evaluated by scripts/verify_certificate.py)")
    w("")
    w("Each line's left side is re-evaluated in exact rational arithmetic and must")
    w("equal the printed right side exactly; the `=` is a claim, not decoration.")
    w("")
    w("F: %s = %s    <!-- step_sum -->"
      % (" + ".join(str(n) for _a, _b, n in L.ROWS), fmt_exact(L.STEP_INSN)))
    w("F: 1 + %d / %d = %s    <!-- step_price -->"
      % (L.STEP_INSN, L.C, fmt_exact(L.STEP_PRICE)))
    w("F: 254 + 11 = %d    <!-- inversion_multiplies -->" % L.INV_MULTS)
    w("F: %d * %d = %d    <!-- inversion_instructions -->" % (L.INV_MULTS, L.M_insn, L.INV_INSN))
    w("F: %d / %d = %s    <!-- inversion_per_step -->"
      % (L.INV_INSN, L.BATCH, fmt_exact(L.INV_PER_STEP)))
    w("F: 4 * %d + 6 * %d + 3 * 48 + 56 = %d    <!-- doubling_instructions -->"
      % (L.M_insn, L.A_insn, L.dbl_instructions))
    w("F: %d / %d = %s    <!-- doubling_price -->"
      % (L.dbl_instructions, L.C, fmt_exact(L.DOUBLING_PRICE)))
    w("F: %d * %d = %d    <!-- batch_inversion_instructions -->"
      % (L.n_batches + 1, L.INV_INSN, L.INV_INSN * (L.n_batches + 1)))
    w("F: %d * %d * 192 + 6 * %d * 8 + 2 ** 40 = %d    <!-- memory_bytes -->"
      % (L.nodes, L.node_cap_records, L.R_TAB, L.MEMORY_BYTES))
    w("")
    w("<!-- END GENERATED -->")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    import sys as _sys
    _sys.exit(main(_sys.argv[1:]))
