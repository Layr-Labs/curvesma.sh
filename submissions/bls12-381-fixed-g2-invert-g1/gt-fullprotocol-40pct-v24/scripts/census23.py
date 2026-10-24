#!/usr/bin/env python3
"""census22.py -- EXACT pattern census for the Frobenius-class GT walk.

WHY THIS EXISTS (v20 verdict flaw 2).  v20 claimed longer balanced returns were
negligible because "at least one adjacent slot change is forced, so the product is a
product of independent group elements".  That is FALSE.  GT is abelian, so

    M_a^x M_b^y M_a^-x M_b^-y = 1   IDENTICALLY,

and the inverse of the factor applied at slot a is available AT THE SAME SLOT because
p^6 = -1 mod r: a revisit of slot a whose kappa label is shifted by 6 applies exactly
the inverse of the earlier factor.  A slot word whose per-slot exponent sums vanish is
therefore an IDENTICAL product, not a coincidence.

EXACT CRITERION.  zeta = p^{-1} mod r is a primitive 12th root of unity mod r.  Over a
slot word the total exponent at slot s is

    E_s = sum_{i: s_i = s} zeta^{kappa_i}   (an integer combination of 12th roots)

and the product is identically 1 iff E_s = 0 for every slot s.  zeta has minimal
polynomial Phi_12(x) = x^4 - x^2 + 1, so E_s = 0 iff the coefficient vector of
sum_i x^{kappa_i} reduced mod Phi_12 is the zero vector -- exact integer linear
algebra, no probability and no independence assumption.

WHY IT MATTERS.  If the product is 1 the walk returns to the SAME group element y, so
the class sequence becomes periodic and the trail can never report again (or reports
one repeated class).  A trap is therefore not a rare coincidence to be dismissed; it
is a structural property of the slot word, and it must be excluded by the guard or
priced as wasted work.

RATE.  Declared premise H4: for a word of length n with m distinct slots the
slot-recurrence constraints carry probability R^{-(n-m)}; the kappa labels are
uniform on Z/12 independently of the slot (measured histogram in guard21/census21
runs); and all n steps must be multiply steps, (5/8)^n.  With W(k) = #{k-tuples of
Z/12 whose root sum vanishes}, a canonical word with multiplicities (k_1..k_m) has
probability (5/8)^n * R^{-(n-m)} * prod_i W(k_i) / 12^n, and each of the m! slot
labellings contributes the same, so the expected number of traps over S windows is

    E(S; R, d) = S * sum_{canonical words clean at depth d} m! * P(word).

Guard depth d: a multiply in slot h is converted to a cyclotomic squaring when h is
one of the previous d multiply slots, so a trapping word must have every repeat of a
slot at distance >= d+1.
"""
import json
import math
import os
import sys
from fractions import Fraction as F

SQ_NUM, SQ_DEN = 3, 8
EV = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")


def powers():
    """x^k reduced in Z[x]/Phi_12, k = 0..11"""
    pw = [[1, 0, 0, 0]]
    for _ in range(1, 12):
        a = pw[-1]
        pw.append([-a[3], a[0], a[1] + a[3], a[2]])
    return pw


PW = powers()


def reduce_counts(counts):
    v = [0, 0, 0, 0]
    for k, c in enumerate(counts):
        if c:
            for j in range(4):
                v[j] += c * PW[k][j]
    return v


def W_table(mmax=16):
    """W[m] = number of m-tuples of Z/12 whose root sum vanishes (0 for m<2)."""
    out = {0: 1}
    for m in range(1, mmax + 1):
        tot = 0
        counts = [0] * 12

        def rec(i, left):
            nonlocal tot
            if i == 11:
                counts[11] = left
                if reduce_counts(counts) == [0, 0, 0, 0]:
                    num = math.factorial(m)
                    for c in counts:
                        num //= math.factorial(c)
                    tot += num
                counts[11] = 0
                return
            for c in range(left + 1):
                counts[i] = c
                rec(i + 1, left - c)
            counts[i] = 0

        rec(0, m)
        out[m] = tot
    return out


W_TAB = W_table(8)


def canonical_words(n, d):
    """canonical (restricted-growth) slot words, no slot repeated within distance d"""
    out = []

    def rec(w, mx):
        if len(w) == n:
            out.append(tuple(w))
            return
        for s in range(mx + 2):
            if s in w[-d:]:
                continue
            w.append(s)
            rec(w, max(mx, s))
            w.pop()

    rec([], -1)
    return out


def word_prob(w, R, wtab=W_TAB):
    n = len(w)
    counts = {}
    for s in w:
        counts[s] = counts.get(s, 0) + 1
    num = 1
    for c in counts.values():
        num *= wtab.get(c, 0)
    if num == 0:
        return F(0), 0
    m = len(counts)
    # exact probability of a restricted-growth pattern under independent uniform slots:
    # falling factorial (R)_m / R^n (NOT R^-(n-m) * m!, which silently re-adds the m!
    # labellings that the canonical word already removes)
    ff = 1
    for i in range(m):
        ff *= (R - i)
    p = (F(SQ_DEN - SQ_NUM, SQ_DEN) ** n) * F(num, 12 ** n) * F(ff, R ** n)
    return p, m


def census(R, d, nmax):
    tot = F(0)
    per_n = {}
    surv = []
    for n in range(2, nmax + 1):
        s = F(0)
        for w in canonical_words(n, d):
            p, m = word_prob(w, R)
            if p:
                s += p
                surv.append({"n": n, "word": "".join("abcdefghijklmnop"[i] for i in w),
                             "slots": m, "log2_prob": round(math.log2(float(p)), 3)})
        if s:
            per_n[n] = s
        tot += s
    return tot, per_n, surv


def log2f(x):
    return round(math.log2(float(x)), 3) if x else None


def run(R, d, S_log2, nmax, label):
    tot, per_n, surv = census(R, d, nmax)
    traps_log2 = S_log2 + log2f(tot) if tot else None
    return {"label": label, "R": R, "depth": d, "nmax": nmax,
            "sum_word_prob_log2": log2f(tot),
            "expected_traps_log2": traps_log2,
            "per_length_log2": dict((str(k), log2f(v)) for k, v in per_n.items()),
            "survivors": sorted(surv, key=lambda r: -r["log2_prob"])[:12],
            "n_surviving_words_counted": len(surv)}


def four_table():
    """The 5 canonical length-4 slot words: does a depth-1 guard admit them?"""
    rows = []
    for w in canonical_words(4, 1):
        p1, m = word_prob(w, 4096)
        clean2 = all(w[i] not in w[max(0, i - 2):i] for i in range(4))
        rows.append({"word": "".join("abcd"[i] for i in w),
                     "slots": m,
                     "canonical_labellings": math.factorial(m),
                     "vanishing_shift_assignments": (lambda c: __import__("math").prod(
                         W_TAB[k] for k in c.values()))({s: w.count(s) for s in set(w)}),
                     "log2_prob_at_R4096": log2f(p1 * math.factorial(m)) if p1 else None,
                     "survives_depth1": True,
                     "survives_depth2": clean2})
    return rows


if __name__ == "__main__":
    S_LOG2 = float(sys.argv[1]) if len(sys.argv) > 1 else 125.8
    out = {"mode": "patterns", "S_log2": S_LOG2,
           "W_table": {str(k): v for k, v in W_TAB.items() if k <= 10},
           "phi12": "x^4 - x^2 + 1",
           "length4_table": four_table(),
           "configs": [
               run(4096, 1, S_LOG2, 8, "v20_shipped_R4096_d1"),
               run(4096, 2, S_LOG2, 10, "R4096_d2"),
               run(4096, 4, S_LOG2, 14, "R4096_d4"),
               run(65536, 5, S_LOG2, 16, "v21_shipped_R65536_d5"),
               run(65536, 6, S_LOG2, 16, "R65536_d6"),
           ]}
    path = os.path.join(EV, "patterns22.json")
    json.dump(out, open(path, "w"), indent=1)
    print(json.dumps({"length4": out["length4_table"],
                      "configs": [{k: c[k] for k in ("label", "depth", "sum_word_prob_log2",
                                                     "expected_traps_log2")}
                                  for c in out["configs"]]}, indent=1))
