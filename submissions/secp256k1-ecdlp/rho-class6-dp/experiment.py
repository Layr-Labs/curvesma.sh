"""Toy end-to-end run of the class-6 distinguished-point rho of proof.md.

Group: y^2 = x^3 + 2 over F_524347, prime order n = 523969 (n = 1 mod 6, so the order-3
endomorphism (x,y) -> (beta*x, y) acts as multiplication by an order-3 scalar lambda).
Walk points are canonical representatives of the classes {+-lambda^i R}: minimal x of
{x, beta x, beta^2 x}, then even y.  The walk uses the fixed-rule cycle escape of proof.md section 2:
a walker that sees its current class among its last 8 classes (since its own start or last escape)
has found an f-cycle of length L <= 8; it moves to the cycle member with the smallest canonical x
and leaves it by one fixed escape step e(m) = canon([2]m) (doubling), a function of that member only.
The run does two things.  (1) `recover`: end-to-end scalar recovery by distinguished-point collisions.
(2) `trap_test`: counts walks that a cycle escape fails to free (never reach a distinguished point within the cap),
for the doubling escape and for the table-step escape used by an earlier version (control).
(3) `merge_test`: starts a second walker from a random point of a first walker's trail (fresh history,
possibly inside a cycle) and checks both reach the same distinguished point; the same test with the
old rule "step with the next table entry when the class repeats" (control) shows the test can fail.
Toy scope only: it does not validate any production cost or the production walk quality.
"""
import json
import math
import random
import sys

P, B, N = 524347, 2, 523969


def inv(a, m):
    return pow(a, -1, m)


def add(A, C):
    if A is None:
        return C
    if C is None:
        return A
    (x1, y1), (x2, y2) = A, C
    if x1 == x2:
        if (y1 + y2) % P == 0:
            return None
        s = 3 * x1 * x1 * inv(2 * y1, P) % P
    else:
        s = (y2 - y1) * inv(x2 - x1, P) % P
    x3 = (s * s - x1 - x2) % P
    return (x3, (s * (x1 - x3) - y1) % P)


def mul(k, A):
    R = None
    while k:
        if k & 1:
            R = add(R, A)
        A = add(A, A)
        k >>= 1
    return R


def setup():
    G = None
    for x in range(1, P):
        r = (x ** 3 + B) % P
        if pow(r, (P - 1) // 2, P) == 1:
            for y in range(1, P):
                if y * y % P == r:
                    G = (x, y)
                    break
            if G:
                break
    assert mul(N, G) is None
    beta = next(b for b in range(2, P) if pow(b, 3, P) == 1)
    lam = next(l for l in range(2, N) if pow(l, 3, N) == 1)
    if mul(lam, G) != (beta * G[0] % P, G[1]):
        beta = beta * beta % P
    assert mul(lam, G) == (beta * G[0] % P, G[1])
    return G, beta, lam


G, BETA, LAM = setup()
NCLS = (N - 1) // 6


def canon(R, beta=None):
    """Return (rep, s) with rep = s*R, s in {+-1, +-lam, +-lam^2} (scalar mod N)."""
    x, y = R
    cands = [(x, 1), (BETA * x % P, LAM), (BETA * BETA * x % P, LAM * LAM % N)]
    cx, s = min(cands)
    if y % 2:
        y = P - y
        s = (-s) % N
    return (cx, y), s


def make_table(rng, r, Q):
    T = []
    for _ in range(r):
        s, t = rng.randrange(N), rng.randrange(N)
        T.append((add(mul(s, G), mul(t, Q)), s, t))
    return T


class Walk:
    """Class walk f(c) = canon(c + T_h(c)), h(c) = bits [TB, TB+log2 r) of the canonical x of c."""

    def __init__(self, T, tb, ring=8):
        self.T, self.tb, self.r, self.ring = T, tb, len(T), ring
        self.mask = (1 << tb) - 1
        self.reentry = 0

    def move(self, R, a, b, shift=0):
        h = ((R[0] >> self.tb) + shift) % self.r
        Tj, sj, tj = self.T[h]
        R2 = add(R, Tj)
        if R2 is None:
            return None
        R2, s = canon(R2)
        return R2, (a + sj) * s % N, (b + tj) * s % N

    def double(self, R, a, b):
        """Doubling escape: canon([2]R) with coefficients (2a, 2b) rescaled by the unit applied."""
        R2 = add(R, R)
        if R2 is None:
            return None
        R2, s = canon(R2)
        return R2, 2 * a * s % N, 2 * b * s % N

    def run(self, R, a, b, cap, mode="fixed", trace=None):
        """Walk until a distinguished point or `cap` steps.  Returns (status, R, a, b, steps, events).
        mode 'fixed': proof.md rule (escape e(m) = canon([2]m)).  mode 'tadd': control, same cycle rule but
        the earlier escape by a table step e(m) = canon(m + T_(h(m)+1)).  mode 'old': control, history-dependent escape."""
        hist = [R[0]]
        steps = events = 0
        escaped_from = set()  # members m this walker has already escaped from (re-entry statistic)
        while steps < cap:
            if mode == "old":
                nxt = self.move(R, a, b, 1 if R[0] in hist[-(self.ring + 1):-1] else 0)
                if R[0] in hist[-(self.ring + 1):-1]:
                    events += 1
            else:
                nxt = self.move(R, a, b)
            steps += 1
            if nxt is None:
                return "inf", R, a, b, steps, events
            R, a, b = nxt
            if trace is not None:
                trace.append((R, a, b))
            if R[0] & self.mask == 0:
                return "dp", R, a, b, steps, events
            if mode in ("fixed", "tadd"):
                prev = hist[-self.ring:]
                if R[0] in prev:
                    # first repeat: the last L entries of prev are the cycle, prev[-L] is this point
                    L = len(prev) - prev.index(R[0])
                    cyc = prev[-L:]
                    j = cyc.index(min(cyc))
                    events += 1
                    self.reentry += cyc[j] in escaped_from  # escape returned to the same cycle
                    escaped_from.add(cyc[j])
                    for _ in range(j):  # walk to the cycle member with the smallest canonical x
                        nxt = self.move(R, a, b)
                        steps += 1
                        if nxt is None:
                            return "inf", R, a, b, steps, events
                        R, a, b = nxt
                        if trace is not None:
                            trace.append((R, a, b))
                    if mode == "fixed":
                        nxt = self.double(R, a, b)  # escape step e(m) = canon([2]m)
                    else:
                        nxt = self.move(R, a, b, 1)  # control: table step with the next entry
                    steps += 1
                    if nxt is None:
                        return "inf", R, a, b, steps, events
                    R, a, b = nxt
                    if trace is not None:
                        trace.append((R, a, b))
                    if R[0] & self.mask == 0:
                        return "dp", R, a, b, steps, events
                    hist = [R[0]]
                    continue
            hist.append(R[0])
            hist = hist[-(self.ring + 1):]
        return "cap", R, a, b, steps, events


def recover(x_secret, rng, tb=3, r=64, nu=1.03, mode="fixed"):
    Q = mul(x_secret, G)
    W = Walk(make_table(rng, r, Q), tb)
    cap = 40 << tb
    store = {}
    steps = 0
    while True:
        a, b = rng.randrange(N), rng.randrange(N)
        R = add(mul(a, G), mul(b, Q))
        if R is None:
            continue
        R, s = canon(R)
        a, b = a * s % N, b * s % N
        st, R, a, b, n, _ = W.run(R, a, b, cap, mode)
        steps += n
        if st == "dp":
            if R in store:
                a2, b2 = store[R]
                if (b - b2) % N:
                    xr = (a2 - a) * inv(b - b2, N) % N
                    return xr == x_secret and mul(xr, G) == Q, steps
            else:
                store[R] = (a, b)
        if steps > 40 * math.isqrt(NCLS) * 20:
            return False, steps


def merge_test(rng, trials, mode, tb=3, r=64):
    """Walker A from a random start; walker B from a random point of A's trail, with empty history."""
    x_secret = rng.randrange(1, N)
    Q = mul(x_secret, G)
    W = Walk(make_table(rng, r, Q), tb)
    cap = 40 << tb
    res = {"pairs": 0, "same_dp": 0, "different_dp": 0, "one_cut": 0, "with_cycle_event": 0}
    while res["pairs"] < trials:
        a, b = rng.randrange(N), rng.randrange(N)
        R = add(mul(a, G), mul(b, Q))
        if R is None:
            continue
        R, s = canon(R)
        a, b = a * s % N, b * s % N
        trace = [(R, a, b)]
        st, RA, _, _, _, ev = W.run(R, a, b, cap, mode, trace)
        if st != "dp" or len(trace) < 2:
            continue
        t = rng.randrange(len(trace) - 1)  # any trail point before the distinguished point
        RB, aB, bB = trace[t]
        stB, RB2, _, _, _, _ = W.run(RB, aB, bB, cap, mode)
        res["pairs"] += 1
        res["with_cycle_event"] += ev > 0
        if stB != "dp":
            res["one_cut"] += 1
        elif RB2 == RA:
            res["same_dp"] += 1
        else:
            res["different_dp"] += 1
    return res


def trap_test(rng, walks, mode, tb=4, r=16):
    """Fresh walks on a 16-entry table (many short cycles): how many saw a cycle event, and how many of those
    never reached a distinguished point within the 40/theta cap (trapped: the escape returned to its cycle)."""
    Q = mul(rng.randrange(1, N), G)
    W = Walk(make_table(rng, r, Q), tb)
    cap = 40 << tb
    res = {"walks": walks, "with_cycle_event": 0, "escape_reentered_cycle": 0, "trapped_at_cap": 0}
    for _ in range(walks):
        a, b = rng.randrange(N), rng.randrange(N)
        R = add(mul(a, G), mul(b, Q))
        if R is None:
            continue
        R, s = canon(R)
        st, _, _, _, _, ev = W.run(R, a * s % N, b * s % N, cap, mode)
        res["with_cycle_event"] += ev > 0
        res["trapped_at_cap"] += st == "cap"
    res["escape_reentered_cycle"] = W.reentry
    return res


def main():
    raw = sys.stdin.read().strip()
    cfg = json.loads(raw) if raw else {}
    trials = int(cfg.get("trials", 400))
    rng = random.Random(12345)
    ok, tot, within = 0, 0, 0
    k0 = math.sqrt(-2 * math.log(0.61) * NCLS)
    for i in range(trials):
        good, st = recover(rng.randrange(1, N), rng)
        ok += good
        tot += st
        within += good and st <= 1.03 * k0 + 16
    fixed = merge_test(rng, 4 * trials, "fixed", tb=4, r=16)  # r=16: many short cycles
    control = merge_test(rng, 4 * trials, "old", tb=4, r=16)
    trap_fixed = trap_test(rng, 50 * trials, "fixed")
    trap_tadd = trap_test(rng, 50 * trials, "tadd")
    print(json.dumps({"group": "y^2=x^3+2 over F_524347", "order": N, "classes": NCLS, "trials": trials,
                      "recovered_correct": ok, "mean_steps": round(tot / trials, 1),
                      "ideal_mean_steps_sqrt(pi*N'/2)": round(math.sqrt(math.pi * NCLS / 2), 1),
                      "within_1.03k0_plus_2_over_theta_steps": within,
                      "merge_test_table16": "second walker restarted at a random point of the first walker's trail",
                      "merge_test_fixed_rule": fixed, "merge_test_old_rule_control": control,
                      "trap_test_table16_doubling_escape": trap_fixed,
                      "trap_test_table16_table_step_escape_control": trap_tadd,
                      "scope": "toy group only; does not validate production cost or walk quality"}))


if __name__ == "__main__":
    main()
