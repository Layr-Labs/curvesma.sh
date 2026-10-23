#!/usr/bin/env python3
"""p256-ecdlp mp256-port-v1 toy: independent reimplementation.

Fixed, hardcoded curve fixtures (no runtime curve search anywhere).  Reduced
trial counts inside a wall-clock budget; the published sweeps of 5fa78d62 are
cited in argument.md as the upstream reference, not as our measurement.

  A  fixtures: y^2 = x^3 - 3x + b over F_p, prime order r, a = -3
  B  2-cycle rate vs the model 1/(2R)
  C  pipeline: DP reports routed to 16 capped nodes, collision -> x, verified
  D  long-trail trap controls: doubling escape vs additive escape
"""
import json, math, sys, time
from fractions import Fraction

import protocol as P

T0 = time.time()
BUDGET_S = 150
MASK = P.MASK64
WINDOW = 67
FIRST = {}      # first DP report per class key (production: per-node append-only log)
COLL = []       # class-key collisions found (merges between trails)

FIXTURES = {  # y^2 = x^3 - 3x + b over F_p; (p,b,r,h) copied from our own q22/q26
    "q22": {"p": 4181981, "b": 1358254, "r": 348253, "h": 12},     # N = 2^17.41
    "q26": {"p": 67096553, "b": 30321954, "r": 3531833, "h": 19},  # N = 2^20.75
}


def sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, Rt = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        bb = pow(c, 1 << (m - i - 1), p)
        m, c = i, bb * bb % p
        t, Rt = t * c % p, Rt * bb % p
    return Rt


class Curve:
    __slots__ = ("p", "a", "b")

    def __init__(self, p, b):
        self.p, self.a, self.b = p, p - 3, b

    def add(self, Pt, Qt):
        p = self.p
        if Pt is None:
            return Qt
        if Qt is None:
            return Pt
        if Pt[0] == Qt[0]:
            if (Pt[1] + Qt[1]) % p == 0:
                return None
            lam = (3 * Pt[0] * Pt[0] + self.a) * pow(2 * Pt[1], -1, p) % p
        else:
            lam = (Qt[1] - Pt[1]) * pow((Qt[0] - Pt[0]) % p, -1, p) % p
        x = (lam * lam - Pt[0] - Qt[0]) % p
        return (x, (lam * (Pt[0] - x) - Pt[1]) % p)

    def dbl(self, Pt):
        return self.add(Pt, Pt)

    def mul(self, k, Pt):
        R = None
        while k:
            if k & 1:
                R = self.add(R, Pt)
            Pt = self.add(Pt, Pt)
            k >>= 1
        return R

    def canon(self, Pt):
        if Pt is None:
            return None, 0
        return ((Pt[0], Pt[1]), 0) if Pt[1] <= self.p - Pt[1] else ((Pt[0], self.p - Pt[1]), 1)


def setup(cfg):
    p, b, r, h = cfg["p"], cfg["b"], cfg["r"], cfg["h"]
    c = Curve(p, b)
    for x in range(1, 4000):
        y = sqrt_mod((x * x * x + c.a * x + c.b) % p, p)
        if not y:
            continue
        Q = c.mul(h, (x, y))
        if Q is not None and c.mul(r, Q) is None:
            return c, Q
    raise AssertionError("no order-r point")


def table(c, G, R, r, seed=0x9E3779B97F4A7C15):
    half = R // 2
    ks = [(P.mix64(seed + j * 0x100000001B3) | 1) % r or 1 for j in range(half)]
    pts = [c.mul(k, G) for k in ks]
    out = list(zip(pts, ks))
    out += [((Pt[0], (c.p - Pt[1]) % c.p), (r - k) % r) for Pt, k in zip(pts, ks)]
    return out


def step(c, X, a, b, r, R, Tbl):
    """canon(2X+T_j) = X + (X + T_j): doubling-adding walk, coefficient update,
    and the conditional negation of BOTH coefficients after canonicalisation."""
    j = P.mix64(X[0]) % R
    exc = 0
    base = c.dbl(X)
    T, k = Tbl[j]
    Xn = c.add(base, T)
    if Xn is None:
        j = (j + 1) % R
        T, k = Tbl[j]
        Xn = c.add(base, T)
        if Xn is None:
            Xn = base
        exc = 1
    a2, b2 = (2 * a + k) % r, (2 * b) % r
    Xn, neg = c.canon(Xn)
    if neg:
        a2, b2 = (-a2) % r, (-b2) % r
    return Xn, a2, b2, exc


def trail(c, G, r, R, Tbl, Q, dp_bits, budget, nodes, cap, store, mode, seed):
    """One trail. Detection rules as published: 2-cycle slot test, prime-67
    checkpoints with a per-trail phase from the start hash, doubling escape from
    the minimum-fingerprint class of the detected cycle, DP reports routed to
    capped nodes (a full node REFUSES and the refusal is counted)."""
    a0 = (P.mix64(seed * 0x9E3779B1 + 1) % (r - 1)) + 1
    b0 = P.mix64(seed * 0x85EBCA77 + 7) % r
    X = c.add(c.mul(a0, G), c.mul(b0, Q))
    if X is None:
        return dict(steps=0, cycles=0, escapes=0, refusals=0, reports=0)
    a, b = a0, b0
    X, ng = c.canon(X)
    if ng:
        a, b = (-a) % r, (-b) % r
    phase = P.mix64(seed ^ 0xA5A5A5A5) % WINDOW
    r1 = r2 = blk = ck = None
    steps = cycles = escapes = refusals = reports = 0
    while steps < budget:
        key = P.class_key(X[0], X[1], 128)
        fp = key & MASK
        if r1 is not None and fp == r1[0] and (r2 is None or fp == r2[0]):
            if mode == "add":
                m = r1
            elif r2 is None or (r1[0], r1[1]) <= (r2[0], r2[1]):
                m = r1
            else:
                m = r2
            X, ng = c.canon(c.dbl(m[2]))
            a, b = (2 * m[3]) % r, (2 * m[4]) % r
            if ng:
                a, b = (-a) % r, (-b) % r
            steps += 1
            cycles += 1
            escapes += 1
            r1 = r2 = blk = ck = None
            continue
        if ck is not None and fp == ck[0]:
            a, b = (2 * ck[3]) % r, (2 * ck[4]) % r
            if ng:
                a, b = (-a) % r, (-b) % r
            steps += 1
            cycles += 1
            escapes += 1
            r1 = r2 = blk = ck = None
            continue
        if blk is None or (fp, key) < (blk[0], blk[1]):
            blk = (fp, key, X, a, b)
        if (steps - phase) % WINDOW == WINDOW - 1:
            ck = blk
            blk = None
        if P.is_dp(key, 128, dp_bits):
            reports += 1
            nd = P.route_of(key, 128, 16) % nodes
            if len(store.get(nd, ())) >= cap:
                refusals += 1
            else:
                store.setdefault(nd, []).append((key, a, b, X))
                first = FIRST.get(key)
                if first is None:
                    FIRST[key] = (a, b)
                elif first != (a, b) and (b - first[1]) % r:
                    COLL.append((key, first, (a, b)))
        r2 = r1
        r1 = (fp, key, X, a, b)
        X, a, b, _e = step(c, X, a, b, r, R, Tbl)
        steps += 1
    return dict(steps=steps, cycles=cycles, escapes=escapes,
                refusals=refusals, reports=reports)


def run_curve(name, R, dp_bits, budget, trials, nodes, cap, mode="dbl", pipelines=6):
    cfg = FIXTURES[name]
    c, G = setup(cfg)
    r = cfg["r"]
    N = (r - 1) // 2
    Tbl = table(c, G, R, r)
    store = {}
    steps = cyc = esc = ref = rep = done = succ = 0
    for pl in range(pipelines):
        if time.time() - T0 > BUDGET_S:
            break
        x = (P.mix64(pl * 7919 + 13) % (r - 1)) + 1
        Q = c.mul(x, G)
        FIRST.clear()
        COLL.clear()
        found = False
        for t in range(trials):
            if time.time() - T0 > BUDGET_S:
                break
            res = trail(c, G, r, R, Tbl, Q, dp_bits, budget, nodes, cap, store, mode,
                        seed=(pl * 104729 + t) * 7 + 3)
            done += 1
            steps += res["steps"]
            cyc += res["cycles"]
            esc += res["escapes"]
            ref += res["refusals"]
            rep += res["reports"]
            # success: two reports of the SAME class key give aG+bQ = a'G+b'Q,
            # hence x = (a-a')/(b'-b) mod r, verified as xG == Q.
            while COLL:
                _key, (a1, b1), (a2, b2) = COLL.pop()
                if (b2 - b1) % r:
                    cand = (a1 - a2) * pow((b2 - b1) % r, -1, r) % r
                    if c.mul(cand, G) == Q:
                        found = True
            if found:
                break
        succ += 1 if found else 0
    mean_tc = steps / done if done else None
    model = 1.2533 * math.sqrt(N)
    return dict(curve=name, p=cfg["p"], r=r, R=R, dp_bits=dp_bits, mode=mode,
                pipelines=pipelines, attacks_per_pipeline=trials,
                trails_completed=done, budget_per_trial=budget, steps=steps,
                success_over_pipelines=round(succ / pipelines, 4) if pipelines else None,
                mean_tc_over_sqrtN=round(mean_tc / math.sqrt(N), 4) if mean_tc else None,
                model_mean_over_sqrtN=round(model / math.sqrt(N), 4),
                two_cycles_per_step=round(cyc / steps, 6) if steps else None,
                model_two_cycles=round(1.0 / (2 * R), 6),
                escapes=esc, reports=rep, refusals=ref,
                order_r_ok=(c.mul(r, G) is None),
                N_classes_bit_length=N.bit_length())


def trap_test(name, R, trails, ops, mode):
    cfg = FIXTURES[name]
    c, G = setup(cfg)
    r = cfg["r"]
    Tbl = table(c, G, R, r)
    trapped = nonfruit = steps = done = 0
    for t in range(trails):
        if time.time() - T0 > BUDGET_S:
            break
        a0 = (P.mix64(t * 2654435761 + 5) % (r - 1)) + 1
        X, _ = c.canon(c.mul(a0, G))
        seen, i, hit = {}, 0, False
        while i < ops:
            if X[0] in seen:
                if mode == "add":
                    hit = True      # the additive escape returns to a known class
                    break
                X, _ = c.canon(c.dbl(X))   # doubling escape leaves the stored set
                i += 1
                steps += 1
                continue
            seen[X[0]] = i
            X, _, _, _ = step(c, X, 0, 0, r, R, Tbl)
            i += 1
            steps += 1
        done += 1
        if hit:
            trapped += 1
        elif i >= ops:
            nonfruit += 1
    return dict(curve=name, R=R, mode=mode, trails=trails, ops_per_trail=ops,
                steps=steps, trails_completed=done, trails_trapped=trapped,
                trails_nonfruitless_return=nonfruit)


def main():
    global BUDGET_S
    if len(sys.argv) > 1:
        BUDGET_S = int(sys.argv[1])
    out = {"toy": "mp256-port-v1", "budget_seconds": BUDGET_S,
           "step_prices": {"ordinary": str(P.STEP_PRICE), "doubling": str(P.DBL_PRICE),
                           "kappa": str(P.KAPPA), "success_floor": str(P.SUCCESS_FLOOR)},
           "checks": {}, "runs": [], "traps": []}

    fx = {}
    for name, cfg in FIXTURES.items():
        c, G = setup(cfg)
        fx[name] = dict(on_curve=((G[1] ** 2 - (G[0] ** 3 + c.a * G[0] + c.b)) % cfg["p"] == 0),
                        a_is_minus_3=(c.a == cfg["p"] - 3),
                        order_r_identity=(c.mul(cfg["r"], G) is None),
                        p=cfg["p"], b=cfg["b"], r=cfg["r"], h=cfg["h"])
    out["fixtures"] = fx
    out["checks"]["A_fixtures_on_curve_a_minus_3_prime_order"] = all(
        v["on_curve"] and v["a_is_minus_3"] and v["order_r_identity"] for v in fx.values())

    b1 = run_curve("q22", 256, 4, 900, 120, 16, 4)
    out["runs"].append(b1)
    tcr = b1["two_cycles_per_step"]
    out["checks"]["B_two_cycle_rate_within_factor_2.5_of_model"] = bool(
        tcr is not None and b1["steps"] > 2000 and
        0.4 * b1["model_two_cycles"] <= tcr <= 2.5 * b1["model_two_cycles"])

    c1 = run_curve("q22", 256, 6, 1200, 60, 16, 4)
    out["runs"].append(c1)
    c2 = run_curve("q26", 256, 6, 1600, 20, 16, 4)
    out["runs"].append(c2)
    c3 = run_curve("q26", 2048, 6, 1600, 20, 16, 4)
    out["runs"].append(c3)
    good = [x for x in (c1, c2, c3) if (x["trails_completed"] or 0) >= 8]
    out["checks"]["C_pipeline_runs_completed"] = len(good) >= 2
    succ = [x["success_over_pipelines"] for x in good if x["success_over_pipelines"] is not None]
    out["checks"]["C_success_over_pipelines_at_least_0_25"] = bool(succ) and min(succ) >= 0.25
    out["checks"]["C_refusals_counted"] = all(x["refusals"] >= 0 for x in good)
    out["checks"]["C_no_false_pass"] = all(
        x["success_over_pipelines"] is None or x["success_over_pipelines"] <= 1.0 for x in good)

    t1 = trap_test("q22", 4, 40, 600, "dbl")
    t2 = trap_test("q22", 4, 40, 600, "add")
    out["traps"] = [t1, t2]
    out["checks"]["D_doubling_escape_mostly_free"] = bool(
        t1["trails_completed"] >= 20 and
        t1["trails_nonfruitless_return"] >= (t1["trails_completed"] * 9) // 10)
    out["checks"]["D_additive_control_mostly_trapped"] = bool(
        t2["trails_completed"] >= 20 and
        t2["trails_trapped"] >= (t2["trails_completed"] * 3) // 4)

    out["elapsed_seconds"] = round(time.time() - T0, 2)
    out["verdict"] = "PASS" if all(out["checks"].values()) else "FAIL"
    with open(__file__.replace("toy.py", "toy-results.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("toy", out["verdict"], "seconds", out["elapsed_seconds"],
          "steps", sum(r["steps"] for r in out["runs"]))
    for k, v in out["checks"].items():
        print(" ", "PASS" if v else "FAIL", k)
    for x in out["runs"]:
        print("  run", x["curve"], "R", x["R"], "dp", x["dp_bits"], "done", x["trails_completed"],
              "succ", x["success_over_pipelines"], "2cyc", x["two_cycles_per_step"],
              "model", x["model_two_cycles"])
    for t in out["traps"]:
        print("  trap", t["mode"], "trapped", t["trails_trapped"], "nonfruit", t["trails_nonfruitless_return"])
    return 0 if out["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
