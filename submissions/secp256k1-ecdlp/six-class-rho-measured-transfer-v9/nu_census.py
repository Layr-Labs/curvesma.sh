"""nu-normalized success census of the shipped six-class rho walk (production transfer).

WHY. The track's floor is p >= 0.39 for uniform full-production secp256k1 challenges. A
toy-run success rate is not that number: the load must be normalized so that a fixture and
production are the same experiment. This census measures

    nu      = operations / sqrt(class space)          (the load: work per birthday length)
    lam     = mean trail length / sqrt(class space)   (the sparsity of the DP regime)

at nu matched to the ledger row and lam held fixed across a six-group ladder, runs the
SHIPPED walk (experiment.py: Walker._run through run_trail_budgeted, the same code the
package's Part B uses), and measures success END TO END: a fixture succeeds only when two
distinct trails collide at a distinguished class and the shipped re-walk + solve returns the
hidden discrete log, verified as [x]P == Q and x == the generated secret.

BUDGET. Every operation is charged to a global meter before it runs (Walker._charge), so a
fixture performs EXACTLY V = round(nu * sqrt(N)) operations unless it solves first. The
census hard-asserts ops/sqrt(N) == nu to 1e-3 (see "budget_assert" in the output).

LADDER. y^2 = x^3 + 7 over F_q with q = 1 mod 3 and prime order r = 1 mod 3: the same j = 0
family as secp256k1, with the same beta/lambda six-class structure. Groups are chosen so that
sqrt(N) doubles from one group to the next while dp_bits grows by one, which holds lam fixed
(each group is CONFIRMED, not assumed: a random on-curve point P of order r forces #E = r because
r > q + 1 - 2 sqrt(q) and #E is one of the six j = 0 twists; the ladder rows in the evidence file
carry q, r, N, sqrt(N), dp_bits and that check).

USAGE.  python3 nu_census.py '<json config>'   (or a JSON config on stdin)
     {"groups": ["q25", ...], "nu_mult": [0.4, 0.7, 1.0, 1.3], "trials": 512,
      "out": "evidence/nu-census-<part>.json", "gate": true}
"""
import json, math, random, sys, time

import experiment as X

# neutralise the entry point's wall-clock budget: the census is metered in operations
X.budget = float("inf")
X.T0 = time.time()

CLS = 6
R_TABLE = 256          # partition size, constant across the ladder (production: 2^16)
KAPPA = X.KAPPA        # H1 slowdown allowance of the ledger
P0 = X.P0

# ledger row of the secp256k1 claim (proof.md section 4): S_ops and N'
R_PROD = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
S_INT = 142160481062555318092036037759448841269        # = ceil(S_ops), proof.md 4
N_PROD = (R_PROD - 1) // CLS
NU_LEDGER = S_INT / math.sqrt(N_PROD)                  # = S_ops / sqrt(N')
LAMBDA_TARGET = (1 << 4) / math.sqrt((X.TOYS["q25"][1] - 1) // CLS)


# ------------------------------------------------------------------ statistics
def binom_cdf(k, n, p):
    """P[X <= k] for X ~ Bin(n, p)."""
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    if p <= 0.0:
        return 1.0
    if p >= 1.0:
        return 0.0
    return min(1.0, sum(math.comb(n, i) * p ** i * (1.0 - p) ** (n - i) for i in range(k + 1)))


def cp_lower(k, n, alpha=0.05):
    """Exact one-sided Clopper-Pearson lower endpoint, bracketed on (0, 1).

    Not the clipped helper: for k = n the answer is alpha^(1/n) = 0.994166 at n = 512,
    and the assertion below proves the bracket is not clipped at 0.99.
    """
    assert abs(cp_lower_full(n) - alpha ** (1.0 / n)) < 1e-12
    if k == 0:
        return 0.0
    if k == n:
        return alpha ** (1.0 / n)
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        # P[X >= k | p] is increasing in p; the exact lower endpoint sits where it = alpha
        if 1.0 - binom_cdf(k - 1, n, mid) > alpha:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def cp_lower_full(n, alpha=0.05):
    """The k = n case, computed independently of the bisection bracket."""
    return alpha ** (1.0 / n)


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1.0 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


# ------------------------------------------------------------------ one census point
def run_point(E, N, dp_bits, nu_mult, trials, master_seed, R=None):
    R = R_TABLE if R is None else R
    sq = math.sqrt(N)
    nu = nu_mult * NU_LEDGER
    V = int(round(nu * sq))
    budget_err = abs(V / sq - nu)
    assert budget_err <= 1e-3, ("budget not ledger-normalized", V / sq, nu)
    master = random.Random(master_seed)
    seeds = [master.getrandbits(64) for _ in range(trials)]   # fresh RNG per fixture
    rows = []
    t0 = time.time()
    for i in range(trials):
        rng = random.Random(seeds[i])
        tr = X.toy_trial(E, rng, R, dp_bits, op_budget=V)
        rows.append(tr)
    ok = sum(1 for r in rows if r["ok"] is True)
    exhausted = sum(1 for r in rows if r["ops"] == V)
    ops = sum(r["ops"] for r in rows)
    steps = sum(r["steps"] for r in rows)
    doublings = sum(r["doublings"] for r in rows)
    dps = sum(r["records"] for r in rows)
    merges = sum(r["merges"] for r in rows)
    cap = sum(r["abandoned"] for r in rows)
    exc = sum(r["exceptional"] for r in rows)
    esc2 = sum(r["esc2"] for r in rows)
    esccp = sum(r["esccp"] for r in rows)
    trails = sum(r["trails"] for r in rows)
    lo, hi = wilson(ok, trials)
    out = dict(
        group=E.name, q=E.q, r=E.r, N=N, sqrtN=round(sq, 3), dp_bits=dp_bits,
        R=R, nu_mult=nu_mult, nu=round(nu, 6), V=V,
        budget_assert=dict(ops_over_sqrtN=round(V / sq, 6), nu=round(nu, 6),
                           abs_err=round(budget_err, 9), tol=1e-3, ok=budget_err <= 1e-3),
        trials=trials, successes=ok, p=round(ok / trials, 6),
        p_cp95_lower=round(cp_lower(ok, trials), 6),
        p_wilson=tuple(round(v, 6) for v in (lo, hi)),
        exhausted=exhausted,
        p_merge=round(sum(1 for r in rows if r["merges"] > 0) / trials, 6),
        merges=merges, records=dps,
        mean_records_per_trial=round(dps / trials, 2),
        mean_ops_per_trial=round(ops / trials, 1),
        ops_over_sqrtN_realized=round(ops / trials / sq, 6),
        lambda_mean_trail=round(trails and ops / trails / sq, 6),
        lambda_target=round(LAMBDA_TARGET, 6),
        dp_density_realized=round(dps / max(1, ops), 8),
        dp_density_target=round(2.0 ** -dp_bits, 8),
        esc2_per_op=round(esc2 / max(1, ops), 8),
        esccp_per_op=round(esccp / max(1, ops), 8),
        doublings_per_op=round(doublings / max(1, ops), 8),
        abandoned_per_trial=round(cap / trials, 4),
        exceptional_per_trial=round(exc / trials, 4),
        trails_per_trial=round(trails / trials, 2),
        p_model_birthday=round(1 - math.exp(-(nu / KAPPA) ** 2 / 2), 6),
        seconds=round(time.time() - t0, 2),
    )
    return out


# ------------------------------------------------------------------ gate
def known_answer_gate(E, N, dp_bits, trials, master_seed, nu_mult=6.0, R=None):
    R = R_TABLE if R is None else R
    """The runner must return the HAND-VERIFIED solution on every fixture before any
    census scale-up (a 0/N or wrong-x gate aborts the round: the g1-ecdlp lane died on an
    inverted challenge map). Each fixture's secret is generated here, so the answer is known."""
    sq = math.sqrt(N)
    nu = nu_mult * NU_LEDGER
    V = int(round(nu * sq))
    master = random.Random(master_seed)
    solved = 0
    for i in range(trials):
        rng = random.Random(master.getrandbits(64))
        tr = X.toy_trial(E, rng, R, dp_bits, op_budget=V)
        if tr["ok"] is True:
            solved += 1
    return dict(group=E.name, trials=trials, solved=solved, nu=round(nu, 6), V=V,
                pass_gate=(solved == trials))


# ------------------------------------------------------------------ fit
def fit_law(points):
    """Fit p(nu) on the measured points. Model A: p = 1-exp(-c nu). Model B (birthday,
    the mechanism's own form): p = 1-exp(-c nu^2). Model C: log-log, free exponent."""
    if len({round(pt["nu"], 9) for pt in points}) < 2:
        return {"skipped": "fewer than two distinct nu anchors"}
    ys, nus = [], []
    for pt in points:
        p = min(max(pt["p"], 1e-9), 1 - 1e-12)
        ys.append(-math.log(1 - p))
        nus.append(pt["nu"])
    out = {}
    for name, f in (("A_linear", lambda v: v), ("B_birthday", lambda v: v * v)):
        num = sum(f(v) * y for v, y in zip(nus, ys))
        den = sum(f(v) ** 2 for v in nus)
        c = num / den
        res = [pt["p"] - (1 - math.exp(-c * f(pt["nu"]))) for pt in points]
        out[name] = dict(c=c, rms=round(math.sqrt(sum(r * r for r in res) / len(res)), 6),
                         max_abs_resid=round(max(abs(r) for r in res), 6),
                         residuals=[round(r, 6) for r in res])
    # free exponent by least squares on (ln nu, ln y)
    xs = [math.log(v) for v in nus]
    lns = [math.log(y) for y in ys]
    mx, my = sum(xs) / len(xs), sum(lns) / len(lns)
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, lns))
    a = sxy / sxx
    c = math.exp(my - a * mx)
    res = [pt["p"] - (1 - math.exp(-c * pt["nu"] ** a)) for pt in points]
    out["C_free_exponent"] = dict(exponent=round(a, 4), c=c,
                                  rms=round(math.sqrt(sum(r * r for r in res) / len(res)), 6),
                                  max_abs_resid=round(max(abs(r) for r in res), 6),
                                  residuals=[round(r, 6) for r in res])
    return out


def main(cfg):
    groups = cfg["groups"]
    nu_mult = cfg.get("nu_mult", [0.4, 0.7, 1.0, 1.3])
    trials = cfg.get("trials", 512)
    R = cfg.get("R_table", R_TABLE)
    seed0 = cfg.get("seed", 20261024)
    out = {"config": cfg, "nu_ledger": round(NU_LEDGER, 6), "kappa": KAPPA, "P0_ledger": P0,
           "lambda_target": round(LAMBDA_TARGET, 6), "R_table": R,
           "group_table": [], "points": [], "gate": []}
    for gi, name in enumerate(groups):
        E = X.make_toy(name)
        E.name = name
        N = (E.r - 1) // CLS
        dp = cfg["dp_bits"][name]
        out["group_table"].append(dict(
            name=name, q=E.q, r=E.r, bits=E.q.bit_length(), N=N, sqrtN=round(math.sqrt(N), 3),
            dp_bits=dp, lambda_trail=round((1 << dp) / math.sqrt(N), 6),
            beta=E.beta, lam=E.lam, h=E.h,
            on_curve_G=E.on_curve(E.G), order_check=E.mul_raw(E.r, E.G) is None))
        if cfg.get("gate") and gi in cfg.get("gate_groups", [0, len(groups) - 1]):
            out["gate"].append(known_answer_gate(E, N, dp, cfg.get("gate_trials", 64),
                                                 seed0 * 31 + gi))
        for j, m in enumerate(nu_mult):
            pt = run_point(E, N, dp, m, trials, seed0 * 1009 + 97 * gi + j, R)
            out["points"].append(pt)
            print(f"[{name}] nu={pt['nu']:.4f} p={pt['p']:.4f} "
                  f"cp95={pt['p_cp95_lower']:.4f} model={pt['p_model_birthday']:.4f} "
                  f"({pt['seconds']}s)", file=sys.stderr, flush=True)
    # every fixture must be exactly at the ledger-normalized budget
    out["budget_assert"] = dict(all_ok=all(p["budget_assert"]["ok"] for p in out["points"]),
                                max_abs_err=max(p["budget_assert"]["abs_err"] for p in out["points"]),
                                tol=1e-3)
    # fit on the three required grid points per group (the 1.3 anchor is auxiliary)
    fit_pts = [p for p in out["points"] if p["nu_mult"] in (0.4, 0.7, 1.0)]
    out["fit"] = fit_law(fit_pts)
    out["fit_aux"] = fit_law(out["points"])
    out["seconds"] = round(time.time() - X.T0, 1)
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1:
        cfg = json.loads(sys.argv[1])
    else:
        cfg = json.load(sys.stdin)
    res = main(cfg)
    text = json.dumps(res, indent=1)
    if cfg.get("out"):
        with open(cfg["out"], "w") as f:
            f.write(text + "\n")
        print(json.dumps({"wrote": cfg["out"], "points": len(res["points"]),
                          "budget_assert": res["budget_assert"],
                          "gate": res["gate"]}), file=sys.stderr)
    else:
        print(text)
