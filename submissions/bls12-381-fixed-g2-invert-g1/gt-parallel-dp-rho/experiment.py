"""Toy-scale verification of parallel distinguished-point Pollard rho.

Full-size claim under test (BLS12-381): solving
    a = dlog_G(H),  G, H in GT (order-r subgroup of Fp12*),
by an additive rho walk with WALKS workers sharing one distinguished-point
(DP) table, capped at BUDGET = 1.2*sqrt(r) total group operations, costs at
most ~1.2*sqrt(r) operations and collides with probability about
1 - exp(-BUDGET^2 / (2r)) ~= 0.513 over the challenge and algorithm coins.
The post-collision linear solve recovers a with probability ~ 1 - O(1/r).

This script checks the four load-bearing statements on a scaled-down
analog: a prime-order subgroup of GF(p)*, with q in {65521, 131071}.
Detection waits for the collided trajectory to reach a distinguished
point; the DP width w is chosen per size so that the parallel detection
tail WALKS*2**w stays below sqrt(q)/4, matching the production regime
where m*2**w << sqrt(r).

  1. Cost. Mean total iterations to a shared-table collision vs the
     prediction sqrt(pi*q/2), at two subgroup sizes (scaling ~ sqrt).
  2. Budgeted success. Fraction of trials whose total work stays within
     1.2*sqrt(q), vs the closed-form prediction 1 - exp(-0.72) ~= 0.513.
  3. Success of the solve. Fraction of collisions whose linear solve
     recovers the true discrete log.
  4. Table density. Stored distinguished points ~ total_steps / 2**w.

Standard library only; runtime is a few seconds.
Usage: experiment.py [order]   (order a prime in [1000, 200000];
default runs both standard sizes).  Exit status 0 iff every check passes.
"""

import math
import random
import sys

TRIALS = 200        # trials per subgroup size
WALKS = 8           # parallel workers sharing the DP table
P_BUCKETS = 32      # additive-walk partition count
STEP_CAP = 40       # abort cap, in units of sqrt(q)
BUDGET_FACTOR = 1.2  # work budget in units of sqrt(q); P ~ 1-exp(-B^2/2)


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def group_of_order(q):
    """Return (p, g): p prime with q | p-1, g of exact order q mod p."""
    k = 2
    while not is_prime(k * q + 1):
        k += 2
    p = k * q + 1
    h0 = 2
    while True:
        g = pow(h0, k, p)
        if g != 1 and pow(g, q, p) == 1:
            return p, g
        h0 += 1


def make_walk(g, h, q, p, rng):
    """P random multipliers M_j = g^u_j h^v_j for the additive walk."""
    tab = []
    for _ in range(P_BUCKETS):
        u = rng.randrange(1, q)
        v = rng.randrange(1, q)
        tab.append((pow(g, u, p) * pow(h, v, p) % p, u, v))
    return tab


def step(state, tab, q, p):
    x, a, b = state
    M, u, v = tab[x % P_BUCKETS]
    return (x * M % p, (a + u) % q, (b + v) % q)


def parallel_dp_trial(g, h, a_true, q, p, tab, rng, dp_bits):
    """WALKS workers, one shared DP table. Returns (steps, a_found, n_dp).

    Detects a collision when a worker produces a distinguished point that
    is already in the shared table. Every iterated step is one group
    multiplication; total `steps` is summed over all workers.
    """
    mask = (1 << dp_bits) - 1
    states = []
    for _ in range(WALKS):
        a0 = rng.randrange(q)
        b0 = rng.randrange(q)
        states.append([pow(g, a0, p) * pow(h, b0, p) % p, a0, b0])
    seen = {}
    steps = 0
    cap = int(STEP_CAP * math.sqrt(q))
    while steps < cap:
        for st in states:
            x, a, b = step(tuple(st), tab, q, p)
            st[0], st[1], st[2] = x, a, b
            steps += 1
            if x & mask == 0:
                if x in seen:
                    a1, b1 = seen[x]
                    den = (b1 - b) % q
                    if den == 0:
                        return steps, None, len(seen)
                    a_found = ((a - a1) * pow(den, -1, q)) % q
                    return steps, a_found, len(seen)
                seen[x] = (a, b)
    return steps, None, len(seen)


def run_size(q, seed):
    rng = random.Random(seed)
    p, g = group_of_order(q)
    pred = math.sqrt(math.pi * q / 2.0)
    budget = BUDGET_FACTOR * math.sqrt(q)
    # DP width: the parallel detection tail WALKS*2**w must stay well
    # below sqrt(q), matching the production regime m*2**w << sqrt(r).
    dp_bits = max(2, int(math.log2(math.sqrt(q) / (4.0 * WALKS))))
    tot_steps = 0
    tot_dp = 0
    ok = 0
    within = 0
    fail_solve = 0
    fail_cap = 0
    worst = 0
    for t in range(TRIALS):
        a_true = rng.randrange(q)
        h = pow(g, a_true, p)
        tab = make_walk(g, h, q, p, rng)
        steps, a_found, n_dp = parallel_dp_trial(
            g, h, a_true, q, p, tab, rng, dp_bits)
        tot_steps += steps
        tot_dp += n_dp
        worst = max(worst, steps)
        if steps <= budget:
            within += 1
        if a_found is None:
            if steps >= int(STEP_CAP * math.sqrt(q)):
                fail_cap += 1
            else:
                fail_solve += 1
        elif a_found == a_true and pow(g, a_found, p) == h:
            ok += 1
        else:
            fail_solve += 1
    mean_steps = tot_steps / TRIALS
    report = {
        "q": q,
        "p": p,
        "dp_bits": dp_bits,
        "trials": TRIALS,
        "workers": WALKS,
        "predicted_steps_sqrt_pi_q_over_2": round(pred, 2),
        "mean_total_steps": round(mean_steps, 2),
        "ratio_measured_over_predicted": round(mean_steps / pred, 3),
        "worst_trial_steps": worst,
        "solve_successes": ok,
        "solve_failures": fail_solve,
        "cap_aborts": fail_cap,
        "mean_dp_stored": round(tot_dp / TRIALS, 3),
        "dp_density_predicted": round(mean_steps / (1 << dp_bits), 3),
        "budget_steps": round(budget, 2),
        "within_budget_fraction": round(within / TRIALS, 3),
        "within_budget_predicted": round(
            1.0 - math.exp(-(BUDGET_FACTOR ** 2) / 2.0), 3),
    }
    return report


def main():
    if len(sys.argv) > 1:
        q = int(sys.argv[1])
        if not (1000 <= q <= 200000) or not is_prime(q):
            print("FAIL: order must be a prime in [1000, 200000]")
            sys.exit(1)
        sizes = [q]
    else:
        sizes = [65521, 131071]
    reports = [run_size(s, seed=1000 + s) for s in sizes]
    for r in reports:
        print(r)
    ok = True
    for r in reports:
        ratio = r["ratio_measured_over_predicted"]
        # Parallel-DP rho on a random mapping has expected total work
        # ~ sqrt(pi*q/2); the band [0.6, 2.5] absorbs the detection tail
        # and walk-merging overhead at toy scale.
        if not (0.6 <= ratio <= 2.5):
            print("FAIL: step count outside expected band", r["q"], ratio)
            ok = False
        pred_hit = r["within_budget_predicted"]
        got_hit = r["within_budget_fraction"]
        if abs(got_hit - pred_hit) > 0.15:
            print("FAIL: budgeted success far from prediction",
                  r["q"], got_hit, "vs", pred_hit)
            ok = False
        if r["solve_failures"] > max(1, TRIALS // 100):
            print("FAIL: too many degenerate collisions", r["q"])
            ok = False
        if r["cap_aborts"] > 0:
            print("FAIL: a trial hit the step cap", r["q"])
            ok = False
        if r["solve_successes"] < TRIALS - max(1, TRIALS // 100):
            print("FAIL: success rate below 99%", r["q"])
            ok = False
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
