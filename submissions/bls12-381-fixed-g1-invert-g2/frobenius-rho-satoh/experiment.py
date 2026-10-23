"""Toy check of the production construction (v3): W parallel walks, mixed
branch elements, orbit canonicalization, a single shared CHAINING hash
table checked on EVERY step (not a probe-capped table -- see proof.md for
why), a SQUARING escape (V <- V^2, (a,b) <- (2a,2b) mod n) on a fruitless
collision -- independently verified safer than an additive escape for this
order-6 automorphism (toy trap-rate comparison in proof.md) -- and linear
solve. Uses the SAME branch count L as production (not a smaller toy-only
L), and runs for EXACTLY the birthday-scale step budget T (scaled to the
toy group's size), reporting the measured success rate AT that budget, not
after a generous "run until recovery" timeout. Models the GT discrete
log's exponent group as Z_n under addition with a multiply-by-unit
automorphism (the unit standing in for GT's Frobenius eigenvalue). Does
not validate production field/pairing arithmetic, instruction costs, or
full-size timing -- see proof.md."""
import json
import sys
import math
import random


def birthday_budget(n, w, p_target):
    m = n / w
    ln_factor = math.log(1 / (1 - p_target))
    T = int(math.isqrt(int(2 * m * ln_factor)))
    while T * T < 2 * m * ln_factor:
        T += 1
    return T


def run(n, w, unit, L, W, trials, p_target, seed):
    half = w // 2
    assert pow(unit, half, n) == 1
    powers = [pow(unit, i, n) for i in range(half)]
    group = sorted(set(powers) | set((n - g) % n for g in powers))
    assert len(group) == w, (n, w, unit, group)

    T_budget = birthday_budget(n, w, p_target)

    rng = random.Random(seed)
    branch_c = [rng.randrange(1, n) for _ in range(L)]
    branch_d = [rng.randrange(1, n) for _ in range(L)]

    def canon(value, a, b):
        best = None
        for g in group:
            cand = ((g * value) % n, (g * a) % n, (g * b) % n)
            if best is None or cand[0] < best[0]:
                best = cand
        return best

    def step(value, a, b, x):
        h = value % L
        return canon((value + branch_c[h] + branch_d[h] * x) % n,
                      (a + branch_c[h]) % n, (b + branch_d[h]) % n)

    def escape_square(value, a, b, x):
        return canon((2 * value) % n, (2 * a) % n, (2 * b) % n)

    correct = wrong = within_budget_fail = 0
    total_steps = 0
    escapes_triggered = 0
    for _ in range(trials):
        x = rng.randrange(1, n)
        walks = []
        table = {}
        for _ in range(W):
            a0 = rng.randrange(1, n)
            b0 = rng.randrange(1, n)
            v0 = (a0 + b0 * x) % n
            s = canon(v0, a0, b0)
            walks.append(list(s))
            table[s[0]] = (s[1], s[2])

        recovered = None
        steps_this_trial = 0
        # Run for EXACTLY T_budget total steps (summed across all W walks),
        # matching the production charged budget -- not "until recovery."
        while steps_this_trial < T_budget and recovered is None:
            for wi in range(W):
                if steps_this_trial >= T_budget:
                    break
                value, a, b = walks[wi]
                value, a, b = step(value, a, b, x)
                steps_this_trial += 1
                total_steps += 1
                if value in table:
                    a0, b0 = table[value]
                    denom = (b0 - b) % n
                    if denom == 0:
                        escapes_triggered += 1
                        value, a, b = escape_square(value, a, b, x)
                        table[value] = (a, b)
                        walks[wi] = [value, a, b]
                        continue
                    cand = ((a - a0) * pow(denom, n - 2, n)) % n
                    recovered = (cand == x % n) or ((n - cand) % n == x % n)
                    walks[wi] = [value, a, b]
                    break
                else:
                    table[value] = (a, b)
                    walks[wi] = [value, a, b]
        if recovered is True:
            correct += 1
        elif recovered is False:
            wrong += 1
        else:
            within_budget_fail += 1  # budget exhausted with no collision yet
    return correct, wrong, within_budget_fail, total_steps, escapes_triggered, T_budget


if __name__ == "__main__":
    params = json.load(sys.stdin)
    correct, wrong, budget_fail, steps, esc, T_budget = run(
        params["order"], params["w"], params["unit"], params["L"],
        params["W"], params["trials"], params.get("p_target", 0.40),
        params.get("seed", 1))
    print(json.dumps({
        "group": "integers modulo n under addition, with multiply-by-unit automorphism",
        "order": params["order"], "w": params["w"], "L": params["L"],
        "W": params["W"], "trials": params["trials"],
        "T_budget_per_trial": T_budget,
        "correct": correct, "wrong_answers": wrong,
        "budget_exhausted_no_collision": budget_fail,
        "measured_success_rate_at_budget": correct / params["trials"],
        "total_steps": steps, "escapes_triggered": esc,
        "all_correct_no_wrong_answers": wrong == 0,
        "scope": "toy group only; tests the mixed-branch/canonicalization/"
                 "every-step-table-check/squaring-escape/linear-solve logic "
                 "exactly as used in production, including the SAME branch "
                 "count L and running for exactly the birthday-scale step "
                 "budget (not an unrestricted timeout); not production "
                 "field/pairing arithmetic or instruction-level cost bounds",
    }))
