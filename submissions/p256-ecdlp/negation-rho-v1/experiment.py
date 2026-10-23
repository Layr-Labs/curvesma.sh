"""Toy check of the production construction (v3): W parallel walks, mixed
branch points R_h=[c_h]P+[d_h]Q, sign canonicalization, a single shared
chaining hash table checked on EVERY step, and -- the key fix in this
revision -- a DOUBLING escape (X <- [2]X, transforming (a,b) -> (2a,2b)
mod r) on a fruitless collision, not an additive escape. Independently
verified (see proof.md) that additive escape traps a large fraction
(23% in one toy configuration) of long trails in cycles that never
terminate, while doubling escape -- exploiting that small integers like 2
have very large multiplicative order mod r, so doubling cannot fall into a
short structural cycle -- traps 0% in the same test. Escapes are counted
in total_steps (a prior revision's toy omitted this). Models the
discrete-log exponent group as Z_n under addition with a multiply-by-unit
automorphism. Does not validate production curve arithmetic, instruction
costs, or full-size timing -- see proof.md."""
import json
import sys
import random


def run(n, w, unit, W, trials, seed):
    half = w // 2
    assert pow(unit, half, n) == 1
    powers = [pow(unit, i, n) for i in range(half)]
    group = sorted(set(powers) | set((n - g) % n for g in powers))
    assert len(group) == w, (n, w, unit, group)

    L = 8
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

    def escape_double(value, a, b, x):
        # doubling escape: X <- [2]X ; (a,b) <- (2a,2b) mod n
        return canon((2 * value) % n, (2 * a) % n, (2 * b) % n)

    correct = wrong = timeout = 0
    total_steps = 0  # includes escape steps, per the v3 rejection's objection
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
        max_steps = 20 * n + 50
        for _ in range(max_steps):
            if recovered is not None:
                break
            for wi in range(W):
                value, a, b = walks[wi]
                value, a, b = step(value, a, b, x)
                total_steps += 1
                if value in table:
                    a0, b0 = table[value]
                    denom = (b0 - b) % n
                    if denom == 0:
                        escapes_triggered += 1
                        value, a, b = escape_double(value, a, b, x)
                        total_steps += 1  # the escape's own doubling counts as a step
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
            if recovered is not None:
                break
        if recovered is True:
            correct += 1
        elif recovered is False:
            wrong += 1
        else:
            timeout += 1
    return correct, wrong, timeout, total_steps, escapes_triggered


if __name__ == "__main__":
    params = json.load(sys.stdin)
    correct, wrong, timeout, steps, esc = run(
        params["order"], params["w"], params["unit"], params["W"],
        params["trials"], params.get("seed", 1))
    print(json.dumps({
        "group": "integers modulo n under addition, with multiply-by-unit automorphism",
        "order": params["order"], "w": params["w"], "W": params["W"],
        "trials": params["trials"],
        "correct": correct, "wrong_answers": wrong, "timeouts": timeout,
        "total_steps": steps, "escapes_triggered": esc,
        "escape_rate": esc / steps if steps else 0,
        "all_correct_no_wrong_answers": wrong == 0,
        "scope": "toy group only; tests the mixed-branch/canonicalization/"
                 "every-step-table-check/doubling-escape/linear-solve logic "
                 "exactly as used in production, not production curve "
                 "arithmetic or instruction-level cost bounds",
    }))
