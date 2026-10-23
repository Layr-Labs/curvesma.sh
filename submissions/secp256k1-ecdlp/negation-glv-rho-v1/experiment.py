"""Toy check of the production construction (v3): W parallel walks, mixed
branch points R_h=[c_h]P+[d_h]Q, orbit canonicalization, a SINGLE SHARED
hash table checked on EVERY step (not just distinguished points -- v2's
distinguished-point design let a fruitless cycle that never touches a
distinguished point trap a walk undetected; checking every step catches
any cycle, fruitless or not, the instant it recurs), fruitless-cycle
escape, and linear solve. Models the discrete-log exponent group as Z_n
under addition with a multiply-by-unit automorphism. Reports
b_nonzero_step_count (confirms b is not stuck, the v1 bug) and
escapes_triggered (confirms fruitless collisions are being caught and
escaped, not silently trapping a walk). Does not validate production
curve/field arithmetic, instruction costs, or full-size timing -- see
proof.md."""
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
    escape_c, escape_d = rng.randrange(1, n), rng.randrange(1, n)

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

    def escape(value, a, b, x):
        return canon((value + escape_c + escape_d * x) % n,
                      (a + escape_c) % n, (b + escape_d) % n)

    correct = wrong = timeout = 0
    total_steps = 0
    b_nonzero_step_count = 0
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
                if b != 0:
                    b_nonzero_step_count += 1
                if value in table:
                    a0, b0 = table[value]
                    denom = (b0 - b) % n
                    if denom == 0:
                        escapes_triggered += 1
                        value, a, b = escape(value, a, b, x)
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
    return correct, wrong, timeout, total_steps, b_nonzero_step_count, escapes_triggered


if __name__ == "__main__":
    params = json.load(sys.stdin)
    correct, wrong, timeout, steps, bnz, esc = run(
        params["order"], params["w"], params["unit"], params["W"],
        params["trials"], params.get("seed", 1))
    print(json.dumps({
        "group": "integers modulo n under addition, with multiply-by-unit automorphism",
        "order": params["order"], "w": params["w"], "W": params["W"],
        "trials": params["trials"],
        "correct": correct, "wrong_answers": wrong, "timeouts": timeout,
        "total_steps": steps, "b_nonzero_step_count": bnz,
        "escapes_triggered": esc,
        "all_correct_no_wrong_answers": wrong == 0,
        "scope": "toy group only; tests the mixed-branch/canonicalization/"
                 "every-step-table-check/escape/linear-solve logic exactly "
                 "as used in production, not production curve/field "
                 "arithmetic or instruction-level cost bounds",
    }))
