"""Toy check of the production construction: W parallel walks, mixed branch
points R_h = [c_h]P + [d_h]Q (both P- and Q-multiples, matching production
exactly -- this is what the v1 package got wrong: it tested this correct
mixed-branch mechanism in the toy but mis-described the production walk as
using only P-multiples, which would leave b permanently at 0 and make the
construction unable to ever solve for x), orbit canonicalization,
distinguished-point storage (not full per-step history, matching the
production memory model), fruitless-cycle escape, and linear solve.
Models <P> as Z_n under addition with a multiply-by-unit automorphism.
Reports b_nonzero_step_count explicitly so a reviewer can confirm b is not
stuck at 0. Does not validate production curve arithmetic, instruction
costs, or full-size timing -- see proof.md."""
import json
import sys
import random


def run(n, w, unit, W, dp_mod, trials, seed):
    half = w // 2
    assert pow(unit, half, n) == 1
    powers = [pow(unit, i, n) for i in range(half)]
    group = sorted(set(powers) | set((n - g) % n for g in powers))
    assert len(group) == w, (n, w, unit, group)

    L = 8
    rng = random.Random(seed)
    branch_c = [rng.randrange(1, n) for _ in range(L)]
    branch_d = [rng.randrange(1, n) for _ in range(L)]   # mixed: P- AND Q-multiples
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

    def is_distinguished(value):
        return value % dp_mod == 0

    correct = wrong = timeout = 0
    total_steps = 0
    b_nonzero_step_count = 0
    for _ in range(trials):
        x = rng.randrange(1, n)
        walks = []
        for _ in range(W):
            a0 = rng.randrange(1, n)
            b0 = rng.randrange(1, n)           # mixed start: b0 != 0
            v0 = (a0 + b0 * x) % n
            walks.append(list(canon(v0, a0, b0)))

        dp_table = {}
        recovered = None
        max_steps = 20 * n + 50
        for _ in range(max_steps):
            if recovered is not None:
                break
            for wi in range(W):
                value, a, b = walks[wi]
                value, a, b = step(value, a, b, x)
                walks[wi] = [value, a, b]
                total_steps += 1
                if b != 0:
                    b_nonzero_step_count += 1
                if is_distinguished(value):
                    if value in dp_table:
                        a0, b0 = dp_table[value]
                        denom = (b0 - b) % n
                        if denom == 0:
                            v2, a2, b2 = escape(value, a, b, x)
                            walks[wi] = [v2, a2, b2]
                            continue
                        cand = ((a - a0) * pow(denom, n - 2, n)) % n
                        recovered = (cand == x % n) or ((n - cand) % n == x % n)
                        break
                    else:
                        dp_table[value] = (a, b)
            if recovered is not None:
                break
        if recovered is True:
            correct += 1
        elif recovered is False:
            wrong += 1
        else:
            timeout += 1
    return correct, wrong, timeout, total_steps, b_nonzero_step_count


if __name__ == "__main__":
    params = json.load(sys.stdin)
    correct, wrong, timeout, steps, bnz = run(
        params["order"], params["w"], params["unit"], params["W"],
        params["dp_mod"], params["trials"], params.get("seed", 1))
    print(json.dumps({
        "group": "integers modulo n under addition, with multiply-by-unit automorphism",
        "order": params["order"], "w": params["w"], "W": params["W"],
        "dp_mod": params["dp_mod"], "trials": params["trials"],
        "correct": correct, "wrong_answers": wrong, "timeouts": timeout,
        "total_steps": steps, "b_nonzero_step_count": bnz,
        "all_correct_no_wrong_answers": wrong == 0,
        "scope": "toy group only; tests the mixed-branch/canonicalization/"
                 "distinguished-point/escape/linear-solve logic exactly as "
                 "used in production, not production curve arithmetic or "
                 "instruction-level cost bounds",
    }))
