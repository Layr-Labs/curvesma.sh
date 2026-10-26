import importlib.util
import json
import math
import multiprocessing
import os
import random
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "walk_ref", os.path.join(HERE, "..", "measurement", "walk_ref.py"))
W = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(W)
p, G = W.p, W.G

TASKS = [
    ("tau", 3, 200000, 703),
    ("tau", 6, 400000, 701),
    ("tau", 6, 400000, 702),
    ("fresh", 2, 60000, 31),
    ("fresh", 3, 60000, 32),
    ("fresh", 3, 60000, 101),
    ("fresh", 3, 60000, 102),
    ("fresh", 3, 60000, 103),
    ("fresh", 2, 60000, 104),
    ("fresh", 2, 60000, 105),
]


def neg(P):
    return None if P is None else (P[0], (-P[1]) % p)


def parity(P):
    return W.M(P[1]) & 1


def key(P):
    return W.M(P[0])


def canon(P):
    return P if parity(P) == 0 else neg(P)


def cells():
    return {"00": 0, "01": 0, "10": 0, "11": 0}


def tau_run(A_bits, steps, seed):
    rnd = random.Random(seed)
    A = 1 << A_bits
    R = [W.mul(rnd.randrange(1, W.r), G) for _ in range(A)]
    Wk = W.mul(rnd.randrange(1, W.r), G)
    all_pairs, p2_pairs, art_pairs = cells(), cells(), cells()
    eq_pairs, neq_pairs = cells(), cells()
    seg_canon, seg_t, seg_idx, seen = [], [], [], {}
    n_escape = 0

    def finish_segment(last_t_revisit):
        for q in range(len(seg_t) - 1):
            cell = f"{seg_t[q]}{seg_t[q + 1]}"
            all_pairs[cell] += 1
            if seg_idx[q] == seg_idx[q + 1]:
                eq_pairs[cell] += 1
            else:
                neq_pairs[cell] += 1
            if q + 2 < len(seg_t):
                p2_pairs[cell] += 1
        if last_t_revisit is not None and seg_t:
            art_pairs[f"{seg_t[-1]}{last_t_revisit}"] += 1

    for _ in range(steps):
        i = key(Wk) & (A - 1)
        s = parity(Wk)
        c = canon(Wk)
        Pp = W.add(c, R[i])
        t = parity(Pp)
        if c[0] in seen:
            finish_segment(t)
            j = seen[c[0]]
            M = min(seg_canon[j:], key=lambda P: key(P))
            Wk = W.add(M, M)
            seg_canon, seg_t, seg_idx, seen = [], [], [], {}
            n_escape += 1
            continue
        seen[c[0]] = len(seg_canon)
        seg_canon.append(c)
        seg_t.append(t)
        seg_idx.append(i)
        Wk = W.add(Wk, R[i] if s == 0 else neg(R[i]))
    return dict(A=A, steps=steps, seed=seed, escapes=n_escape, all=all_pairs, p2=p2_pairs,
                artifact_incl_revisit_tau=art_pairs, eq_index=eq_pairs, neq_index=neq_pairs)


def fresh_run(A_bits, steps, seed):
    rnd = random.Random(seed)
    A = 1 << A_bits
    R = [W.mul(rnd.randrange(1, W.r), G) for _ in range(A)]
    Wk = W.mul(rnd.randrange(1, W.r), G)
    seg_canon, seg_idx, seg_t, seen = [], [], [], {}
    arr = []
    flips_fresh = [0, 0]
    pairs = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0}
    for _ in range(steps):
        i = key(Wk) & (A - 1)
        s = parity(Wk)
        c = canon(Wk)
        if c[0] in seen:
            j = seen[c[0]]
            M = min(seg_canon[j:], key=lambda P: key(P))
            Wk = W.add(M, M)
            seg_canon, seg_idx, seg_t, seen = [], [], [], {}
            continue
        m = len(seg_canon)
        Pp = W.add(c, R[i])
        t = parity(Pp)
        fresh_arrival = canon(Pp)[0] not in seen
        seen[c[0]] = m
        seg_canon.append(c)
        seg_idx.append(i)
        seg_t.append(t)
        if fresh_arrival:
            flips_fresh[t] += 1
        arr.append((i, t, fresh_arrival))
        Wk = W.add(Wk, R[i] if s == 0 else neg(R[i]))
    for q in range(len(arr) - 1):
        (i0, t0, f0), (i1, t1, f1) = arr[q], arr[q + 1]
        if f0 and f1 and i0 == i1:
            pairs[(t0, t1)] += 1
    return dict(A=A, steps=steps, seed=seed, fresh_arrival_flips=flips_fresh,
                pairs={f"{a}{b}": v for (a, b), v in pairs.items()})


def run_task(task):
    kind, A_bits, steps, seed = task
    return tau_run(A_bits, steps, seed) if kind == "tau" else fresh_run(A_bits, steps, seed)


def zscore(k, n, q):
    return (k - n * q) / math.sqrt(n * q * (1 - q))


def share_block(c, pred=None):
    n = c["10"] + c["11"]
    k = c["11"]
    out = {"n": n, "minus_then_minus": k, "share": round(k / n, 6), "z_vs_half": round(zscore(k, n, 0.5), 2)}
    if pred is not None:
        out["predicted"] = pred
        out["z_vs_predicted"] = round(zscore(k, n, pred), 2)
    return out


def poisson_pmf(k, lam):
    return math.exp(-lam + k * math.log(lam) - math.lgamma(k + 1))


def poisson_cdf(k, lam):
    return sum(poisson_pmf(i, lam) for i in range(k + 1))


def poisson_two_sided(k, lam):
    return min(1.0, 2 * min(poisson_cdf(k, lam), 1 - poisson_cdf(k - 1, lam)))


def production_steps():
    seed1 = json.load(open(os.path.join(HERE, "measurement", "seed1.json"), encoding="utf-8"))
    names = ("input_set_1", "input_set_2")
    c63 = [seed1[n]["events"]["batch_steps_63"]["two_cycles"] for n in names]
    c64 = [seed1[n]["events"]["batch_steps_64"]["two_cycles"] for n in names]
    walks, A = 65536, 65536
    rate = 1 / (2 * A)
    exp63, exp64 = 63 * walks * rate, 64 * walks * rate
    inc = [c64[i] - c63[i] for i in range(2)]
    return {
        "purpose": "Two-cycle counts copied from measurement/seed1.json (production build, 2^16 adders, 65536 walks); "
                   "expected per input set = steps * walks / (2A), with A = 2^16.",
        "source": "evidence/measurement/seed1.json",
        "two_cycles": {names[i]: {"batch_steps_63": c63[i], "batch_steps_64": c64[i]} for i in range(2)},
        "expected_per_set": {"batch_steps_63": exp63, "batch_steps_64": exp64},
        "increment_64": {"input_set_1": inc[0], "input_set_2": inc[1], "pooled": sum(inc),
                         "expected_pooled": 2 * rate * walks,
                         "one_sided_p": round(1 - poisson_cdf(sum(inc) - 1, 1.0), 4)},
        "cumulative_63": {"pooled": sum(c63), "expected": 2 * exp63,
                          "two_sided_p": round(poisson_two_sided(sum(c63), 2 * exp63), 4)},
        "cumulative_64": {"pooled": sum(c64), "expected": 2 * exp64,
                          "two_sided_p": round(poisson_two_sided(sum(c64), 2 * exp64), 4)},
    }


if __name__ == "__main__":
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(min(len(TASKS), os.cpu_count() or 1)) as pool:
        results = pool.map(run_task, TASKS, chunksize=1)
    tau = [r for r, t in zip(results, TASKS) if t[0] == "tau"]
    fresh = [r for r, t in zip(results, TASKS) if t[0] == "fresh"]
    prod = production_steps()
    summary = {}
    tau_summary = {}
    for r in tau:
        q = 0.5 + 1 / (2 * r["A"])
        tau_summary[f"A{r['A']}_s{r['seed']}"] = {
            "all_pairs": share_block(r["all"], q),
            "equal_index": share_block(r["eq_index"]),
            "non_equal_index": share_block(r["neq_index"]),
        }
    summary["tau_runs"] = tau_summary
    a64 = [r for r in tau if r["A"] == 64]
    pooled64 = {k: {c: sum(r[k][c] for r in a64) for c in ("00", "01", "10", "11")}
                for k in ("all", "eq_index", "neq_index")}
    summary["A64_two_seeds"] = {
        "seeds": [r["seed"] for r in a64],
        "all_pairs": share_block(pooled64["all"], 0.5 + 1 / 128),
        "equal_index": share_block(pooled64["eq_index"]),
        "non_equal_index": share_block(pooled64["neq_index"]),
    }
    n_fresh = sum(sum(r["fresh_arrival_flips"]) for r in fresh)
    k_fresh = sum(r["fresh_arrival_flips"][1] for r in fresh)
    summary["fresh_arrival"] = {
        "runs": [{"A": r["A"], "seed": r["seed"], "n": sum(r["fresh_arrival_flips"]),
                  "share_tau_plus_1": round(r["fresh_arrival_flips"][1] / sum(r["fresh_arrival_flips"]), 4),
                  "z_vs_half": round(zscore(r["fresh_arrival_flips"][1], sum(r["fresh_arrival_flips"]), 0.5), 2)}
                 for r in fresh],
        "pooled": {"n": n_fresh, "share_tau_plus_1": round(k_fresh / n_fresh, 4),
                   "z_vs_half": round(zscore(k_fresh, n_fresh, 0.5), 2)},
    }
    out = {
        "purpose": "Sign-law replica counts for proof.md section 2.2. tau_transition_runs and fresh_arrival_runs are raw "
                   "counts from small builds (A = 8 and 64, production escape rule); production_steps copies the production "
                   "two-cycle counts; summary holds the derived shares, counts and z values. Written and checked by sign_law.py.",
        "tau_transition_runs": {f"A{r['A']}_s{r['seed']}": r for r in tau},
        "fresh_arrival_runs": fresh,
        "production_steps": prod,
        "summary": summary,
    }
    path = os.path.join(HERE, "sign-law.json")
    text = json.dumps(out, separators=(",", ":"), ensure_ascii=False) + "\n"
    stored = open(path, encoding="utf-8").read() if os.path.exists(path) else None
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    with open(path, encoding="utf-8") as f:
        back = f.read()
    assert back == text and json.loads(back) == json.loads(text)
    print("checked:", os.path.basename(path), len(text.encode("utf-8")), "bytes;",
          "unchanged from stored copy" if stored == text else ("no stored copy" if stored is None else "CHANGED from stored copy"))
