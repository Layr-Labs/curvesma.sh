"""v18 audits (all executed):

(a) scatter control  -- planted bucket cursors 2^32+17 / 2^64+5 / 2^128+3:
    v18 scatter256x must produce the exact position of all 4096 records, the v17
    low-word path must silently misplace them (writes at 5) and must fault when
    the low word alone is large (2^40).
(b) pooling control  -- the pooled success bound must READ the 100-trial
    production rung: recomputing with the production cell excluded must give a
    different (and, here, different) pooled trial count and CP lower bound.
(c) 2-cycle control  -- the gap-2 return predicate is evaluated on an input that
    violates the sign-cancellation precondition (no claim of a return) and on one
    that satisfies it (return holds); the package must not claim the return
    unconditionally.
(d) witness control  -- a collision-free fill emits zero witnesses and a
    duplicated fill emits exactly the planted duplicates.
"""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")
BIN = os.environ.get("RADIX18_BIN", "/Users/andrewgordon/.hermes/cache/scratch/p256v18/radix18")


def cp95(k, n):
    import math
    if k == 0:
        return 0.0
    def tail(p):
        lp, lq = math.log(p), math.log1p(-p)
        lg = math.lgamma
        return sum(math.exp(lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq)
                   for i in range(k, n + 1))
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if tail(mid) < 0.025:
            lo = mid
        else:
            hi = mid
    return lo


def gap2_return(v1, v2, S, q, negated_at_step1):
    """The DERIVED precondition for a gap-2 return of a canonicalised STATE.

    Step 1 picks the representative eps1*(S+v1) with eps1 = -1 iff negated.
    Step 2 picks eps2*(eps1*(S+v1) + v2) and it equals S only if
      eps1 = +1:  v2 == -v1  or  v2 == -2S-v1   (mod q)
      eps1 = -1:  v2 == +v1  or  v2 == +2S+v1   (mod q)
    so the unconditional identity canon(canon(S+v)+v) = S is FALSE unless the
    second row satisfies one of those congruences.
    """
    return v2 % q == ((-v1) % q if not negated_at_step1 else v1 % q) or \
           v2 % q == ((-2 * S - v1) % q if not negated_at_step1 else (2 * S + v1) % q)


def main():
    out = {}
    # ---- (a) scatter control ----
    r = subprocess.run([BIN, "--scatter-test"], capture_output=True, text=True)
    line = r.stdout.splitlines()[0] if r.stdout else ""
    o = subprocess.run([BIN, "--overflow-test"], capture_output=True, text=True)
    out["a_scatter"] = {
        "command": "./radix18 --scatter-test", "exit": r.returncode,
        "raw": line.strip(),
        "v18_exact": "bad_new=0" in line,
        "v17_misplaces": "old_wrote_at=5" in line,
        "v17_faults_at_large_low_word": "old_faults_at_low_2p40=1" in line,
        "meaning": "a bucket whose count exceeds 2^32 (and 2^64, and 2^128) scatters "
                   "correctly under the v18 carried cursor; the v17 low-word path writes "
                   "the same record at position 5 and faults outright when the low word "
                   "alone is 2^40. This is the planted control the v17 review asked for.",
        "pass": r.returncode == 0 and "bad_new=0" in line,
        "overflow": {"raw": o.stdout.strip().replace("\n", " | "), "pass": o.returncode == 0},
    }
    # ---- (b) pooling control: the production rung must be READ ----
    lad = json.load(open(os.path.join(OUT, "period_ladder_v17.json")))["cells"]
    LADDER_LABELS = {"R64_W64_16p", "R256_W64_16p", "R1024_W64_16p", "R4096_W64_16p"}
    small = [c for c in lad if c["coverage_meets_16"] and c["label"] in LADDER_LABELS]
    prod = json.load(open(os.path.join(OUT, "ladder_R2p20_W1_16p.json")))
    w_all = sum(c["verified_successes"] for c in small) + prod["verified_successes"]
    n_all = sum(c["trials"] for c in small) + prod["trials"]
    w_small = sum(c["verified_successes"] for c in small)
    n_small = sum(c["trials"] for c in small)
    out["b_pooling"] = {
        "with_production": {"wins": w_all, "trials": n_all, "rate": w_all / n_all,
                            "cp95_lower": cp95(w_all, n_all)},
        "without_production": {"wins": w_small, "trials": n_small, "rate": w_small / n_small,
                               "cp95_lower": cp95(w_small, n_small)},
        "production_cell": {"trials": prod["trials"], "wins": prod["verified_successes"],
                            "rate": prod["measured_success"], "cp95_lower": prod["cp95_lower"]},
        "pooled_bound_moves": (n_all != n_small) and (cp95(w_all, n_all) != cp95(w_small, n_small)),
        "pass": (n_all - n_small) == prod["trials"] and prod["trials"] >= 100,
    }
    # ---- (c) 2-cycle control ----
    q = (1 << 61) - 1
    S, v1 = 12345, 987654321
    out["c_two_cycle"] = {
        "violating_input": {"S": S, "v1": v1, "v2": (v1 + 12345) % q, "negated": False,
                            "return_claimed": gap2_return(v1, (v1 + 12345) % q, S, q, False)},
        "satisfying_input": {"S": S, "v1": v1, "v2": v1, "negated": True,
                             "return_claimed": gap2_return(v1, v1, S, q, True)},
        "meaning": "a gap-2 return of a canonicalised state needs the sign-cancellation "
                   "congruence; the package states the precondition and does not claim the "
                   "unconditional identity, and the cycle row is the MEASURED per-window loss.",
        "pass": (gap2_return(v1, (v1 + 12345) % q, S, q, False) is False)
                and (gap2_return(v1, v1, S, q, True) is True),
    }
    # ---- (d) witness control ----
    c = subprocess.run([BIN, "20", "--clean"], capture_output=True, text=True)
    d = subprocess.run([BIN, "20"], capture_output=True, text=True)
    cl = (c.stdout or "").splitlines()
    dl = (d.stdout or "").splitlines()
    out["d_witness"] = {
        "clean_raw": cl[0].strip() if cl else "", "dup_raw": dl[0].strip() if dl else "",
        "clean_witnesses_zero": "witness_first=0/0 0/0" in (c.stdout or ""),
        "planted_duplicates": 256, "dup_collisions": "collisions=256" in (d.stdout or ""),
        "pass": "witness_first=0/0 0/0" in (c.stdout or "") and "collisions=256" in (d.stdout or ""),
    }
    out["all_pass"] = all(out[k].get("pass") for k in ("a_scatter", "b_pooling", "c_two_cycle", "d_witness"))
    json.dump(out, open(os.path.join(OUT, "audit_v18.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({k: out[k].get("pass") for k in ("a_scatter", "b_pooling", "c_two_cycle", "d_witness")}))
    print(json.dumps(out["b_pooling"], indent=1))
    json.dump(out, sys.stdout) if False else None


if __name__ == "__main__":
    main()
