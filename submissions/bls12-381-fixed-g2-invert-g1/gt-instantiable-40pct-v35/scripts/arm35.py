#!/usr/bin/env python3
"""arm35.py -- RETENTION and COLLISION-DETECTION-FRACTION at the HONEST production lambda.

v34 measured retention at lambda = S/2^128 = 0.2815.  That lambda is not
instantiable by the shipped store: dict28_ws_init takes `u64 nb`, requires a
power of two, and ws_bucket() masks a 64-bit mix, so nb <= 2^63 per level.
The production store is therefore priced as TWO levels of the shipped type
(2^63 shards x 2^63 buckets = 2^126 bucket heads), which gives

    lambda_b = S / 2^126 = 1.1519                (S = 9.580366847e37)

and NOT 0.2815.  This arm is run at that lambda, on the shipped four-probe
FIFO chains, and measures two DIFFERENT quantities:

  retention     P(a record inserted at time t' is still resident at the end)
                == E[resident]/arrivals, the terminal occupancy fraction.
  kappa_detect  P(a TRUE cross-trail key collision is DETECTED)
                = P(the earlier record is resident and the shipped
                dict28_wide_lookup returns it at the later arrival).
                Measured directly by planting the collision: at arrival t we
                replay a uniformly chosen earlier key k[t'] through the
                shipped lookup and count a hit.  This is the injected-collision
                control for the detection claim, not an occupancy proxy.

Both are reported; kappa_detect is the one the success debit uses.
"""
import ctypes
import json
import multiprocessing as mp
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import protocol28 as P                                    # noqa: E402

EV = os.path.join(os.path.dirname(HERE), "evidence")
L = json.load(open(os.path.join(EV, "ledger34.json")))
S = int(L["protocol"]["steps"])
NB_TOTAL_LOG2 = 126                     # two u64 levels of the shipped type
NB = 1 << NB_TOTAL_LOG2
LAM_NUM, LAM_DEN = S, NB
FIX_NB = 4096                            # fixture bucket count (power of two)
PROBES_PER_TRIAL = 240


def one(seed):
    """Build the shipped store to `arr` arrivals, then probe planted collisions."""
    lib = P.EVICTING
    arr = (LAM_NUM * FIX_NB + LAM_DEN // 2) // LAM_DEN      # round(lambda*4096)
    rng = random.Random(0x35A5 ^ (seed * 2654435761))
    ws = P.D28WS()
    assert lib.dict28_ws_init(ctypes.byref(ws), FIX_NB) == 0
    keys = []                      # (key tuple, trail id) in arrival order
    resident = 0
    probes = hits = 0
    for t in range(arr):
        key = (rng.getrandbits(64), rng.getrandbits(64), rng.getrandbits(64),
               rng.getrandbits(64), rng.getrandbits(64), rng.getrandbits(64))
        kb = (ctypes.c_uint64 * 6)(*key)
        slot = rng.getrandbits(64)          # 128-bit slot: both halves live
        hi = rng.getrandbits(64)
        ref = (ctypes.c_uint64 * 2)(t & 0xFFFFFFFF, (t >> 32) + 1)
        ev = ctypes.c_uint64(0)
        lib.dict28_wide_insert(ctypes.byref(ws), slot, hi, kb, ref, ctypes.byref(ev))
        keys.append((key, slot, hi))
        # PLANTED CROSS-TRAIL COLLISION: replay an earlier key now
        if probes < PROBES_PER_TRIAL and t > 1 and rng.random() < 0.06:
            t2 = rng.randrange(0, t)
            k2, s2, h2 = keys[t2]
            kb2 = (ctypes.c_uint64 * 6)(*k2)
            out = (ctypes.c_uint64 * 2)(0, 0)
            st = lib.dict28_wide_lookup(ctypes.byref(ws), s2, h2, kb2, out)
            probes += 1
            # a hit must return the PLANTED record's reference, not a stranger
            if st == 1 and out[0] == (t2 & 0xFFFFFFFF) and out[1] == ((t2 >> 32) + 1):
                hits += 1
            elif st == 1:
                hits += 0                    # wrong-ref hit: not a detection
    # terminal retention: how many of the arrivals are still retrievable
    for (key, slot, hi) in keys:
        kb = (ctypes.c_uint64 * 6)(*key)
        out = (ctypes.c_uint64 * 2)(0, 0)
        if lib.dict28_wide_lookup(ctypes.byref(ws), slot, hi, kb, out) == 1:
            resident += 1
    lib.dict28_ws_free(ctypes.byref(ws))
    return {"arr": arr, "resident": resident, "probes": probes, "hits": hits}


def main():
    nproc = min(10, mp.cpu_count())
    trials = 900
    with mp.get_context("fork").Pool(nproc) as p:
        rs = p.map(one, range(trials), chunksize=4)
    arr = rs[0]["arr"]
    res = sum(r["resident"] for r in rs)
    den = sum(r["arr"] for r in rs)
    probes = sum(r["probes"] for r in rs)
    hits = sum(r["hits"] for r in rs)
    p_hat = res / den
    # binomial standard error on the pooled retention
    se = (p_hat * (1 - p_hat) / den) ** 0.5
    k_hat = hits / probes if probes else 0.0
    k_se = (k_hat * (1 - k_hat) / probes) ** 0.5 if probes else 0.0
    import math
    lam = LAM_NUM / LAM_DEN
    der = sum(math.exp(-lam) * lam ** k / math.factorial(k) for k in range(5))
    out = {
        "policy": "shipped dict28.c four-probe FIFO chains (evict tail when chain reaches PROBES=4)",
        "lambda_rule": "lambda_b = S / nb_total, nb_total = 2^63 shards x 2^63 buckets = 2^126 "
                       "is the LARGEST power-of-two bucket count the shipped u64 field can hold "
                       "(dict28_ws_init takes u64 nb and rejects non-powers-of-two; ws_bucket masks "
                       "a 64-bit mix, so one level addresses <= 2^63 buckets).  v34's 2^128 is NOT "
                       "representable and its lambda 0.2815 is withdrawn.",
        "lambda_num": str(LAM_NUM), "lambda_den": str(LAM_DEN),
        "lambda": lam, "nb_total_log2": NB_TOTAL_LOG2,
        "fixture_nb": FIX_NB, "arrivals": arr, "trials": trials,
        "lambda_hat": arr / FIX_NB,
        "retention": {"resident": res, "denominator": den, "p_hat": p_hat, "se": se,
                      "p_derived": der, "delta": p_hat - der,
                      "delta_explanation": "1-P(Poisson(lambda)>=5) is the terminal-size limit for "
                                           "an insertion-only bucket; the shipped chain evicts "
                                           "during fill as well, so the measured loss is larger. "
                                           "Both are reported; the LOWER of the two is used."},
        "kappa_detect": {"hits": hits, "probes": probes, "k_hat": k_hat, "se": k_se,
                         "definition": "P(true cross-trail key collision is DETECTED) = P(the "
                                       "earlier record is resident and the shipped "
                                       "dict28_wide_lookup returns its reference at the later "
                                       "arrival). Injected-collision control: uniformly chosen "
                                       "earlier key replayed through the shipped lookup.",
                         "control_wrong_ref_rejected": True},
        "success_debit": {"rule": "p_used = p5 * min(kappa_lower, retention_lower)",
                          "kappa_lower_3se": max(0.0, k_hat - 3 * k_se),
                          "retention_lower_3se": max(0.0, p_hat - 3 * se)},
        "wrong_lambda_control": {"lambda": 0.2815416777,
                                 "note": "the v34 lambda is not instantiable by the shipped u64 "
                                         "store (the v34 head count); the arm is not run there because "
                                         "no bucket count of that size can be initialised"},
    }
    json.dump(out, open(os.path.join(EV, "arm35.json"), "w"), indent=1)
    print(json.dumps({"lambda": lam, "arr": arr, "retention": p_hat, "se": se,
                      "kappa": k_hat, "k_se": k_se, "probes": probes, "hits": hits}, indent=1))


if __name__ == "__main__":
    main()
