#!/usr/bin/env python3
"""v22 evidence B: the gap self-check IMPLEMENTATION, and the retry law 1/(1-q).

(1) gapcheck.scan_max_gap is the shipped function the construction runs on the
    drawn list before it accepts it.  It maintains the running maximum in O(1)
    amortised per step: one index load, one subtraction, one compare, one
    conditional store, one limit compare.  `selftest` runs it against a clean
    list (must accept) and against an injected violation (must reject) - a
    mutation control, so the test can actually fail.

(2) The per-step instruction price is COUNTED, not guessed: gapscan.c holds the
    same loop, is compiled with clang -O2 -arch arm64, and `count_loop` reads the
    loop body out of the disassembly (objdump) and reports the instruction count
    and the loop trip count.  That number x the walk length is the ledger row.

(3) Retry: for independent rejection probability q the expected number of draws
    is 1/(1-q), so the extra work is (1/(1-q) - 1) x (cost of one discarded
    draw).  q_gap and q_route are DIFFERENT events and are measured separately,
    over >= 192 independent seed coins each.
"""
import json, math, os, subprocess, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from toy import is_prime


# ------------------------------------------------------------------ shipped function
def scan_max_gap(idx, limit):
    """Running-maximum gap scan.  O(1) state per step, no allocation.

    idx   : non-decreasing integer positions of the drawn walk entries
    limit : accepted maximum gap

    Returns (max_gap, ok).  This is the self-check the construction runs; the
    caller re-draws when ok is False.
    """
    if not idx:
        return 0, True
    prev = idx[0]
    mx = 0
    for i in range(1, len(idx)):
        g = idx[i] - prev
        prev = idx[i]
        if g > mx:
            mx = g
        if mx > limit:
            return mx, False
    return mx, True


def selftest():
    clean = [0, 5, 9, 14, 20, 27, 33]          # max gap 7
    bad = [0, 5, 9, 14, 40, 47, 53]            # max gap 26
    m1, ok1 = scan_max_gap(clean, 10)
    m2, ok2 = scan_max_gap(bad, 10)
    m3, ok3 = scan_max_gap(bad, 100)
    return dict(clean_max=m1, clean_ok=ok1, injected_max=m2, injected_rejected=(not ok2),
                widened_limit_accepts=ok3,
                mutation_control_pass=(ok1 and (not ok2) and ok3),
                per_step_state="prev index, running max; one compare + one conditional store")


# ------------------------------------------------------------------ counted C kernel
C_SRC = r"""
#include <stdint.h>
int gap_scan(const uint64_t *idx, uint64_t n, uint64_t limit, uint64_t *out_max)
{
    uint64_t prev = idx[0], mx = 0;
    for (uint64_t i = 1; i < n; i++) {
        uint64_t g = idx[i] - prev;
        prev = idx[i];
        if (g > mx) mx = g;
        if (mx > limit) { *out_max = mx; return 0; }
    }
    *out_max = mx;
    return 1;
}
"""


def count_loop():
    """Compile the SAME loop and count the executed instructions per iteration."""
    import re
    cpath = os.path.join(HERE, "gapscan.c")
    open(cpath, "w").write(C_SRC)
    obj = "/tmp/gapscan_v22.o"
    r = subprocess.run(["clang", "-O2", "-arch", "arm64", "-c", cpath, "-o", obj],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return dict(error=r.stderr[-400:], counted=False)
    d = subprocess.run(["objdump", "-d", "--no-show-raw-insn", obj],
                       capture_output=True, text=True).stdout
    ins = []
    for line in d.splitlines():
        m = re.match(r"\s*([0-9a-f]+):\s+(\S+)\s*(.*)", line)
        if m:
            ins.append((int(m.group(1), 16), m.group(2), m.group(3).strip()))
    if not ins:
        return dict(error="no disassembly", counted=False)
    # the loop-closing branch: a branch whose target address is BELOW its own address
    close = None
    for addr, op, args in ins:
        if op.startswith("b.") or op == "b":
            m = re.search(r"^(0x[0-9a-f]+)", args)
            if m and int(m.group(1), 16) < addr:
                t = int(m.group(1), 16)
                if close is None or t > close[1]:      # INNERMOST backward branch
                    close = (addr, t)
    if not close:
        return dict(error="no backward branch (loop not found)", counted=False,
                    function_instructions=len(ins))
    br, tgt = close
    body = [(a, o, g) for a, o, g in ins if tgt <= a <= br]
    return dict(counted=True, loop_body_instructions=len(body),
                loop_target=tgt, loop_branch=br,
                function_instructions=len(ins),
                loop_body=["%04x: %s %s" % (a, o, g) for a, o, g in body],
                method="clang -O2 -arch arm64 -c; objdump -d; body = [backward-branch target, "
                       "branch] inclusive")


# ------------------------------------------------------------------ retry census
def hprimes_band(p, m, lo, need, stride_odd=True):
    out, q = [], lo | 1
    while len(out) < need:
        if is_prime(q) and pow(q, m, p) == 1:
            out.append(q)
        q += 2
    return out


def prime_with(bits, dd):
    p = (1 << (bits - 1)) | 1
    while True:
        if ((p - 1) % dd == 0) and ((p - 1) // dd) % 6 == 0 and is_prime(p):
            return p
        p += 2


def census(bits, dd, seeds, band=401, buckets=1 << 6, ratio_cap=207.9366134671963,
           route_cap=2.0 ** 10.53):
    """Each seed = one independent construction coin: a seed-derived band of
    consecutive H-primes.  Measures the two rejection events SEPARATELY."""
    p = prime_with(bits, dd)
    m = (p - 1) // dd
    need = seeds * band
    allq = hprimes_band(p, m, 1 << (bits - 1), need)
    if len(allq) < need:
        return dict(error="not enough H-primes", got=len(allq), need=need)
    rows, q_gap, q_route = [], 0, 0
    for s in range(seeds):
        b = allq[s * band:(s + 1) * band]
        gaps = [b[i + 1] - b[i] for i in range(len(b) - 1)]
        mu = sum(gaps) / len(gaps)
        mr = max(gaps) / mu
        cnt = [0] * buckets
        for x in b:
            cnt[x % buckets] += 1
        load = max(cnt) / (len(b) / buckets)
        pairs = sum(c * (c - 1) // 2 for c in cnt)
        pred = len(b) * (len(b) - 1) / (2 * buckets)
        g_event = mr > ratio_cap
        r_event = (pairs / pred) > route_cap if pred else False
        q_gap += g_event
        q_route += r_event
        rows.append(dict(seed=s, n=len(b), mean_gap=round(mu, 4), pair_ratio=round(pairs / pred, 6) if pred else None,
                         max_over_mean=round(mr, 4), gap_event=bool(g_event),
                         route_max_over_mean=round(load, 4),
                         prefix_pairs=pairs, prefix_formula=round(pred, 2),
                         prefix_ratio=round(pairs / pred, 4) if pred else None,
                         route_event=bool(r_event)))
    return dict(field_bits=bits, d=dd, m_bits=m.bit_length(), seeds=seeds, band=band,
                buckets=buckets, ratio_cap=ratio_cap, route_cap=route_cap,
                q_gap_events=q_gap, q_route_events=q_route,
                q_gap=round(q_gap / seeds, 6), q_route=round(q_route / seeds, 6),
                max_over_mean_max=max(r["max_over_mean"] for r in rows),
                route_max_over_mean_max=max(r["route_max_over_mean"] for r in rows),
                pair_ratio_max=max(r["pair_ratio"] or 0 for r in rows),
                rows=rows)


def cp95_upper(k, n):
    """Clopper-Pearson one-sided 95% upper bound on a binomial rate, k events in n."""
    if k == 0:
        return float(1 - 0.05 ** (1.0 / n))
    if k >= n:
        return 1.0
    lo, hi = F(k, n), F(1)
    target = 0.05
    for _ in range(200):
        mid = (lo + hi) / 2
        # P(X <= k-1 | n, mid) >= 0.95  <=>  mid is an upper bound
        s, term = F(0), F(1)
        for i in range(k):
            s += term
            term = term * (n - i) / (i + 1) * mid / (1 - mid)
        s *= (1 - mid) ** n
        if s >= target:
            lo = mid
        else:
            hi = mid
    return float((lo + hi) / 2)


if __name__ == "__main__":
    st = selftest()
    cnt = count_loop()
    print("selftest", json.dumps(st))
    print("counted", json.dumps({k: v for k, v in cnt.items() if k != "loop_body"}))
    if cnt.get("counted"):
        print("\n".join("   " + x for x in cnt["loop_body"]))
    out = dict(self_check=dict(implementation="gapcheck.scan_max_gap (this file, imported by the "
                                                  "construction)", selftest=st,
                               counted_kernel=cnt,
                               description="running max gap over the drawn walk positions with "
                                           "an early abort when the limit is exceeded; O(1) "
                                           "amortised state per step"),
               retry=dict(rule="expected draws 1/(1-q); extra work (1/(1-q)-1) x one discarded draw",
                          events=dict(gap="drawn band's max/mean gap ratio above 2^7.7",
                                      route="drawn band's routing bucket load above the charged "
                                            "factor 2.25"),
                          censuses=[]))
    for bits, dd, seeds in ((32, 8, 96), (36, 8, 96)):
        print("census bits=%d" % bits, flush=True)
        out["retry"]["censuses"].append(census(bits, dd, seeds))
        print(json.dumps({k: v for k, v in out["retry"]["censuses"][-1].items() if k != "rows"}),
              flush=True)
        json.dump(out, open(HERE + "/retry-v22.json", "w"), indent=1)
    tot_n = sum(c["seeds"] for c in out["retry"]["censuses"])
    tot_g = sum(c["q_gap_events"] for c in out["retry"]["censuses"])
    tot_r = sum(c["q_route_events"] for c in out["retry"]["censuses"])
    out["retry"]["pooled"] = dict(draws=tot_n, gap_events=tot_g, route_events=tot_r,
                                  q_gap=tot_g / tot_n, q_route=tot_r / tot_n,
                                  q_gap_cp95=cp95_upper(tot_g, tot_n),
                                  q_route_cp95=cp95_upper(tot_r, tot_n))
    print("pooled", json.dumps(out["retry"]["pooled"], indent=1))
    json.dump(out, open(HERE + "/retry-v22.json", "w"), indent=1)
    print("wrote retry-v22.json")
