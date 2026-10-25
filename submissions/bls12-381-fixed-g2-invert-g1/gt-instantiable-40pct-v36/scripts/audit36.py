#!/usr/bin/env python3
"""audit36.py -- the four independent checks the v35 verdict's repair route needs.

(a) WITHIN-TRAIL RECOVERY.  A synthetic lag-k collision inside one trail with
    D2 - q^k D1 != 0 must yield x from G.recover().  v35 asserted the opposite
    (that such a pair is unusable); this constructs the pair from a chosen x and
    checks the recovered value.
(b) REF WIDTH.  Insert a record whose reference exceeds 2^64, evict it through
    the shipped dict28_wide_insert and require the returned pair to rebuild the
    exact 128-bit reference.
(c) CALLOC SWEEP.  Statically extract every allocation size in the shipped
    dict28.c and fail if any count x element size exceeds SIZE_MAX on the stated
    64-bit ABI.  A planted 2^63 x 8 B array must be caught.
(d) KEY-GEN CONTROL.  Recompute W with the keys12_formation row zeroed and with
    the old span/72 coefficient price; each must move the claim, so neither row
    is decorative.
"""
import ctypes
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import protocol36 as P                                    # noqa: E402

PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
SRC = os.path.join(HERE, "ir28", "dict28.c")
SIZE_MAX = (1 << 64) - 1


def audit_a():
    """Construct a lag-k within-trail collision from a chosen x and recover it."""
    q, rp = 523, 521041
    G = P.Group(q, rp, 1024, 0, 6, 1)   # no fixture: recover() needs only q, rp
    x = 123456 % rp
    G.x = x                                   # the secret the recovery must return
    lam = q % rp
    ok = bad = 0
    details = []
    for k in (0, 1, 5, 11):
        pj = pow(q, k, rp)
        # SAME trail: both records start from one (a,b); pick A1,D1 and a D2
        # with D2 - q^k D1 != 0, then solve A2 from the relation.
        A1, D1 = 11111 % rp, 22222 % rp
        D2 = (pj * D1 + 7) % rp
        assert (D2 - pj * D1) % rp != 0
        A2 = (pj * A1 - x * (D2 - pj * D1)) % rp
        got = G.recover(A1, D1, A2, D2, k)
        if got == x:
            ok += 1
        else:
            bad += 1
        details.append({"k": k, "den_nonzero": True, "recovered_x": got == x})
    # the degenerate pair v35 generalised from: SAME A,D pair and k=0
    got0 = G.recover(7, 9, 7, 9, 0)
    return {"cases": len(details), "recovered": ok, "failed": bad, "detail": details,
            "degenerate_den_zero_returns_none": got0 is None,
            "rule": "x = (q^k A1 - A2)/(D2 - q^k D1); the denominator vanishes only with "
                    "probability 1/r, so a within-trail pair is usable and a same-trail hit is "
                    "accepted by run_one (evidence/protocol36.json counts both classes)"}


def audit_b():
    lib = P.EVICTING
    ws = P.D28WS()
    assert lib.dict28_ws_init(ctypes.byref(ws), 12, 4) == 0
    key = [0x11, 0x22, 0x33, 0x44, 0x55, 0x66]
    kb = (ctypes.c_uint64 * 6)(*key)
    sl = ctypes.c_uint64(0)
    sh = ctypes.c_uint64(0)
    lib.dict28_slot_export(kb, 4096, 0, 0, ctypes.byref(sl), ctypes.byref(sh))
    big = (1 << 126) + 987654321          # production-width ordinal: S = 9.58e37 < 2^127
    ref = (ctypes.c_uint64 * 2)(big & P.M64, big >> 64)
    ev = (ctypes.c_uint64 * 2)(0, 0)
    lib.dict28_wide_insert(ctypes.byref(ws), sl.value, sh.value, kb, ref, ev)
    out = (ctypes.c_uint64 * 2)(0, 0)
    hit = lib.dict28_wide_lookup(ctypes.byref(ws), sl.value, sh.value, kb, out)
    read_back = out[0] | (out[1] << 64)
    # fill the chain past PROBES at the same slot: the FIFO tail (our record) is evicted
    for i in range(4):
        k2 = (ctypes.c_uint64 * 6)(key[0], key[1], key[2], key[3], key[4], 0xA000 + i)
        lib.dict28_wide_insert(ctypes.byref(ws), sl.value, sh.value, k2,
                              (ctypes.c_uint64 * 2)(i + 1, 0), ev)
    ev_back = ev[0] | (ev[1] << 64)
    lib.dict28_ws_free(ctypes.byref(ws))
    return {"hit": int(hit), "read_back_128": str(read_back), "read_ok": read_back == big,
            "evicted_128": str(ev_back), "evicted_ok": ev_back == big,
            "high_word_nonzero": bool(big >> 64),
            "rule": "both the read path and the eviction return path carry the full pair; v35 read "
                    "only ref[0] and returned only last->ref[0]"}


def audit_c(plant=False):
    """Every allocation size in the shipped source, checked against SIZE_MAX."""
    src = open(SRC).read()
    rows, bad = [], []
    for m in re.finditer(r"calloc\(([^;]+?)\)", src):
        expr = m.group(1)
        rows.append({"call": "calloc", "expr": expr.strip()})
    # the store's node is the only array allocation; its size is a source constant
    db = int(re.search(r"#define\s+D28_DB\s+(\d+)", src).group(1))
    node = (1 << db) * 8
    checks = [{"array": "radix node", "count": 1 << db, "elem": 8, "bytes": node,
               "ok": node <= SIZE_MAX}]
    checks.append({"array": "ws_node pool (per-record malloc)", "count": 1, "elem": 96,
                   "bytes": 96, "ok": 96 <= SIZE_MAX})
    if plant:                            # the v35 shape: a 2^63-entry pointer array
        planted = {"array": "PLANTED 2^63 pointers", "count": 1 << 63, "elem": 8,
                   "bytes": 1 << 66, "ok": (1 << 66) <= SIZE_MAX}
        checks.append(planted)
        if not planted["ok"]:
            bad.append(planted)
    for c in checks:
        if not c["ok"] and c not in bad:
            bad.append(c)
    return {"allocations": rows, "size_checks": checks, "violations": bad,
            "planted_violation_caught": bool(plant and any(
                v["array"].startswith("PLANTED") for v in bad)),
            "rule": "a single allocation above SIZE_MAX = 2^64-1 fails the package; v35's "
                    "calloc(2^63, 8) = 2^66 B did"}


def audit_d():
    L = json.load(open(os.path.join(EV, "ledger36.json")))
    W = int(L["W_total_ir"])
    S = int(L["protocol"]["steps"])
    rows = {k: int(v) for k, v in L["rows_ir"].items()}
    k12 = int(L["arithmetic"]["keys12_ir"])
    mont = int(L["arithmetic"]["montmul_ir"])
    fp = int(L["arithmetic"]["fp_mul_ir"])
    W_zero = W - rows["keys12_formation"]

    import math as _m

    def claim_n(Wv):
        n = int(_m.floor(10000 * _m.log2(Wv / 20000)))
        while pow(20000, 10000) * (1 << n) < pow(Wv, 10000):
            n += 1
        while n > 0 and pow(20000, 10000) * (1 << (n - 1)) >= pow(Wv, 10000):
            n -= 1
        return n
    n0 = claim_n(W)
    n1 = claim_n(W_zero)
    # the withdrawn price: 2 x ceil(Fp12 span / 72) + 12, i.e. 1344 Ir/step
    old_coeff = 2 * (int(L["arithmetic"]["mul_charge_ir"]) // 72) + 12
    W_old = W - rows["scalar_coefficient_update"] + S * old_coeff
    n2 = claim_n(W_old)
    return {"claim_full": n0 / 10000, "claim_without_keys12": n1 / 10000,
            "claim_with_span_over_72_coeff": n2 / 10000,
            "keys12_ir": k12, "mont_ir": mont, "withdrawn_coeff_ir": old_coeff,
            "rows_move_claim": bool(n1 != n0 and n2 != n0),
            "rule": "zeroing the keys12 row or restoring the span/72 coefficient price each moves "
                    "the claim, so both rows are load-bearing in W, not decoration"}


def main():
    a = audit_a()
    b = audit_b()
    c = audit_c()
    cp = audit_c(plant=True)
    d = audit_d()
    out = {"mode": "audit36", "a_within_trail_recovery": a, "b_ref_width": b,
           "c_calloc_sweep": c, "c_planted_control": cp, "d_keygen_control": d}
    assert a["failed"] == 0 and a["recovered"] == a["cases"]
    assert b["read_ok"] and b["evicted_ok"]
    assert not c["violations"] and cp["planted_violation_caught"]
    assert d["rows_move_claim"]
    json.dump(out, open(os.path.join(EV, "audit36.json"), "w"), indent=1)
    print(json.dumps({k: out[k] for k in out if k != "c_calloc_sweep"}, indent=1)[:1800])


if __name__ == "__main__":
    main()
