#!/usr/bin/env python3
"""probe_v15n.py -- executes the SHIPPED alignment object (scripts/align_step_v15n.c,
compiled clang -O2 -arch arm64) through ctypes, on genuine secp256k1 points.

Arms (all counts emitted to evidence/equivariance-v15n.json):
  arith        fmul/fsub of the shipped object vs Python bigint           (fixed)
  class        class_slot(key,slot) of the shipped object vs bigint min   (fixed)
  orbit        six orbit partners of 600 random points: the shipped
               class_slot must put every partner in ONE class (0 splits)  (fixed)
  raw_control  the previously shipped rule (index = raw x & 2047, no
               alignment): partners must split (mutation control)
Nothing here re-implements the alignment: every decision comes from the .dylib.
"""
import ctypes, json, os, random, subprocess, sys, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
BETA = 0x7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE
G = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
     0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)
U64 = ctypes.c_uint64


def build():
    out = os.path.join(EV, "libalign_v15n.dylib")
    subprocess.run(["clang", "-O2", "-arch", "arm64", "-fPIC", "-shared", "-o", out,
                    os.path.join(HERE, "align_step_v15n.c"),
                    os.path.join(HERE, "driver_align_v15n.c")], check=True)
    return out


def arr(v):
    a = (U64 * 4)()
    for i in range(4):
        a[i] = (v >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
    return a


def gv(a):
    return sum(int(a[i]) << (64 * i) for i in range(4))


def add(Pt, Qt):
    if Pt is None: return Qt
    if Qt is None: return Pt
    x1, y1 = Pt; x2, y2 = Qt
    if x1 == x2:
        if (y1 + y2) % P == 0: return None
        lam = 3 * x1 * x1 * pow(2 * y1, -1, P) % P
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, P) % P
    x3 = (lam * lam - x1 - x2) % P
    return (x3, (lam * (x1 - x3) - y1) % P)


def main():
    lib = ctypes.CDLL(build())
    N = 600
    random.seed(20261024)

    bad = 0
    for _ in range(200):
        a_, b_ = random.randrange(P), random.randrange(P)
        r = (U64 * 4)()
        lib.api_fmul(r, arr(a_), arr(b_))
        if gv(r) != (a_ * b_) % P: bad += 1
        s = (U64 * 4)()
        lib.api_fsub(s, arr(a_), arr(b_))
        if gv(s) != (a_ - b_) % P: bad += 1
    arith = {"cases": 400, "mismatches": bad}

    bad = 0
    for _ in range(N):
        x = random.randrange(1, P)
        key = (U64 * 4)(); slot = ctypes.c_uint64()
        lib.api_class(arr(x), key, ctypes.byref(slot))
        vals = [x, x * BETA % P, x * BETA % P * BETA % P]
        i = min(range(3), key=lambda k: vals[k])
        if gv(key) != vals[i] or slot.value != i: bad += 1
    klass = {"cases": N, "mismatches": bad}

    # orbit arm: partners sigma^k S, k = 0,1,2 (beta) times sign
    pts = []
    for _ in range(N):
        k = random.randrange(1, 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141)
        pts.append(k)
    orbit_split = 0
    raw_split = 0
    for k in pts:
        Q = mul(k, G)
        base_key = None; base_raw = None
        for j in range(3):
            x = Q[0]
            for _ in range(j):
                x = x * BETA % P
            for sgn in (0, 1):
                y = Q[1] if sgn == 0 else P - Q[1]
                key = (U64 * 4)(); slot = ctypes.c_uint64()
                lib.api_class(arr(x), key, ctypes.byref(slot))
                kk = gv(key)
                if base_key is None: base_key = kk
                elif kk != base_key: orbit_split += 1
                raw = x & 2047
                if base_raw is None: base_raw = raw
                elif raw != base_raw: raw_split += 1
    res = {"arith": arith, "class": klass,
           "orbit_partners": 6 * N, "orbit_pairs_checked": 5 * N,
           "orbit_splits": orbit_split,
           "raw_control_pairs": 5 * N, "raw_control_expected_splits": 4 * N, "raw_control_splits": raw_split,
           "object": "scripts/align_step_v15n.c (clang -O2 -arch arm64) via ctypes",
           "sha256_lib": hashlib.sha256(open(os.path.join(EV, "libalign_v15n.dylib"), "rb").read()).hexdigest()}
    open(os.path.join(EV, "equivariance-v15n.json"), "w").write(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
    ok = (arith["mismatches"] == 0 and klass["mismatches"] == 0
          and orbit_split == 0 and raw_split == 4 * N)  # negation keeps x: 4 of 5 partners split
    print("PROBE", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def mul(k, Pt):
    R = None
    while k:
        if k & 1: R = add(R, Pt)
        Pt = add(Pt, Pt); k >>= 1
    return R


if __name__ == "__main__":
    sys.exit(main())
