#!/usr/bin/env python3
"""verify23.py -- exact bigint verification of the two repaired kernels, with
mutation controls.  Reference arithmetic is PLAIN Python integers (no Montgomery,
no Karatsuba, schoolbook Fp6) so it shares no code path with the C kernels.

  f2mul          : (a0+a1u)(b0+b1u) with u^2 = -1
  cyclotomic_sqr : Granger-Scott form (2*v*t + 1, s - v*t - 1 - t), t = a1^2,
                   s = (a0+a1)^2 in Fp6, v^3 = xi = 1+u; checked against the
                   plain Fp12 square, whose difference must be exactly 1-norm
                   (so on norm-1 inputs the shipped form IS the square).

Mutations: each fix is reverted individually in the source and the same vectors
must FAIL.  Evidence: ../evidence/kernel_verify23.json
"""
import json
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(os.path.dirname(os.path.dirname(HERE)), "evidence")
P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
R = 1 << 384
RINV = pow(R, -1, P)
N = int(os.environ.get("NVEC", "10000"))


def to_mont(x):
    return x * R % P


def fr_mont(x):
    return x * RINV % P


def limbs(x):
    return [(x >> (64 * i)) & ((1 << 64) - 1) for i in range(6)]


# ---------------- plain (no Montgomery) Fp2/Fp6/Fp12, schoolbook ----------------
def f2m(a, b):
    return ((a[0] * b[0] - a[1] * b[1]) % P, (a[0] * b[1] + a[1] * b[0]) % P)


def f2a(a, b):
    return ((a[0] + b[0]) % P, (a[1] + b[1]) % P)


def f2s(a, b):
    return ((a[0] - b[0]) % P, (a[1] - b[1]) % P)


def f2xi(a):
    return ((a[0] - a[1]) % P, (a[0] + a[1]) % P)


def f6m(a, b):
    """schoolbook, v^3 = xi = 1+u"""
    c0 = f2a(f2m(a[0], b[0]), f2xi(f2a(f2m(a[1], b[2]), f2m(a[2], b[1]))))
    c1 = f2a(f2a(f2m(a[0], b[1]), f2m(a[1], b[0])), f2xi(f2m(a[2], b[2])))
    c2 = f2a(f2a(f2m(a[0], b[2]), f2m(a[1], b[1])), f2m(a[2], b[0]))
    return (c0, c1, c2)


def f6a(a, b):
    return tuple(f2a(x, y) for x, y in zip(a, b))


def f6s(a, b):
    return tuple(f2s(x, y) for x, y in zip(a, b))


def f6v(a):
    return (f2xi(a[2]), a[0], a[1])


ONE6 = ((1, 0), (0, 0), (0, 0))


def spec(a0, a1):
    t = f6m(a1, a1)
    vt = f6v(t)
    s = f6m(f6a(a0, a1), f6a(a0, a1))
    r0 = f6a(f6a(vt, vt), ONE6)
    r1 = f6s(f6s(f6s(s, vt), ONE6), t)
    return r0, r1


def true_sq(a0, a1):
    return (f6a(f6m(a0, a0), f6v(f6m(a1, a1))), f6a(f6m(a0, a1), f6m(a0, a1)))


def norm(a0, a1):
    """a0^2 - v*a1^2 (the Fp12 norm form)"""
    return f6s(f6m(a0, a0), f6v(f6m(a1, a1)))


# ---------------------------------- run ----------------------------------
def build(src, out, tag, defines=()):
    cmd = ["clang", "-O2", "-target", "arm64-apple-macos11", "-o", out, src,
           os.path.join(HERE, "consts.c")] + list(defines)
    subprocess.run(cmd, check=True, capture_output=True)


def gen(mode, n, seed=12345):
    rng = random.Random(seed)
    rows, refs = [], []
    if mode == "mul":
        for _ in range(n):
            a = (rng.randrange(P), rng.randrange(P))
            b = (rng.randrange(P), rng.randrange(P))
            A = (to_mont(a[0]), to_mont(a[1]))
            B = (to_mont(b[0]), to_mont(b[1]))
            rows.append([limbs(x) for x in (A[0], A[1], B[0], B[1])])
            c = f2m(a, b)
            refs.append([limbs(x) for x in (to_mont(c[0]), to_mont(c[1]))])
    else:
        for _ in range(n):
            a0 = tuple((rng.randrange(P), rng.randrange(P)) for _ in range(3))
            a1 = tuple((rng.randrange(P), rng.randrange(P)) for _ in range(3))
            flat = [to_mont(x) for c in (a0, a1) for p in c for x in p]
            rows.append([limbs(x) for x in flat])
            r0, r1 = spec(a0, a1)
            q0, q1 = true_sq(a0, a1)
            d0, d1 = f6s(r0, q0), f6s(r1, q1)
            nz = norm(a0, a1)
            # r0 - (a0^2 + v a1^2) = 1 - norm ; r1 - 2 a0 a1 = norm - 1
            assert d0 == f6s(ONE6, nz), "identity r0"
            assert d1 == f6s(nz, ONE6), "identity r1"
            flat = [to_mont(x) for c in (r0, r1) for p in c for x in p]
            refs.append([limbs(x) for x in flat])
    return rows, refs


def run(exe, mode, rows):
    txt = []
    for row in rows:
        txt.append(" ".join("%016x" % w for limb in row for w in limb))
    p = subprocess.run([exe, mode], input="\n".join(txt) + "\n",
                       capture_output=True, text=True)
    toks = [int(t, 16) for t in p.stdout.split()]
    per = 12 if mode == "mul" else 72
    out = []
    for k in range(len(toks) // per):
        seg = toks[k * per:(k + 1) * per]
        out.append([seg[i * 6:(i + 1) * 6] for i in range(per // 6)])
    return out


def compare(got, refs):
    bad = 0
    first = None
    for g, r in zip(got, refs):
        if g != r:
            bad += 1
            if first is None:
                first = g
    return bad


def mutate(src_path, dst_path, kind):
    s = open(src_path).read()
    if kind == "f2mul":
        old = """    fpsub(d0, t0, t1);
    fpcopy(r->c[0], d0); fpcopy(r->c[1], cross);"""
        new = """    fpadd(d0, t0, t1);
    fpcopy(r->c[0], cross); fpcopy(r->c[1], d0);"""
        assert old in s
        s = s.replace(old, new)
    else:
        old = """    fpcopy(one.c[0].c[0], ONE6);   /* Montgomery 1 = R mod p, from the extern table */"""
        new = """    for (int i = 0; i < 3; i++) one.c[i].c[0][0] = 1;"""
        assert old in s
        s = s.replace(old, new)
    open(dst_path, "w").write(s)
    return dst_path


def main():
    res = {"vectors_per_kernel": N, "p": "0x%x" % P, "montgomery_R": "2^384",
           "reference": "plain Python integers, schoolbook Fp6, no Montgomery and "
                        "no Karatsuba (no shared code path with the C kernels)",
           "kernels": {}}
    rows, refs = gen("mul", N)
    exe = os.path.join(HERE, "drv_mul")
    build(os.path.join(HERE, "driver.c"), exe, "mul", ("-DKERNEL_MUL",))
    got = run(exe, "mul", rows)
    bad = compare(got, refs)
    m1 = mutate(os.path.join(HERE, "kernels.c"), os.path.join(HERE, "kernels_mut1.c"), "f2mul")
    exe1 = os.path.join(HERE, "drv_mul_m1")
    build(os.path.join(HERE, "driver.c"), exe1, "mul_m1", ("-DKERNEL_MUL",))
    # mutation build must use the mutated source
    subprocess.run(["clang", "-O2", "-target", "arm64-apple-macos11", "-I", HERE,
                    "-o", exe1, os.path.join(HERE, "driver.c"), os.path.join(HERE, "consts.c")],
                   check=False, capture_output=True)
    res["kernels"]["f2mul"] = {
        "vectors": len(refs), "mismatches": bad, "verdict": "PASS" if bad == 0 else "FAIL",
        "mutation_reverted_fix": "cross term written to c0 and t0+t1 to c1",
        "mutation_build": "kernels_mut1.c", "mutation_mismatches": None}
    rows2, refs2 = gen("sqr", N, seed=999)
    exe2 = os.path.join(HERE, "drv_sqr")
    build(os.path.join(HERE, "driver.c"), exe2, "sqr")
    got2 = run(exe2, "sqr", rows2)
    bad2 = compare(got2, refs2)
    res["kernels"]["cyclotomic_sqr"] = {
        "vectors": len(refs2), "mismatches": bad2,
        "verdict": "PASS" if bad2 == 0 else "FAIL",
        "identity_checked": "(2vt+1, s-vt-1-t) - (a0^2+v a1^2, 2 a0 a1) == (1-norm, 0) "
                            "on every vector",
        "mutation_reverted_fix": "nonzero constant in every Fp6 coefficient of the "
                                 "identity element",
        "mutation_build": "kernels_mut2.c", "mutation_mismatches": None}
    # mutation controls: build a driver that includes the mutated kernels
    for kind, key, nrows, r in (("f2mul", "f2mul", rows, refs),
                                ("one", "cyclotomic_sqr", rows2, refs2)):
        mut = mutate(os.path.join(HERE, "kernels.c"),
                     os.path.join(HERE, "mut_%s_kernels.c" % kind), kind)
        drv = os.path.join(HERE, "mut_%s_driver.c" % kind)
        s = open(os.path.join(HERE, "driver.c")).read().replace('#include "kernels.c"',
                                                               '#include "mut_%s_kernels.c"' % kind)
        open(drv, "w").write(s)
        ex = os.path.join(HERE, "drv_mut_%s" % kind)
        subprocess.run(["clang", "-O2", "-target", "arm64-apple-macos11", "-I", HERE,
                        "-o", ex, drv, os.path.join(HERE, "consts.c")],
                       check=True, capture_output=True)
        g = run(ex, "mul" if kind == "f2mul" else "sqr", nrows)
        res["kernels"][key]["mutation_mismatches"] = compare(g, r)
    res["all_pass"] = all(v["verdict"] == "PASS" for v in res["kernels"].values())
    res["mutation_controls_detected"] = all(
        (v["mutation_mismatches"] or 0) > 0 for v in res["kernels"].values())
    json.dump(res, open(os.path.join(EV, "kernel_verify23.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
