#!/usr/bin/env python3
"""verify25.py -- exact bigint verification of EVERY priced kernel, with per-kernel
mutation controls.

Reference arithmetic is PLAIN Python integers (no Montgomery, no Karatsuba,
schoolbook tower) so it shares no code path with the C kernels.  Every input is a
Montgomery residue a*2^384 mod p; the C output is compared limb-for-limb with the
Montgomery form of the plain-bigint product/inverse.

Kernels covered (>=10000 random vectors each):
  f2mul f2sqr f2inv f6mul f6sqr f6inv f12mul f12sqr f12inv cyclotomic_sqr
cyclotomic_sqr is checked twice: against the Granger-Scott closed form AND against
the true Fp12 square through the exact identity
      r0 - (a0^2 + v a1^2) = 1 - norm,   r1 - 2 a0 a1 = norm - 1,
so on norm-1 inputs the shipped form IS the square.

Mutation controls: for every kernel one targeted corruption of the source is built
and the same vectors must FAIL.  Evidence: ../evidence/kernel_verify25.json
"""
import json
import os
import random
import subprocess

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


def f2inv(a):
    n = (a[0] * a[0] + a[1] * a[1]) % P
    ni = pow(n, -1, P)
    return (a[0] * ni % P, (-a[1]) * ni % P)


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


def f6inv(a):
    A = f2s(f2m(a[0], a[0]), f2xi(f2m(a[1], a[2])))
    B = f2s(f2xi(f2m(a[2], a[2])), f2m(a[0], a[1]))
    C = f2s(f2m(a[1], a[1]), f2m(a[0], a[2]))
    n = f2a(f2m(a[0], A), f2xi(f2a(f2m(a[2], B), f2m(a[1], C))))
    ni = f2inv(n)
    return (f2m(A, ni), f2m(B, ni), f2m(C, ni))


def f12m(a, b):
    c0 = f6a(f6m(a[0], b[0]), f6v(f6m(a[1], b[1])))
    c1 = f6a(f6m(a[0], b[1]), f6m(a[1], b[0]))
    return (c0, c1)


def f12inv(a):
    n = f6s(f6m(a[0], a[0]), f6v(f6m(a[1], a[1])))
    ni = f6inv(n)
    return (f6m(a[0], ni), tuple(f2s((0, 0), x) for x in f6m(a[1], ni)))


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
    return f6s(f6m(a0, a0), f6v(f6m(a1, a1)))


# ------------------------------- vector builders -------------------------------
def flat_f2(a):
    return [limbs(to_mont(x)) for x in a]


def flat_f6(a):
    return [limbs(to_mont(x)) for p in a for x in p]


def flat_f12(a):
    return [limbs(to_mont(x)) for c in a for p in c for x in p]


def rand_f2(rng):
    return (rng.randrange(P), rng.randrange(P))


def rand_f6(rng):
    return tuple(rand_f2(rng) for _ in range(3))


def rand_f12(rng):
    return (rand_f6(rng), rand_f6(rng))


def mont_out_f2(r):
    return [limbs(to_mont(x)) for x in r]


def mont_out_f6(r):
    return [limbs(to_mont(x)) for p in r for x in p]


def mont_out_f12(r):
    return [limbs(to_mont(x)) for c in r for p in c for x in p]


def gen_f2mul(rng):
    a, b = rand_f2(rng), rand_f2(rng)
    return flat_f2(a) + flat_f2(b), mont_out_f2(f2m(a, b))


def gen_f2sqr(rng):
    a = rand_f2(rng)
    return flat_f2(a), mont_out_f2(f2m(a, a))


def gen_f2inv(rng):
    a = rand_f2(rng)
    return flat_f2(a), mont_out_f2(f2inv(a))


def gen_f6mul(rng):
    a, b = rand_f6(rng), rand_f6(rng)
    return flat_f6(a) + flat_f6(b), mont_out_f6(f6m(a, b))


def gen_f6sqr(rng):
    a = rand_f6(rng)
    return flat_f6(a), mont_out_f6(f6m(a, a))


def gen_f6inv(rng):
    a = rand_f6(rng)
    return flat_f6(a), mont_out_f6(f6inv(a))


def gen_f12mul(rng):
    a, b = rand_f12(rng), rand_f12(rng)
    return flat_f12(a) + flat_f12(b), mont_out_f12(f12m(a, b))


def gen_f12sqr(rng):
    a = rand_f12(rng)
    return flat_f12(a), mont_out_f12(f12m(a, a))


def gen_f12inv(rng):
    a = rand_f12(rng)
    return flat_f12(a), mont_out_f12(f12inv(a))


def gen_csqr(rng):
    a = rand_f12(rng)
    r = spec(*a)
    q = true_sq(*a)
    nz = norm(*a)
    assert f6s(r[0], q[0]) == f6s(ONE6, nz), "Granger-Scott r0 identity"
    assert f6s(r[1], q[1]) == f6s(nz, ONE6), "Granger-Scott r1 identity"
    return flat_f12(a), mont_out_f12(r)


# ------------------------------- mutations -------------------------------------
MUT = {
    "f2mul": ("""    fpsub(d0, t0, t1);
    fpcopy(r->c[0], d0); fpcopy(r->c[1], cross);""",
              """    fpadd(d0, t0, t1);
    fpcopy(r->c[0], cross); fpcopy(r->c[1], d0);"""),
    "f2sqr": ("""    fpcopy(r->c[0], t0); fpadd(r->c[1], t1, t1);""",
              """    fpcopy(r->c[0], t0); fpcopy(r->c[1], t1);"""),
    "f2inv": ("""    fpmul(t, a->c[1], ni); fpneg(r->c[1], t);""",
              """    fpmul(t, a->c[1], ni); fpcopy(r->c[1], t);"""),
    "f6mul": ("""    f2mulxi(&k, &t2); f2add(&r->c[1], &u2, &k);
    f2add(&r->c[2], &u1, &t1);""",
              """    f2mulxi(&k, &t2); f2add(&r->c[1], &u1, &k);
    f2add(&r->c[2], &u2, &t1);"""),
    "f6sqr": ("""AI void f6sqr(fp6 *r, const fp6 *a) { f6mul(r, a, a); }""",
              """AI void f6sqr(fp6 *r, const fp6 *a) { f6mul(r, a, a); fpsub(r->c[2].c[0], r->c[2].c[0], ONE6); }"""),
    "f6inv": ("""    f2mul(&r->c[1], &B, &ni);""",
              """    f2sub(&r->c[1], &A, &ni);"""),
    "f12mul": ("""    f6mul(&t3, &a->c[1], &b->c[0]);
    f6add(&s1, &t2, &t3);""",
               """    f6mul(&t3, &a->c[1], &b->c[0]);
    f6sub(&s1, &t2, &t3);"""),
    "f12sqr": ("""AI void f12sqr(fp12 *r, const fp12 *a) { f12mul(r, a, a); }""",
               """AI void f12sqr(fp12 *r, const fp12 *a) { f12mul(r, a, a); f6add(&r->c[0], &r->c[0], &r->c[0]); }"""),
    "f12inv": ("""    f6mul(&t, &a->c[1], &ni); f6neg(&r->c[1], &t);""",
               """    f6mul(&t, &a->c[1], &ni); r->c[1] = t;"""),
    "cyclotomic_sqr": ("""    fpcopy(one.c[0].c[0], ONE6);""",
                       """    for (int i = 0; i < 3; i++) one.c[i].c[0][0] = 1;"""),
}

MODES = {
    "f2mul": ("f2mul", gen_f2mul, 24, 12),
    "f2sqr": ("f2sqr", gen_f2sqr, 12, 12),
    "f2inv": ("f2inv", gen_f2inv, 12, 12),
    "f6mul": ("f6mul", gen_f6mul, 72, 36),
    "f6sqr": ("f6sqr", gen_f6sqr, 36, 36),
    "f6inv": ("f6inv", gen_f6inv, 36, 36),
    "f12mul": ("f12mul", gen_f12mul, 144, 72),
    "f12sqr": ("f12sqr", gen_f12sqr, 72, 72),
    "f12inv": ("f12inv", gen_f12inv, 72, 72),
    "cyclotomic_sqr": ("csqr", gen_csqr, 72, 72),
}


def build(driver, kernels_name, out):
    s = open(os.path.join(HERE, driver)).read().replace(
        '#include "kernels.c"', '#include "%s"' % kernels_name)
    drv = os.path.join(HERE, "_drv_%s.c" % out)
    open(drv, "w").write(s)
    subprocess.run(["clang", "-O0", "-target", "arm64-apple-macos11", "-I", HERE,
                    "-o", os.path.join(HERE, out), drv,
                    os.path.join(HERE, "consts.c")], check=True, capture_output=True)


def run(exe, mode, rows, per):
    txt = "\n".join(" ".join("%016x" % w for limb in row for w in limb) for row in rows)
    p = subprocess.run([os.path.join(HERE, exe), mode], input=txt + "\n",
                       capture_output=True, text=True)
    toks = [int(t, 16) for t in p.stdout.split()]
    assert len(toks) % per == 0, (len(toks), per)
    return [[toks[i * per:(i + 1) * per][j * 6:(j + 1) * 6]
             for j in range(per // 6)]
            for i in range(len(toks) // per)]


def mismatches(got, refs):
    bad = 0
    for g, r in zip(got, refs):
        if g != r:
            bad += 1
    return bad


def main():
    src = open(os.path.join(HERE, "kernels.c")).read()
    res = {"mode": "bigint-tower-sweep", "vectors_per_kernel": N,
           "p": "0x%x" % P, "montgomery_R": "2^384",
           "reference": "plain Python integers, schoolbook Fp6/Fp12, no Montgomery and "
                        "no Karatsuba (no shared code path with the C kernels)",
           "kernels": {}}
    rng = random.Random(20251024)
    build("driver.c", "kernels.c", "drv_all")
    all_pass, all_mut = True, True
    for name, (mode, gen, nin, nout) in MODES.items():
        rows, refs = [], []
        for _ in range(N):
            row, ref = gen(rng)
            assert sum(len(x) for x in row) == nin and sum(len(x) for x in ref) == nout, \
                (name, len(row), len(ref))
            rows.append(row)
            refs.append(ref)
        got = run("drv_all", mode, rows, nout)
        bad = mismatches(got, refs)
        old, new = MUT[name]
        assert old in src, "mutation anchor missing for " + name
        mut_src = src.replace(old, new)
        open(os.path.join(HERE, "mut_kernels.c"), "w").write(mut_src)
        build("driver.c", "mut_kernels.c", "drv_mut")
        gotm = run("drv_mut", mode, rows, nout)
        badm = mismatches(gotm, refs)
        res["kernels"][name] = {
            "vectors": len(refs), "mismatches": bad,
            "verdict": "PASS" if bad == 0 else "FAIL",
            "mutation": MUT[name][1].strip().replace("\n", " ")[:160],
            "mutation_mismatches": badm,
            "mutation_control": "PASS" if badm > 0 else "FAIL"}
        all_pass &= bad == 0
        all_mut &= badm > 0
        print(name, "mismatches", bad, "mutation", badm, flush=True)
    # remove temporary mutated source and binaries from the shipped tree
    for f in ("mut_kernels.c", "drv_all", "drv_mut"):
        p = os.path.join(HERE, f)
        if os.path.exists(p):
            os.remove(p)
    for f in os.listdir(HERE):
        if f.startswith("_drv_") and f.endswith(".c"):
            os.remove(os.path.join(HERE, f))
    res["all_pass"] = all_pass
    res["mutation_controls_detected"] = all_mut
    json.dump(res, open(os.path.join(EV, "kernel_verify25.json"), "w"), indent=1)
    print(json.dumps({k: v["verdict"] for k, v in res["kernels"].items()}, indent=1))
    print("all_pass", all_pass, "all_mut", all_mut)


if __name__ == "__main__":
    main()
