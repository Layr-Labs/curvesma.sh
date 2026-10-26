"""Run mul.s and sqr.s in emu.py against a Python big-integer reference of blst's tower in Montgomery
form. Every output limb must match; the executed count must be the same on every call.
Usage: python3 -I check.py <mul|sqr> <n_gt> <n_extra> <seed>"""
import json, os, random, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import emu

P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
Rr = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
RM = (1 << 384) % P


# reference tower, schoolbook formulas (independent of the kernels' Karatsuba)
def f2mul(a, b):
    return ((a[0] * b[0] - a[1] * b[1]) % P, (a[0] * b[1] + a[1] * b[0]) % P)


def f2add(a, b):
    return ((a[0] + b[0]) % P, (a[1] + b[1]) % P)


def f2xi(a):
    return ((a[0] - a[1]) % P, (a[0] + a[1]) % P)


def f6mul(a, b):
    c = [(0, 0)] * 5
    for i in range(3):
        for j in range(3):
            c[i + j] = f2add(c[i + j], f2mul(a[i], b[j]))
    return (f2add(c[0], f2xi(c[3])), f2add(c[1], f2xi(c[4])), c[2])


def f6add(a, b):
    return tuple(f2add(a[i], b[i]) for i in range(3))


def f6v(a):
    return (f2xi(a[2]), a[0], a[1])


def f12mul(a, b):
    t0, t1 = f6mul(a[0], b[0]), f6mul(a[1], b[1])
    c1 = f6add(f6mul(a[0], b[1]), f6mul(a[1], b[0]))
    return (f6add(t0, f6v(t1)), c1)


def f12pow(a, e):
    r = ((( 1, 0), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0)))
    for bit in bin(e)[2:]:
        r = f12mul(r, r)
        if bit == "1":
            r = f12mul(r, a)
    return r


def flat(a):
    return [a[i][j][k] for i in range(2) for j in range(3) for k in range(2)]


def unflat(v):
    return tuple(tuple((v[6 * i + 2 * j], v[6 * i + 2 * j + 1]) for j in range(3)) for i in range(2))


def rand12(rng):
    return unflat([rng.randrange(P) for _ in range(12)])


def limbs(v):
    return [(v >> (64 * i)) & ((1 << 64) - 1) for i in range(6)]


def put(mem, base, vals, n=6):
    for i, v in enumerate(vals):
        for j in range(n):
            mem[base + 8 * (n * i + j)] = (v >> (64 * j)) & ((1 << 64) - 1)


OUT, XB, EB, PB, SP = 0x10000, 0x20000, 0x30000, 0x40000, 0x80000


def setup(kind):
    lay = json.load(open(os.path.join(HERE, f"{kind}.layout.json")))
    ops = emu.compile_prog(emu.parse(open(os.path.join(HERE, f"{kind}.s")).read()))
    pool = {}
    pl = lay["pool"]
    put(pool, PB, [P])
    pool[PB + 48] = int(pl["n0"])
    put(pool, PB + 56, [2 * P, 4 * P])
    for k, K in lay["K"].items():
        put(pool, PB + int(pl["K_base"]) + 96 * int(k), [int(K)], n=12)
    return lay, ops, pool


def run_one(ops, pool, lay, xv, mv, rng):
    m = emu.Machine()
    for i in range(31):
        m.r[i] = rng.getrandbits(64)
    saved = list(m.r[19:31])
    m.mem.update(pool)
    xm = [(v * RM) % P for v in flat(xv)]
    put(m.mem, XB, xm)
    if mv is not None:
        mm = [(v * RM) % P for v in flat(mv)]
        env = {f"m{i}": mm[i] for i in range(12)}
        ent = [sum(c * env[n] for n, c in L.items()) for L in lay["entry"]]
        assert all(0 <= e < (1 << 384) for e in ent)
        put(m.mem, EB, ent)
    m.r[0], m.r[1], m.r[2], m.r[3], m.r[31] = OUT, XB, EB, PB, SP
    n = emu.run(ops, m)
    assert m.r[31] == SP and m.r[19:31] == saved, "callee-saved state or sp not restored"
    got = [sum(m.mem[OUT + 8 * (6 * i + j)] << (64 * j) for j in range(6)) for i in range(12)]
    return n, got


def main():
    kind, n_gt, n_extra, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    rng = random.Random(seed)
    lay, ops, pool = setup(kind)
    t = time.time()
    if kind == "mul":
        e = (P ** 12 - 1) // Rr
        gens = [f12pow(rand12(rng), e) for _ in range(4)]
        one = unflat([1] + [0] * 11)
        assert all(f12pow(g, Rr) == one for g in gens), "GT generator is not an r-th root of unity"
        cases = []
        x, y = gens[0], gens[1]
        for i in range(n_gt):
            x = f12mul(x, gens[2])
            y = f12mul(y, gens[3 if i % 2 else 0])
            cases.append(("gt", x, y))
        for i in range(n_extra):
            if i % 3 == 0:     # coordinates at the extremes of [0, p)
                cases.append(("edge", unflat([rng.choice((0, 1, P - 1, P - 2)) for _ in range(12)]),
                              unflat([rng.choice((0, 1, P - 1, P - 2)) for _ in range(12)])))
            else:
                cases.append(("rand", rand12(rng), rand12(rng)))
    else:
        e = (P ** 6 - 1) * (P ** 2 + 1)
        gens = [f12pow(rand12(rng), e) for _ in range(4)]
        for g in gens:     # cyclotomic: g^(p^4 - p^2 + 1) * ... check norm property g^(p^6+1) = 1
            assert f12pow(g, P ** 6 + 1) == unflat([1] + [0] * 11)
        cases = []
        x = gens[0]
        for i in range(n_gt):
            x = f12mul(x, gens[1 + i % 3])
            cases.append(("cyc", x, None))
    tgen = time.time() - t
    counts, fails, kinds = set(), 0, {}
    t = time.time()
    for kname, xv, mv in cases:
        n, got = run_one(ops, pool, lay, xv, mv, rng)
        ref = f12mul(xv, mv if mv is not None else xv)
        want = [(v * RM) % P for v in flat(ref)]
        counts.add(n)
        kinds[kname] = kinds.get(kname, 0) + 1
        if got != want:
            fails += 1
            if fails <= 3:
                print("MISMATCH", kname, [i for i in range(12) if got[i] != want[i]])
    res = {"kernel": kind, "cases": kinds, "fails": fails, "executed_counts": sorted(counts),
           "static_count": lay["executed"], "seconds_gen": round(tgen, 1), "seconds_run": round(time.time() - t, 1)}
    print(json.dumps(res))
    assert len(counts) == 1 and counts == {lay["executed"]}, "executed count differs between calls"
    assert fails == 0
    print("PASS")


if __name__ == "__main__":
    main()
