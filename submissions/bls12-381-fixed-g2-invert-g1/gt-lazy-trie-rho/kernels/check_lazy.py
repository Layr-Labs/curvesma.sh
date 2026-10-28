"""Run lazy.s and root.s (gen_lazy.py) in emu.py against a big-integer reference on real BLS12-381 GT
elements: for a walk element x and a trie-node product Pi, lazy.s must return exactly the stored
(Montgomery) a, b, B, a zeta, B zeta of y' = x * Pi, the minimum K of the twelve stored keys and its index
kappa, and kappa must select the conjugate with the minimum key among the twelve Frobenius conjugates of y'
computed explicitly; root.s the same for y' = x, the six pre-sums and xi_x. Every call must execute the same count, and every memory
operand of both listings must be an encodable AArch64 immediate.
Usage: python3 -I check_lazy.py <n_cases> <seed>"""
import json, os, random, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import emu
import check as C
import gen_lazy as L

P = C.P
RM = C.RM
OUT, XB, EB, PB, SP = 0x10000, 0x20000, 0x30000, 0x40000, 0x80000


def wflat(a):
    """blst flat vector -> w^m coefficients (Fp2 pairs), m = 0..5."""
    v = C.flat(a)
    return [(v[L.fidx(m, 0)], v[L.fidx(m, 1)]) for m in range(6)]


def frob_keys_direct(a):
    """Re(w^1 coefficient) of a^(p^k), k = 0..11, by explicit Frobenius maps (the definition)."""
    out = []
    y = a
    for k in range(12):
        out.append(wflat(y)[1][0])
        y = C.f12pow(y, P)
    return out


def gamma(m):
    """Node vector of Pi = m: gamma_0 = m_1, gamma_1 = m_0, gamma_k = xi m_{7-k} (k = 2..5)."""
    w = wflat(m)
    xi = lambda c: ((c[0] - c[1]) % P, (c[0] + c[1]) % P)
    g = [w[1], w[0]] + [xi(w[7 - k]) for k in range(2, 6)]
    vec = []
    for re, im in g:
        re, im = re * RM % P, im * RM % P
        vec += [re, im, re + im]
    xr, xi = xi_sum(g)
    c = [4 * P - xr * RM % P, 2 * P - xi * RM % P]
    assert 3 * P < c[0] <= 4 * P and P < c[1] <= 2 * P
    return vec + c


ENC = re.compile(r"\s*(ldp|stp|ldr|str)\s+.*\[(\w+)(?:,\s*#(-?\d+))?\]")


def encodable(text):
    """Every memory operand is an encodable AArch64 immediate: ldp/stp 64-bit in [-512, 504], a multiple of 8;
    ldr/str scaled in [0, 32760], a multiple of 8, or unscaled in [-256, 255]. Returns the operand count."""
    n = 0
    for line in text.splitlines():
        m = ENC.match(line)
        if m:
            op, off = m.group(1), int(m.group(3) or 0)
            n += 1
            assert ((-512 <= off <= 504 and off % 8 == 0) if op in ("ldp", "stp")
                    else (0 <= off <= 32760 and off % 8 == 0) or -256 <= off <= 255), line
    return n


def xi_sum(c):
    """sum_i c_{2i} c_{2i+1} in Fp2 (u^2 = -1), plain residues."""
    re = sum(c[2 * i][0] * c[2 * i + 1][0] - c[2 * i][1] * c[2 * i + 1][1] for i in range(3)) % P
    im = sum(c[2 * i][0] * c[2 * i + 1][1] + c[2 * i][1] * c[2 * i + 1][0] for i in range(3)) % P
    return re, im


def run(ops, pool, xv, entry, rng):
    m = emu.Machine()
    for i in range(31):
        m.r[i] = rng.getrandbits(64)
    saved = list(m.r[19:31])
    m.mem.update(pool)
    xm = [(v * RM) % P for v in C.flat(xv)]
    C.put(m.mem, XB, xm)
    w = wflat(xv)
    C.put(m.mem, XB + L.XS, [(re * RM % P) + (im * RM % P) for re, im in w])
    C.put(m.mem, XB + L.XI, [v * RM % P for v in xi_sum(w)])
    if entry is not None:
        C.put(m.mem, EB, entry)
    m.r[0], m.r[1], m.r[2], m.r[3], m.r[31] = OUT, XB, EB, PB, SP
    n = emu.run(ops, m)
    assert m.r[31] == SP and m.r[19:31] == saved, "callee-saved state or sp not restored"
    rd = lambda off: sum(m.mem[OUT + off + 8 * j] << (64 * j) for j in range(6))
    res = {k: rd(L.OUT[k]) for k in ("a", "b", "B", "az", "Bz", "K")}
    res["kap"] = m.mem[OUT + L.OUT["kap"]]
    xs = [sum(m.mem[XB + L.XS + 48 * k + 8 * j] << (64 * j) for j in range(6)) for k in range(8)]
    return n, res, xs


def reference(y):
    w = wflat(y)
    a, b = w[1]
    d1 = L.D1
    B = (a * d1[0] + b * d1[1]) % P
    st = lambda v: v * RM % P
    keys_int = frob_keys_direct(y)
    keys = [st(v) for v in keys_int]
    K = min(keys)
    return {"a": st(a), "b": st(b), "B": st(B), "az": st(a * L.ZETA % P), "Bz": st(B * L.ZETA % P),
            "K": K, "kap_set": [k for k in range(12) if keys[k] == K]}


def main():
    n_cases, seed = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(seed)
    pool = L.pool_words()
    pool = {PB + k: v for k, v in pool.items()}
    progs, enc = {}, {}
    for name in ("lazy", "root"):
        text = open(os.path.join(HERE, f"{name}.s")).read()
        enc[name] = encodable(text)
        progs[name] = emu.compile_prog(emu.parse(text))
    e = (P ** 12 - 1) // C.Rr
    gens = [C.f12pow(C.rand12(rng), e) for _ in range(3)]
    counts = {"lazy": set(), "root": set()}
    fails = {"lazy": 0, "root": 0}
    x, m = gens[0], gens[1]
    kinds = {}
    for i in range(n_cases):
        kind = ("gt", "gt", "rand", "edge")[i % 4]
        kinds[kind] = kinds.get(kind, 0) + 1
        if kind == "gt":
            x = C.f12mul(x, gens[2])
            m = C.f12mul(m, gens[1])
        elif kind == "rand":
            x, m = C.rand12(rng), C.rand12(rng)
        else:
            x = C.unflat([rng.choice((0, 1, P - 1, P - 2, rng.randrange(P))) for _ in range(12)])
            m = C.unflat([rng.choice((0, 1, P - 1, P - 2, rng.randrange(P))) for _ in range(12)])
        for name in ("lazy", "root"):
            if name == "lazy":
                n, got, _ = run(progs[name], pool, x, gamma(m), rng)
                y = C.f12mul(x, m)
            else:
                n, got, xs = run(progs[name], pool, x, None, rng)
                y = x
                w = wflat(x)
                if xs != [(re * RM % P) + (im * RM % P) for re, im in w] + [v * RM % P for v in xi_sum(w)]:
                    fails[name] += 1
            ref = reference(y)
            counts[name].add(n)
            ok = all(got[k] == ref[k] for k in ("a", "b", "B", "az", "Bz", "K")) and got["kap"] in ref["kap_set"]
            if not ok:
                fails[name] += 1
                if fails[name] <= 3:
                    print("MISMATCH", name, {k: (got[k], ref.get(k)) for k in got if got[k] != ref.get(k)})
    res = {"cases": kinds, "encodable_memory_operands": enc, "executed_counts": {k: sorted(v) for k, v in counts.items()}, "fails": fails}
    print(json.dumps(res))
    assert all(len(v) == 1 for v in counts.values()) and not any(fails.values())
    print("PASS")


if __name__ == "__main__":
    main()
