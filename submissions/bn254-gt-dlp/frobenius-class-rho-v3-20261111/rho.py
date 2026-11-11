"""Executable mathematical reference for the Frobenius-class rho. Python is NOT
the scored machine; the priced scalar implementation and its finite instruction
bounds are in machine.md and arithmetic.c.

sample_scalar, fp12_mul, fp12_pow and montgomery_word_reference are carried over
unchanged from the organizer's reviewed baseline rho-baseline-v2-20261110
(credited in proof.md). New here: Frobenius matrices, the canonical class
representative, the stack (Nivasch) cycle detector and the fruitless-cycle
escape. Instrumentation counters (escapes, wasted steps, depth) do not change
the walk.
"""
import hashlib

R_T = 1024          # increment-table size (10-bit index)
D_MAX = 4096        # stack-depth cap (abort = failure)
TIE_MAX = 1024      # cap on full-key tie resolutions (abort = failure)
ORDER = (1, 7, 0, 2, 3, 4, 5, 6, 8, 9, 10, 11)  # coefficient priority of the order


def sample_scalar(r, coins):
    for _ in range(1024):
        v = coins.getrandbits(r.bit_length())
        if v < r:
            return v
    raise RuntimeError("bounded randomness rejection")


def fp12_mul(a, b, p):
    c = [0] * 23
    for i in range(12):
        for j in range(12):
            c[i + j] = (c[i + j] + a[i] * b[j]) % p
    for k in range(22, 11, -1):
        c[k - 6] = (c[k - 6] + 18 * c[k]) % p
        c[k - 12] = (c[k - 12] - 82 * c[k]) % p
    return tuple(c[:12])


def fp12_pow(a, n, p):
    out = (1,) + (0,) * 11
    while n:
        if n & 1:
            out = fp12_mul(out, a, p)
        a = fp12_mul(a, a, p)
        n >>= 1
    return out


def montgomery_word_reference(a, b, modulus):
    """Exact fixed-loop word algorithm matching arithmetic.c, including carries."""
    mask = (1 << 64) - 1
    aa = [(a >> (64 * i)) & mask for i in range(4)]
    bb = [(b >> (64 * i)) & mask for i in range(4)]
    mm = [(modulus >> (64 * i)) & mask for i in range(4)]
    n0 = (-pow(mm[0], -1, 1 << 64)) & mask
    t = [0] * 9
    for i in range(4):
        carry = 0
        for j in range(4):
            v = aa[i] * bb[j] + t[i + j] + carry
            assert v < 1 << 128
            t[i + j], carry = v & mask, v >> 64
        for k in range(i + 4, 9):
            v = t[k] + carry
            t[k], carry = v & mask, v >> 64
        assert carry == 0
    for i in range(4):
        q = (t[i] * n0) & mask
        carry = 0
        for j in range(4):
            v = q * mm[j] + t[i + j] + carry
            assert v < 1 << 128
            t[i + j], carry = v & mask, v >> 64
        for k in range(i + 4, 9):
            v = t[k] + carry
            t[k], carry = v & mask, v >> 64
        assert carry == 0 and t[i] == 0
    u = sum(t[i + 4] << (64 * i) for i in range(5))
    assert 0 <= u < 2 * modulus < 1 << 255
    return u - modulus if u >= modulus else u


# ---------------------------------------------------------------- Frobenius

def frobenius_matrices(p):
    """M[k][i][j]: coefficient of w^i in Frob^k(w^j), Frob(y)=y^p, k=0..11."""
    w = (0, 1) + (0,) * 10
    w1 = fp12_pow(w, p, p)
    cols, cur = [], (1,) + (0,) * 11
    for _ in range(12):
        cols.append(cur)
        cur = fp12_mul(cur, w1, p)
    m1 = [[cols[j][i] for j in range(12)] for i in range(12)]
    mats = [[[int(i == j) for j in range(12)] for i in range(12)]]
    for _ in range(11):
        prev = mats[-1]
        mats.append([[sum(m1[i][t] * prev[t][j] for t in range(12)) % p
                      for j in range(12)] for i in range(12)])
    return mats


def apply(mat, y, p):
    return tuple(sum(mat[i][j] * y[j] for j in range(12) if mat[i][j]) % p
                 for i in range(12))


def order_key(y, p):
    """Total order on Fp12: compare canonical Montgomery residues yR mod p in
    coefficient priority ORDER (the internal representation that is stored)."""
    rm = pow(2, 256, p)
    return tuple(y[i] * rm % p for i in ORDER)


def canon_fp12(z, mats, p):
    """Return (Frob^k(z), k, tie) with Frob^k(z) minimal under order_key.
    Keys use only rows 1 and 7 (which read only coefficients 1 and 7);
    a tie on that two-coefficient key is resolved by full keys."""
    rm = pow(2, 256, p)
    keys = []
    for k in range(12):
        m = mats[k]
        c1 = (m[1][1] * z[1] + m[1][7] * z[7]) % p
        c7 = (m[7][1] * z[1] + m[7][7] * z[7]) % p
        keys.append((c1 * rm % p, c7 * rm % p))
    kmin = min(keys)
    ks = [k for k in range(12) if keys[k] == kmin]
    if len(ks) == 1:
        return apply(mats[ks[0]], z, p), ks[0], False
    best = min(ks, key=lambda k: (order_key(apply(mats[k], z, p), p), k))
    return apply(mats[best], z, p), best, True


# ---------------------------------------------------------------- the walk

def class_rho(r, g, z, ops, coins, cap, table_size=R_T):
    """One bounded Frobenius-class walk with stack cycle detection.

    ops: mul(y,m), sqr(y), pw(base,e), canon(y)->(rep,k,tie), key(y),
    enc(y)->bytes, one, lam (list: rep = y^lam[k]).
    Returns a dict with outcome in {identity, success, cap, stack_cap,
    tie_budget, verification_failed}; every produced element counts as a step.
    """
    st = {"outcome": None, "x": None, "steps": 0, "escapes": 0, "two_cycles": 0,
          "wasted_steps": 0, "ties": 0, "max_depth": 0, "final": None}
    if z == ops["one"]:
        st.update(outcome="identity", x=0)
        return st
    lam, mul, pw, canon, key = ops["lam"], ops["mul"], ops["pw"], ops["canon"], ops["key"]
    salt = coins.getrandbits(256).to_bytes(32, "big")
    table = []
    for _ in range(table_size):
        c, d = sample_scalar(r, coins), sample_scalar(r, coins)
        table.append((mul(pw(g, c), pw(z, d)), c, d))
    a, b = sample_scalar(r, coins), sample_scalar(r, coins)
    y, k, tie = canon(mul(pw(g, a), pw(z, b)))
    a, b = a * lam[k] % r, b * lam[k] % r
    ties = int(tie)
    stack = [(key(y), y, a, b, 0)]
    escape = False
    steps = 0
    while True:
        if steps >= cap:
            st["outcome"] = "cap"
            break
        steps += 1
        if escape:                       # fruitless-cycle escape: square the cycle minimum
            nz, na, nb = ops["sqr"](y), 2 * a, 2 * b
        else:
            j = int.from_bytes(hashlib.sha256(salt + ops["enc"](y)).digest()[:2], "big") % table_size
            m, c, d = table[j]
            nz, na, nb = mul(y, m), a + c, b + d
        escape = False
        y, k, tie = canon(nz)
        a, b = na * lam[k] % r, nb * lam[k] % r
        ties += tie
        if ties > TIE_MAX:
            st["outcome"] = "tie_budget"
            break
        kk = key(y)
        while stack and stack[-1][0] > kk:
            stack.pop()
        if stack and stack[-1][0] == kk:
            _, _, A, B, t0 = stack[-1]
            if (B - b) % r:
                # g^a z^b = g^A z^B  =>  x = (a - A)/(B - b) mod r
                x = (a - A) * pow((B - b) % r, r - 2, r) % r
                ok = pw(g, x) == z
                st.update(outcome="success" if ok else "verification_failed", x=x if ok else None)
                break
            assert (A - a) % r == 0      # B=b forces a=A since g has order r
            escape = True
            st["escapes"] += 1
            st["wasted_steps"] += steps - t0
            st["two_cycles"] += (steps - t0 == 2)
            continue
        stack.append((kk, y, a, b, steps))
        st["max_depth"] = max(st["max_depth"], len(stack))
        if len(stack) > D_MAX:
            st["outcome"] = "stack_cap"
            break
    st.update(steps=steps, ties=ties, final=(y, a, b))
    return st


def toy_ops(q, lam):
    """Additive prime-order group Z/q, generator 1, automorphisms y -> lam[k]*y."""
    def canon(y):
        rep, k = min((lam[k] * y % q, k) for k in range(len(lam)))
        return rep, k, False
    return {"mul": lambda u, v: (u + v) % q, "sqr": lambda u: 2 * u % q,
            "pw": lambda base, e: base * e % q, "canon": canon, "key": lambda u: u,
            "enc": lambda u: u.to_bytes(384, "big"), "one": 0, "lam": lam}


def bn254_ops(p, r, mats):
    rm = pow(2, 256, p)
    return {"mul": lambda u, v: fp12_mul(u, v, p), "sqr": lambda u: fp12_mul(u, u, p),
            "pw": lambda base, e: fp12_pow(base, e, p),
            "canon": lambda u: canon_fp12(u, mats, p), "key": lambda u: order_key(u, p),
            "enc": lambda u: b"".join((c * rm % p).to_bytes(32, "big") for c in u),
            "one": (1,) + (0,) * 11, "lam": [pow(p, k, r) for k in range(12)]}
