"""Native test of secp_walk.c against Python big integers and affine secp256k1 arithmetic.

usage: python3 test_native.py ./libsecpw.so
"""
import ctypes, random, sys

P = 2**256 - 2**32 - 977
N = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
G = (0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
     0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8)
M64 = 2**64 - 1

lib = ctypes.CDLL(sys.argv[1])
L4 = ctypes.c_uint64 * 4


class PT(ctypes.Structure):
    _fields_ = [("x", L4), ("y", L4)]


class FE(ctypes.Structure):
    _fields_ = [("v", L4)]


def tl(v):
    return L4(*[(v >> (64 * i)) & M64 for i in range(4)])


def fl(l):
    return sum(int(l[i]) << (64 * i) for i in range(4))


def add(A, B):
    if A is None:
        return B
    if B is None:
        return A
    if A[0] == B[0]:
        if (A[1] + B[1]) % P == 0:
            return None
        lam = 3 * A[0] * A[0] * pow(2 * A[1], -1, P) % P
    else:
        lam = (B[1] - A[1]) * pow(B[0] - A[0], -1, P) % P
    x = (lam * lam - A[0] - B[0]) % P
    return (x, (lam * (A[0] - x) - A[1]) % P)


def mul(k, A):
    R = None
    while k:
        if k & 1:
            R = add(R, A)
        A = add(A, A)
        k >>= 1
    return R


def ucall(fn, *a):
    out = L4()
    fn(out, *[tl(v) for v in a])
    return fl(out)


random.seed(20261024)
edge = [0, 1, 2, P - 1, P, P + 1, 2**256 - 1, 2**256 - 2, 2**255, 2**192, 2**64 - 1, 2**128 - 1,
        2**33, 2**32 + 977, P - 2**33, 2**256 - 2**33, 2**128, 2**64, (2**64 - 1) << 192]
vals = (edge + [random.getrandbits(256) for _ in range(20000)]
        + [P - random.getrandbits(40) for _ in range(2000)]
        + [2**256 - 1 - random.getrandbits(40) for _ in range(2000)]
        + [random.getrandbits(64) << (64 * random.randrange(4)) for _ in range(2000)])
bad = 0
n = 0
for a in edge:
    for b in edge:
        r = ucall(lib.fe_mul_f, a, b)
        n += 1
        if not (r < 2**256 and r % P == a * b % P):
            bad += 1
for _ in range(60000):
    a, b = random.choice(vals), random.choice(vals)
    r = ucall(lib.fe_mul_f, a, b)
    n += 1
    if not (r < 2**256 and r % P == a * b % P):
        bad += 1
print("fe_mul_f", n, "checked")
n2 = 0
for a in vals:
    r = ucall(lib.fe_sqr_f, a)
    n2 += 1
    if not (r < 2**256 and r % P == a * a % P):
        bad += 1
print("fe_sqr_f", n2, "checked")
n3 = 0
for a in edge[1:] + [random.randrange(1, P) for _ in range(300)]:
    if a % P == 0:
        continue
    r = ucall(lib.fe_inv, a)
    n3 += 1
    if not (r < 2**256 and r * a % P == 1):
        bad += 1
print("fe_inv", n3, "checked")

lib.batch_add.argtypes = [ctypes.c_size_t, ctypes.POINTER(PT), ctypes.POINTER(ctypes.POINTER(PT)), ctypes.POINTER(FE)]
tab = [mul(random.randrange(1, N), G) for _ in range(12)] + [(G[0], G[1]), (G[0], P - G[1])]
tot = 0
bb = 0
for nn in [2, 3, 17, 257, 1000]:
    for rep in range(3):
        W, T, tps, want = (PT * nn)(), (PT * nn)(), (ctypes.POINTER(PT) * nn)(), []
        for i in range(nn):
            a = mul(random.randrange(1, N), G) if random.random() < 0.3 else G
            b = tab[random.randrange(len(tab))]
            if a[0] == b[0]:
                b = tab[0] if tab[0][0] != a[0] else tab[1]
            if random.random() < 0.5:
                b = (b[0], (P - b[1]) % P)
            W[i].x, W[i].y, T[i].x, T[i].y = tl(a[0]), tl(a[1]), tl(b[0]), tl(b[1])
            tps[i] = ctypes.pointer(T[i])
            want.append(add(a, b))
        acc = (FE * nn)()
        lib.batch_add(nn, W, tps, acc)
        for i in range(nn):
            got = (fl(W[i].x), fl(W[i].y))
            tot += 1
            if got != want[i]:
                bb += 1
print("batch_add", tot, "walks checked,", bb, "mismatches (outputs compared as canonical integers)")
print("RESULT", "PASS" if bad == 0 and bb == 0 else "FAIL", "field mismatches:", bad)
