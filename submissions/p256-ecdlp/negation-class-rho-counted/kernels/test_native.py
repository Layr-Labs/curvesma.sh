"""Native correctness test of p256_walk.c against Python big integers.

  cc -O2 -DWITH_SLOW_PATH_STUB -shared -fPIC -o libp256w.so p256_walk.c   (.dylib on macOS)
  python3 test_native.py ./libp256w.so [seed]

Checks (all against NIST P-256 affine arithmetic in Python):
  fe_mul / fe_add / fe_sub on edge cases and random inputs, fe_inv,
  one full batch round as scheduled in proof.md section 5.1: pass1 over the addition
  list, pass1_dbl over the doubling queue, one shared inversion, pass2_dbl, pass2;
  the history fields, the arrival events, and that an event leaves the history alone.
"""
import ctypes, random, sys
lib = ctypes.CDLL(sys.argv[1] if len(sys.argv) > 1 else "./libp256w.so")
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 7
p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
a = p - 3
b = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
G = (0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,
     0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)
r = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
M64 = 2**64 - 1
R = 2**256
Rinv = pow(R, -1, p)
mont = lambda x: x * R % p
unm = lambda x: x * Rinv % p


def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        l = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        l = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)


def mul(k, P):
    Q = None
    while k:
        if k & 1: Q = add(Q, P)
        P = add(P, P); k >>= 1
    return Q


class FE(ctypes.Structure): _fields_ = [("l", ctypes.c_uint64 * 4)]
class PT(ctypes.Structure): _fields_ = [("x", FE), ("y", FE)]
class W(ctypes.Structure):
    _fields_ = [("P", PT), ("km1", ctypes.c_uint64 * 3), ("km2", ctypes.c_uint64 * 3),
                ("kcp", ctypes.c_uint64 * 3), ("len", ctypes.c_uint64),
                ("flags", ctypes.c_uint64), ("idx", ctypes.c_uint64)]


def tof(x):
    f = FE()
    for i in range(4): f.l[i] = (x >> (64 * i)) & M64
    return f


def fr(f): return sum(f.l[i] << (64 * i) for i in range(4))


pF = ctypes.POINTER(FE); pW = ctypes.POINTER(W); pP = ctypes.POINTER(PT); u64 = ctypes.c_uint64
for fn in (lib.fe_mul, lib.fe_add, lib.fe_sub): fn.argtypes = [pF] * 3; fn.restype = None
lib.fe_inv.argtypes = [pF, pF]
lib.pass1.argtypes = [pW, pP, pF, pF, pF, u64]; lib.pass1.restype = u64
lib.pass2.argtypes = [pW, pP, pF, pF, pF, pF, u64, u64]; lib.pass2.restype = u64
ppW = ctypes.POINTER(pW)
lib.pass1_dbl.argtypes = [ppW, pF, pF, pF, u64]; lib.pass1_dbl.restype = u64
lib.pass2_dbl.argtypes = [ppW, pF, pF, pF, pF, pF, u64, u64]; lib.pass2_dbl.restype = u64
random.seed(seed)

# ---- field kernels
edge = [0, 1, 2, p - 1, p - 2, p - 3, (p + 1) // 2, 2**255, 2**224, 2**192, 2**96, 2**64, 2**64 - 1,
        2**32, R % p, Rinv, p - 2**224, p - 2**192]
vals = edge + [random.randrange(p) for _ in range(400)]
cases = 0
def chk(x, y):
    global cases
    o = FE()
    lib.fe_mul(ctypes.byref(o), ctypes.byref(tof(x)), ctypes.byref(tof(y))); assert fr(o) == x * y * Rinv % p
    lib.fe_add(ctypes.byref(o), ctypes.byref(tof(x)), ctypes.byref(tof(y))); assert fr(o) == (x + y) % p
    lib.fe_sub(ctypes.byref(o), ctypes.byref(tof(x)), ctypes.byref(tof(y))); assert fr(o) == (x - y) % p
    cases += 1
for x in vals[:40]:
    for y in vals[:40]: chk(x, y)
for _ in range(50000): chk(random.randrange(p), random.randrange(p))
for _ in range(50):
    v = random.randrange(1, p); o = FE()
    lib.fe_inv(ctypes.byref(o), ctypes.byref(tof(mont(v)))); assert fr(o) == mont(pow(v, p - 2, p))
print("field kernels: cases", cases, "fe_inv 50 OK")

# ---- one batch round: additions + doubling queue sharing one inversion
RB = 11; RR = 1 << RB
T = (PT * (2 * RR))()
pt_ = mul(random.randrange(1, r), G)
for h in range(RR):
    pt_ = add(pt_, mul(random.randrange(1, 2**40), G)) if h % 64 == 0 else add(pt_, G)
    for t in (0, 1):
        q = pt_ if t == 0 else (pt_[0], (-pt_[1]) % p)
        T[2 * h + t].x = tof(mont(q[0])); T[2 * h + t].y = tof(mont(q[1]))
n, nd, rounds = 100, 7, 4
Wk = (W * n)(); Wd = (W * nd)(); pos = []; posd = []
def newpt():
    q = mul(random.randrange(1, r), G)
    return q if random.random() < 0.5 else (q[0], (-q[1]) % p)
for i in range(n):
    q = newpt(); pos.append(q); Wk[i].P.x = tof(mont(q[0])); Wk[i].P.y = tof(mont(q[1]))
for i in range(nd):
    q = newpt(); posd.append(q); Wd[i].P.x = tof(mont(q[0])); Wd[i].P.y = tof(mont(q[1]))
    Wd[i].flags = 1; Wd[i].km1[0] = 5; Wd[i].kcp[1] = 9
Q = (pW * nd)(*[ctypes.cast(ctypes.byref(Wd[i]), pW) for i in range(nd)])
DX = (FE * n)(); PRE = (FE * n)(); DXd = (FE * nd)(); PREd = (FE * nd)()
one = tof(mont(1)); hist = [[] for _ in range(n)]
fp = lambda m: (m & M64, (m >> 128) & M64, (m >> 192) & M64)
for rd in range(rounds):
    acc = tof(mont(1))
    assert lib.pass1(Wk, T, DX, PRE, ctypes.byref(acc), n) == 0
    assert lib.pass1_dbl(Q, DXd, PREd, ctypes.byref(acc), nd) == 0
    exp = []
    for i in range(n):
        xm = pos[i][0] * R % p; ym = pos[i][1] * R % p
        h = xm & (RR - 1); t = 1 if (p - ym) < ym else 0
        qq = (unm(fr(T[2 * h + t].x)), unm(fr(T[2 * h + t].y)))
        exp.append(add(pos[i], qq)); assert Wk[i].idx == 2 * h + t
    expd = [add(q, q) for q in posd]
    inv = FE(); lib.fe_inv(ctypes.byref(inv), ctypes.byref(acc))
    assert lib.pass2_dbl(Q, DXd, PREd, ctypes.byref(PRE[n - 1]), ctypes.byref(one), ctypes.byref(inv), nd, 1 << 60) == 0
    assert lib.pass2(Wk, T, DX, PRE, ctypes.byref(one), ctypes.byref(inv), n, 1 << 60) == 0
    for i in range(n):
        got = (unm(fr(Wk[i].P.x)), unm(fr(Wk[i].P.y))); assert got == exp[i]
        pos[i] = got; hist[i].append(got[0])
    for i in range(nd):
        got = (unm(fr(Wd[i].P.x)), unm(fr(Wd[i].P.y))); assert got == expd[i]
        posd[i] = got
        assert Wd[i].flags == 0 and Wd[i].kcp[1] == 0 and Wd[i].len == rd + 1
        assert tuple(Wd[i].km1) == fp(got[0] * R % p)   # new class arrived after the clear
        assert tuple(Wd[i].km2) == (0, 0, 0)
        Wd[i].flags = 1    # next round the same walkers double again
for i in range(n):
    assert tuple(Wk[i].km1) == fp(hist[i][-1] * R % p) and tuple(Wk[i].km2) == fp(hist[i][-2] * R % p)
    assert Wk[i].len == rounds
print("batch round: additions", n * rounds, "doublings", nd * rounds, "OK; history, len OK")

# ---- a round with an empty addition list: pass2_dbl starts from Montgomery one
Wd2 = (W * 3)(); pd = []
for i in range(3):
    q = newpt(); pd.append(q); Wd2[i].P.x = tof(mont(q[0])); Wd2[i].P.y = tof(mont(q[1])); Wd2[i].flags = 1
Q2 = (pW * 3)(*[ctypes.cast(ctypes.byref(Wd2[i]), pW) for i in range(3)])
DX2 = (FE * 3)(); PRE2 = (FE * 3)(); acc = tof(mont(1))
assert lib.pass1_dbl(Q2, DX2, PRE2, ctypes.byref(acc), 3) == 0
inv = FE(); lib.fe_inv(ctypes.byref(inv), ctypes.byref(acc))
assert lib.pass2_dbl(Q2, DX2, PRE2, ctypes.byref(one), ctypes.byref(one), ctypes.byref(inv), 3, 1 << 60) == 0
for i in range(3):
    assert (unm(fr(Wd2[i].P.x)), unm(fr(Wd2[i].P.y))) == add(pd[i], pd[i])
print("doubling-only round OK")

# ---- events: a planted K_{-2} match and a length cap hand the walker to the slow path
def one_step(lmax, plant):
    w = (W * 1)(); w[0].P = Wk[0].P
    pre = FE(); dxx = FE(); acc = tof(mont(1))
    lib.pass1(w, T, ctypes.byref(dxx), ctypes.byref(pre), ctypes.byref(acc), 1)
    inv = FE(); lib.fe_inv(ctypes.byref(inv), ctypes.byref(acc))
    h = fr(w[0].P.x) & (RR - 1); ym = fr(w[0].P.y); t = 1 if (p - ym) < ym else 0
    nxt = add((unm(fr(w[0].P.x)), unm(ym)), (unm(fr(T[2 * h + t].x)), unm(fr(T[2 * h + t].y))))
    if plant:
        m = nxt[0] * R % p
        w[0].km2[0] = m & M64; w[0].km2[1] = (m >> 128) & M64; w[0].km2[2] = (m >> 192) & M64
    w[0].km1[0] = 77
    ev = lib.pass2(w, T, ctypes.byref(dxx), ctypes.byref(pre), ctypes.byref(one), ctypes.byref(inv), 1, lmax)
    return w[0], ev
w0, ev = one_step(1 << 60, True)
assert ev == 1 and (w0.flags & 2) and w0.km1[0] == 77, "2-cycle match: event, history untouched"
w0, ev = one_step(1, False)
assert ev == 1 and (w0.flags & 2) and w0.km1[0] == 77, "length cap: event, history untouched"
w0, ev = one_step(1 << 60, False)
assert ev == 0 and w0.km1[0] != 77
print("events: 2-cycle match OK, length cap OK, no event OK; history untouched on event")
print("ALL OK")
