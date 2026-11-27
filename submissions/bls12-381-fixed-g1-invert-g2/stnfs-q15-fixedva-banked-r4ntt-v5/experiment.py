#!/usr/bin/env python3
"""Standalone Docker-Safe Verifier for stnfs-q15-fixedva-banked-r4ntt-v5 (stdlib only, 0 external files)."""
from __future__ import annotations
import bisect, json, math, random, sys, time
from fractions import Fraction

X_BLS = -0xD201000000010000
ABS_X, LOG2_ABS_X = abs(X_BLS), math.log2(abs(X_BLS))
P_BLS = int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab", 16)
R_BLS = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
P0_X = int("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb", 16)
P0_Y = int("08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1", 16)
Q0_X = (int("024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8", 16), int("13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e", 16))
Q0_Y = (int("0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801", 16), int("0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be", 16))
C_COST, V_A_SH, SLAB_BASE = 20000.0, 0.9351 * (2.0**56), 1 << 56

REF: dict = {
  "pkg": "stnfs-q15-fixedva-banked-r4ntt-v5", "tag": "CAND5_DUSART_PR16_7080",
  "Bf": 70.80, "Bg": 71.16, "eta": 0.88, "p_retry": 0.0625, "af": 5.3207, "ag": -1.0339, "cush": -0.04,
  "qa": 0.65, "rbr": 0.82, "lf": 408.0, "lg": 688.0, "oe": 10.0, "wb3": 16,
  "b": {"pglo": 28.0, "pw1": 54.0, "pw4": 57.0, "pghi0": 61.8, "pghi1": 59.8, "pflo": 32.0, "pfhi0": 61.4, "pfhi1": 59.4, "qg": 3.4, "qf": 3.8},
  "mc": {
    "sl": 1.7614281390842421, "Vk4": 119.18753166306236, "Rret": 62.63014818225149, "Rreq": 62.63117622977014,
    "sh": 0.0018108318276692436, "se": 0.001126283202156246, "R0": 1350.07745547612, "Ra": 1719.4681182270274,
    "fv": 0.08349025, "mnf": 541.8378473462442, "mng": 670.0026212296631, "dnf": 651.0005601859636, "dng": 685.4442160333803,
    "L4": 108.66281711928426, "g4": 0.0016381482873693107, "pex": 0.0036373788824551894, "cpad": 0.005257195488331446, "Vkch": 119.20048833973382,
  },
  "pol": {"lam": 332708.3045915901, "ret": 0.6858739394706085, "sg": -13.966475955810715, "sf": -8.677562400876722, "s1": -22.644038356687435},
  "s1b": {"p1": (127 << 54) + 1, "g1": 3, "p2": (114 << 54) + 1, "g2": 7, "p3": (108 << 54) + 1, "g3": 5, "bfly": 11.125, "pw": 12.0, "sc": 12.0, "c1": 58.0, "c2": 90.0, "saf": 1.10, "dm": 1.67},
  "rep": {"t12": 18.25, "t3": 20.25, "bkt": 5.5, "hdr": 12.0},
  "la": {"m": 2048, "n": 32, "p": 4, "nl": 0.00150, "sw": 0.190, "x4": 12.90, "x1": 16.0, "rv": 940.0, "ov": 1.06},
  "sge": [(16, 0.005753, 0.001233, 0.31206), (18, 0.002721, 0.000415, 0.12060), (20, 0.002397, 0.000403, 0.16816), (22, 0.001577, 0.000490, 0.12874), (23, 0.001005, 0.000097, 0.08283), (24, 0.001045, 0.000104, 0.12454)],
  "cr": {
    "relation_collection": 107.8180, "linear_algebra": 104.2520, "replay_a_tables": 103.6790,
    "snell_resume_tail": 101.9690, "both_individual_logs_exp": 96.6500, "ecm_split_expected": 96.0500,
    "meter_safety_buffer": 93.4001, "band_split_and_trees": 90.0001, "canonical_a_domain": 84.0001,
    "schirokauer_maps": 80.0001, "final_log_backsub_and_gcd": 72.0001, "matrix_generator_sge": 62.5120,
    "stage0_setup_and_iso": 45.0001, "final_verification_Q": 40.0001,
  },
  "Kg": 7194374076228249696, "num": 106677, "den": 1000, "ram_ub": 117.5,
}

class ExactDickman:
  def __init__(self, umax: float = 30.0, N: int = 200) -> None:
    self.umax, self.N = umax, N; h, tot = 1.0 / N, int(round(umax * N)) + 1
    v = [1.0] * (N + 1) + [0.0] * (tot - N - 1)
    for i in range(N + 1, min(2 * N + 1, tot)): v[i] = 1.0 - math.log(i * h)
    for i in range(2 * N + 1, tot): v[i] = (h * (0.5 * v[i - N] + sum(v[i - N + 1 : i]))) / (i * h - 0.5 * h)
    self.v = v
  def rho(self, u: float) -> float:
    if u <= 1.0: return 1.0
    if u >= self.umax: return 1e-300
    pos = u * self.N; idx = int(pos)
    if idx + 1 >= len(self.v): return self.v[-1]
    f = pos - idx; return math.exp((1.0 - f) * math.log(self.v[idx]) + f * math.log(self.v[idx + 1]))

DK = ExactDickman()

def _primes(n: int) -> list[int]:
  ip = bytearray(b"\x01") * (n + 1); ip[0] = ip[1] = 0
  for i in range(2, int(n**0.5) + 1):
    if ip[i]: ip[i * i : n + 1 : i] = bytearray(len(range(i * i, n + 1, i)))
  return [i for i in range(2, n + 1) if ip[i]]

def pi_dus(B: float) -> float:
  lnB = B * math.log(2.0); return ((2.0**B) / lnB) * (1.0 + 1.0 / lnB + 2.334 / (lnB * lnB))

def U_ub(Bf: float, Bg: float) -> float:
  return (pi_dus(Bf) + pi_dus(Bg)) / 12.0 + 36.0 * pi_dus(0.5 * Bf) + 12.0 * pi_dus(0.5 * Bg)

def ppb_dus(lo: float, hi: float) -> float:
  th_hi = (2.0**hi) * (1.0 + 0.01 / max(1.0, hi * math.log(2.0)))
  th_lo = 0.0 if lo <= 1.0 else (2.0**lo) * (1.0 - 0.01 / max(1.0, lo * math.log(2.0)))
  return max(1.0, (th_hi - th_lo) / math.log(2.0))

# 1. Stage 0 Algebra
def _g1a(p1, p2):
  if p1 is None: return p2
  if p2 is None: return p1
  x1, y1, x2, y2 = p1[0], p1[1], p2[0], p2[1]
  if x1 == x2:
    if (y1 + y2) % P_BLS == 0: return None
    lam = (3 * x1 * x1 * pow(2 * y1, -1, P_BLS)) % P_BLS
  else: lam = ((y2 - y1) * pow(x2 - x1, -1, P_BLS)) % P_BLS
  x3 = (lam * lam - x1 - x2) % P_BLS; return (x3, (lam * (x1 - x3) - y1) % P_BLS)

def check_stage0() -> dict:
  assert (X_BLS - 1) ** 2 * (X_BLS**4 - X_BLS**2 + 1) // 3 + X_BLS == P_BLS and X_BLS**4 - X_BLS**2 + 1 == R_BLS
  lam = X_BLS % R_BLS
  assert P_BLS % R_BLS == lam and (pow(lam, 4, R_BLS) - pow(lam, 2, R_BLS) + 1) % R_BLS == 0
  assert pow(lam, 6, R_BLS) == R_BLS - 1 and pow(lam, 12, R_BLS) == 1 and pow(lam, 4, R_BLS) != 1
  r_pt, a_pt, k = None, (P0_X, P0_Y), R_BLS
  while k:
    if k & 1: r_pt = _g1a(r_pt, a_pt)
    a_pt = _g1a(a_pt, a_pt); k >>= 1
  assert (P0_Y * P0_Y - (P0_X**3 + 4)) % P_BLS == 0 and r_pt is None
  m2 = lambda a, b: ((a[0] * b[0] - a[1] * b[1]) % P_BLS, (a[0] * b[1] + a[1] * b[0]) % P_BLS)
  q3 = m2(m2(Q0_X, Q0_X), Q0_X); assert m2(Q0_Y, Q0_Y) == ((q3[0] + 4) % P_BLS, (q3[1] + 4) % P_BLS)
  assert P_BLS % 15 == 7 and [j for j in range(1, 5) if pow(P_BLS, j, 5) == 1][0] == 4 and pow(X_BLS % P_BLS, (P_BLS - 1) // 3, P_BLS) != 1
  assert any((X_BLS**3 - X_BLS**2 + (c * c % P_BLS) * X_BLS - c) % P_BLS == 0 for g in range(2, 20) for cand in (pow(g, (P_BLS - 1) // 3, P_BLS),) if cand != 1 for c in (cand, (cand * cand) % P_BLS))
  assert X_BLS**6 - 2 * X_BLS**5 + 2 * X_BLS**3 + X_BLS + 1 == 3 * P_BLS and (R_BLS - 1) % 12 == 0 and ABS_X == (1 << 16) * 906349 * 254760293
  return {"ok": True, "p_mod_15": 7, "ord_5_p": 4, "n_sm_chi1": 4, "abs_x": "2^16*906349*254760293"}

# 2. Unconditional Algebraic Chebotarev Proof (Theorems 2A & 2B)
def check_chebotarev(p_lim: int = 26000) -> dict:
  zwm = lambda x, y: (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0] - x[1] * y[1])
  zwn = lambda x: x[0] * x[0] - x[0] * x[1] + x[1] * x[1]
  zwdiv = lambda d, n: all(c % zwn(d) == 0 for c in zwm(n, (d[0] - d[1], -d[1])))
  dw, dw2, p19, p67 = (41, 24), (17, -24), (-5, -2), (-9, -2)
  assert zwm(p19, p67) == dw and zwn(dw) == 1273 == 19 * 67 and math.gcd(1273, 15) == 1 and not zwdiv(p19, dw2) and not zwdiv(p67, dw2)
  assert (118 * 118 + 118 + 1) % 151 == 0 and all((t**3 - t**2 + 32 * t - 118) % 151 != 0 for t in range(151))
  pw, ww = 991, 113; assert pw % 45 == 1 and (ww * ww + ww + 1) % pw == 0
  wr = [t for t in range(pw) if (t**3 - t**2 + (ww * ww % pw) * t - ww) % pw == 0]; assert wr == [373, 802, 808]
  kch = pow((wr[0] * pow(wr[1], pw - 2, pw)) % pw, (pw - 1) // 3, pw); assert kch == 877 and pow(kch, 3, pw) == 1
  s3 = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
  G = [((a, b, c), s) for a in range(3) for b in range(3) for c in range(3) for s in s3]
  gmul = lambda a, b: (tuple((a[0][b[1][i]] + b[0][i]) % 3 for i in range(3)), tuple(a[1][b[1][i]] for i in range(3)))
  def ginv(g):
    si = [0, 0, 0]
    for i in range(3): si[g[1][i]] = i
    st = tuple(si); return (tuple((-g[0][st[i]]) % 3 for i in range(3)), st)
  seen, ccl = set(), []
  for g in G:
    if g not in seen: cc = frozenset(gmul(gmul(h, g), ginv(h)) for h in G); seen |= cc; ccl.append(cc)
  assert len(G) == 162 and len(ccl) == 22
  ide, hz9 = ((0, 0, 0), (0, 1, 2)), frozenset(g for g in G if sum(g[0]) % 3 == 0)
  nsubs, q = {frozenset([ide])}, [frozenset([ide])]
  while q:
    N = q.pop()
    for cc in ccl:
      if cc.issubset(N): continue
      M = set(N) | set(cc)
      while True:
        add = {gmul(a, b) for a in M for b in cc} - M
        if not add: break
        M |= add
      Mf = frozenset(M)
      if Mf not in nsubs: nsubs.add(Mf); q.append(Mf)
  assert sorted(len(N) for N in nsubs) == [1, 3, 9, 27, 27, 54, 81, 162]
  sh = [N for N in nsubs if N < hz9]; mx = [N for N in sh if not any(N < M for M in sh)]; assert len(mx) == 1 and len(mx[0]) == 27
  rc, c0a, cnza = {0: 0, 3: 0, 6: 0, 9: 0}, 0, 0
  for v, s in G:
    fx = sum(1 for i in range(3) for c in range(3) if s[i] == i and (c + v[i]) % 3 == c); rc[fx] += 1
    if fx:
      if sum(v) % 3 == 0: c0a += 1
      else: cnza += 1
  p0, p3, p6, p9 = (Fraction(rc[k], 162) for k in (0, 3, 6, 9))
  dok, d0, dnz = 1 - p0, Fraction(c0a, 54), Fraction(cnza, 108)
  df = Fraction(1, 8) * (Fraction(1, 3) * (1 - (1 - d0) ** 2) + Fraction(2, 3) * (1 - (1 - dnz) ** 2))
  assert (p0, p3, p6, p9) == (Fraction(58, 81), Fraction(13, 54), Fraction(1, 27), Fraction(1, 162)) and dok == Fraction(23, 81) and df == Fraction(2131, 34992)
  ps = _primes(p_lim); p15 = [p for p in ps if p % 15 == 1]; assert len(p15) > 300; bm = fp = ca = 0
  for p in p15:
    w = next(pow(g, (p - 1) // 3, p) for g in range(2, p) if pow(g, (p - 1) // 3, p) != 1); ph, e3 = False, (p - 1) // 3
    for wv in (w, (w * w) % p):
      w2, ny = (wv * wv) % p, 0
      for t in range(1, p):
        if (t * (t * (t - 1) + w2) - wv) % p == 0 and pow(t, e3, p) == 1: ny += 3
      if ny % 3: bm += 1
      if ny: ca += 1; ph = True
    if ph: fp += 1
  od, odn = ca / (2.0 * len(p15)), fp / float(len(ps)); assert bm == 0 and abs(od - float(dok)) < 0.04 and abs(odn - float(df)) < 0.01
  return {"ok": True, "delta_Q": "41+24*w", "kummer_p": pw, "galois_order": 162, "d_ok_f": str(dok), "dens_f": str(df), "split_p": len(p15), "live_d_ok": round(od, 5), "live_dens": round(odn, 5)}

# 3. Live Algebraic Murphy alpha (p <= 5000), Continuous Norms & Buchstab f-Cap
def _prc(poly, qe, z, o, ad, sb, ml, iv):
  def tr(a):
    while len(a) > 1 and a[-1] == z: a.pop()
    return a
  def pm(a, b):
    r = [z] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
      if x != z:
        for j, y in enumerate(b):
          if y != z: r[i + j] = ad(r[i + j], ml(x, y))
    return tr(r)
  def md(a, m):
    a, dm, il = list(a), len(m) - 1, iv(m[-1])
    while len(a) - 1 >= dm:
      c, sh = ml(a[-1], il), len(a) - 1 - dm
      if c != z:
        for i in range(dm + 1): a[sh + i] = sb(a[sh + i], ml(c, m[i]))
      a.pop()
    return tr(a)
  res, b, e = [o], md([z, o], poly), qe
  while e:
    if e & 1: res = md(pm(res, b), poly)
    b = md(pm(b, b), poly); e >>= 1
  if len(res) < 2: res += [z] * (2 - len(res))
  res[1] = sb(res[1], o); tr(res)
  if res == [z]: return len(poly) - 1
  a, c = list(poly), res
  while any(x != z for x in c):
    a, c = c, md(a, c)
    if c == [z]: break
  return len(a) - 1

def _live_murphy(Pmax: int = 5000) -> tuple[float, float]:
  af = ag = 0.0
  for p in _primes(Pmax):
    er, lp = 1.0 / (p - 1.0), math.log2(p)
    if p == 2: ef, eg = 0.0, 8.0 * (1.0 / 17.0) * (sum(3.0 * k * (15.0 / (16.0**k)) for k in range(1, 6)) + 16.0 / (16.0**5))
    elif p == 3: ef, eg = 4.0 * (3.0 / 82.0), 4.0 * (1.75 / 82.0)
    elif p == 5: ef, eg = 0.0, 2.0 * 3.0 * 25.0 / (25.0**2 - 1.0)
    else:
      cur, fd = p % 15, 1
      while cur != 1: cur = (cur * p) % 15; fd += 1
      gc, qs = 8 // fd, p**fd; pf, xm = float(qs), X_BLS % p
      rg = 1 if xm == 0 else ((3 if pow(xm, (qs - 1) // 3, p) == 1 else 0) if (qs - 1) % 3 == 0 else 1)
      eg = fd * gc * rg * pf / (pf * pf - 1.0)
      ops = (0, 1, lambda a, b: (a + b) % p, lambda a, b: (a - b) % p, lambda a, b: (a * b) % p, lambda a: pow(a, p - 2, p))
      if fd == 1:
        w1 = next(pow(g, (p - 1) // 3, p) for g in range(2, p) if pow(g, (p - 1) // 3, p) != 1)
        ef = sum(4 * _prc([(-w) % p, 0, 0, (w * w) % p, 0, 0, (-1) % p, 0, 0, 1], p, *ops) for w in (w1, (w1 * w1) % p)) * pf / (pf * pf - 1.0)
      elif (fd == 2 and p <= 300) or (fd == 4 and p <= 13):
        if p % 3 == 1:
          w1 = next(pow(g, (p - 1) // 3, p) for g in range(2, p) if pow(g, (p - 1) // 3, p) != 1)
          rf = sum((4 // fd) * _prc([(-w) % p, 0, 0, (w * w) % p, 0, 0, (-1) % p, 0, 0, 1], qs, *ops) for w in (w1, (w1 * w1) % p))
        else:
          f2i = lambda a: (((a[0] - a[1]) * pow((a[0] * a[0] - a[0] * a[1] + a[1] * a[1]) % p, p - 2, p)) % p, ((-a[1]) * pow((a[0] * a[0] - a[0] * a[1] + a[1] * a[1]) % p, p - 2, p)) % p)
          pf2 = [(0, (-1) % p), (0, 0), (0, 0), ((-1) % p, (-1) % p), (0, 0), (0, 0), ((-1) % p, 0), (0, 0), (0, 0), (1, 0)]
          rf = gc * _prc(pf2, qs, (0, 0), (1, 0), lambda a, b: ((a[0] + b[0]) % p, (a[1] + b[1]) % p), lambda a, b: ((a[0] - b[0]) % p, (a[1] - b[1]) % p), lambda a, b: ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0] - a[1] * b[1]) % p), f2i)
        ef = fd * rf * pf / (pf * pf - 1.0)
      else: ef = 8.0 / pf
    af += (er - ef) * lp; ag += (er - eg) * lp
  return af, ag

def _cof_ex(Bf: float, pfhi: float, mnf: float, lf: float, af: float) -> tuple[float, float]:
  M, dx, h = 80, (Bf - pfhi) / 80.0, 2.0; nb = int(math.ceil(660.0 / h)) + 1; cur = [0.0] * nb; cur[0] = 1.0; tp = ep = 0.0; ne = mnf + af
  for k in range(1, 11):
    nxt = [0.0] * nb
    for m in range(M):
      x = pfhi + (m + 0.5) * dx; w = dx / (x * k); sh = int(round(x / h))
      for b in range(nb - sh):
        if cur[b] > 0.0: nxt[b + sh] += cur[b] * w
    cur = nxt
    for b in range(nb):
      sb = b * h
      if sb > ne: break
      cv = cur[b] * DK.rho((ne - sb) / pfhi); tp += cv
      if sb > lf: ep += cv
  tp += DK.rho(ne / pfhi); return ep / tp, tp / DK.rho(ne / Bf)

def check_murphy_norms() -> dict:
  af, ag = _live_murphy(5000); assert abs(af - REF["af"]) < 0.002 and abs(ag - REF["ag"]) < 0.002
  mc, b = REF["mc"], REF["b"]; R0, Ra, q = mc["R0"], mc["Ra"], REF["qa"]
  dnf = (72.0 / q) * math.log2((Ra**q + R0**q) / 4.0); dng = 8.0 * LOG2_ABS_X + (24.0 / q) * math.log2((R0**q + (Ra**q) / (ABS_X ** (q / 3.0))) / 4.0)
  assert abs(dnf - mc["dnf"]) < 1e-9 and dnf < 656.0 and abs(dng - mc["dng"]) < 1e-9 and dng < 688.0 and 6 * 656 <= 4096 and 5 * 688 <= 4096 and 10 * 408 <= 4096
  tfm = int(math.floor(656.0 / b["qf"])) - 127; assert tfm == 45 and (128 - tfm) + int(math.floor(656.0 / b["qf"])) <= 255
  for ngc in range(610, 689): assert (128 - max(1, int(math.floor((ngc + 1.0) / b["qg"])) - 127)) + int(math.floor(ngc / b["qg"])) <= 255
  w0 = 1.0 / (1.0 + REF["rbr"] ** 8); p0, bn0 = _cof_ex(REF["Bf"], b["pfhi0"], mc["mnf"], REF["lf"], REF["af"]); p1, bn1 = _cof_ex(REF["Bf"], b["pfhi1"], mc["mnf"], REF["lf"], REF["af"])
  assert abs(bn0 - 1.0) < 0.01 and abs(bn1 - 1.0) < 0.01; pex = w0 * p0 + (1.0 - w0) * p1; cpad = math.log2(1.0 / (1.0 - pex))
  assert abs(pex - mc["pex"]) < 1e-9 and abs(cpad - mc["cpad"]) < 1e-9
  return {"ok": True, "live_af": round(af, 5), "live_ag": round(ag, 5), "det_nf": round(dnf, 6), "det_ng": round(dng, 6), "p_ex_f": round(pex, 6), "cof_pad": round(cpad, 6), "bnorm": [round(bn0, 5), round(bn1, 5)]}

# 4. Live 2D Panjer DP + Policy Solver + Out-of-Sample 16D MC + Live Chords
VH = ((0.481211825060, 0.481211825060, -0.481211825060, -0.481211825060), (0.687654016088, -0.206442191029, 0.396282847190, -0.877494672250), (0.206442191029, -0.687654016088, 0.877494672250, -0.396282847190), (0.877494672250, -0.396282847190, -0.687654016088, 0.206442191029), (0.396282847190, -0.877494672250, -0.206442191029, 0.687654016088), (1.083936863278, -1.083936863278, 0.189840656161, -0.189840656161), (0.189840656161, -0.189840656161, -1.083936863278, 1.083936863278))
VHH = tuple(0.5 * sum(c * c for c in v) for v in VH)

def _panjer_dp(Bb: float, plo: float, phi: float, qs: float, ab: float, bmin: int, bstep: int, nb: int):
  KM, TM = 720, 260; mo = [0.0] * (KM + 1)
  for a_i, b_i in ((0.561459483567, plo), (phi, Bb)):
    if b_i <= a_i: continue
    ns = max(20, int((b_i - a_i) * 8.0)); dx = (b_i - a_i) / ns
    for i in range(ns):
      x0 = a_i + i * dx; kb = int(round(x0 + 0.5 * dx))
      if 1 <= kb <= KM: mo[kb] += math.log((x0 + dx) / x0)
  nzo = [(j, j * mo[j]) for j in range(1, KM + 1) if mo[j] > 0.0]; Po = [0.0] * (KM + 1); Po[0] = math.exp(-sum(mo[1:]))
  for k in range(1, KM + 1): Po[k] = sum(jm * Po[k - j] for j, jm in nzo if j <= k) / k
  mi, m1 = {}, [0.0] * (TM + 1); ns = max(80, int((phi - plo) * 16.0)); dx = (phi - plo) / ns
  for i in range(ns):
    x0 = plo + i * dx; xm = x0 + 0.5 * dx; w, kb, tk = math.log((x0 + dx) / x0), int(round(xm)), int(math.floor(xm / qs))
    if kb >= 1 and 1 <= tk <= TM: mi[(kb, tk)] = mi.get((kb, tk), 0.0) + w; m1[tk] += w
  li = sum(m1[1:]); nz1 = [(j, j * m1[j]) for j in range(1, TM + 1) if m1[j] > 0.0]; P1 = [0.0] * (TM + 1); P1[0] = math.exp(-li)
  for k in range(1, TM + 1): P1[k] = sum(jm * P1[k - j] for j, jm in nz1 if j <= k) / k
  sv, run = [0.0] * (TM + 1), 0.0
  for t in range(TM, -1, -1): run += P1[t]; sv[t] = run
  its = sorted([(dk, dt, dk * v) for (dk, dt), v in mi.items()]); P2 = [{} for _ in range(KM + 1)]; P2[0][0] = math.exp(-li)
  for k in range(its[0][0], KM + 1):
    ik, dst = 1.0 / k, {}
    for dk, dt, dv in its:
      if dk > k: break
      pr = P2[k - dk]
      if not pr: continue
      sd = dv * ik
      for tt, p0 in pr.items():
        if tt + dt <= TM: dst[tt + dt] = dst.get(tt + dt, 0.0) + sd * p0
    if dst: P2[k] = {tt: pv for tt, pv in dst.items() if pv > 1e-13}
  nzk = [(k, list(P2[k].items())) for k in range(KM + 1) if P2[k]]; rt = []
  for bi in range(nb):
    kN = max(1, min(KM, int(round(bmin + (bi + 0.5) * bstep + ab)))); sm = [0.0] * (TM + 1)
    for k, ki in nzk:
      if k > kN: break
      po = Po[kN - k]
      if po <= 1e-300: continue
      for tt, pt in ki: sm[tt] += pt * po
    tot, tl, row = sum(sm), 0.0, [0.0] * 129
    for tt in range(TM, -1, -1):
      tl += sm[tt]
      if tt <= 128: row[tt] = (tl / tot) if tot > 0.0 else 0.0
    rt.append(row)
  return sv, rt

def check_live_panjer_mc_chords(nfit: int = 5000, neval: int = 5000, nch: int = 1000) -> dict:
  af, ag = _live_murphy(5000); b, mc, pol = REF["b"], REF["mc"], REF["pol"]
  Bf, Bg, qa, rbr, Ra, R0 = REF["Bf"], REF["Bg"], REF["qa"], REF["rbr"], mc["Ra"], mc["R0"]
  w0 = 1.0 / (1.0 + rbr**8); w1 = 1.0 - w0; G0, GS, GB, F0, FS, FB, N0, NS, NB = 610, 2, 40, 510, 4, 25, 440, 2, 100
  sgL, rgL, sfL, rfL = [], [], [], []
  for pgh, pfh in ((b["pghi0"], b["pfhi0"]), (b["pghi1"], b["pfhi1"])):
    sg, rg = _panjer_dp(Bg, b["pglo"], pgh, b["qg"], ag, G0, GS, GB); sf, rf = _panjer_dp(Bf, b["pflo"], pfh, b["qf"], af, N0, NS, NB)
    sgL.append(sg); rgL.append(rg); sfL.append(sf); rfL.append(rf)
  iqa, agd, il2, tp = 1.0 / qa, 2.0 / qa, 1.0 / math.log(2.0), 2.0 * math.pi
  W, W2 = complex(-0.5, 0.8660254037844386), complex(-0.5, -0.8660254037844386)
  Ap = [(a**6, a**12, a**18) for a in ((Ra / (4.0**iqa)) * ((0.25 if k == 0 else 0.75) ** 0.125) for k in range(2))]
  def batch(nd: int, sd: int):
    rng, n0, vt, vo, pts = random.Random(sd), int(round(nd * w0)), 0, 0, []
    for i in range(nd):
      sh = 0 if i < n0 else 1; r0s, r1s, Rb, (A6, A12, A18) = (0.0, 0.5, R0, Ap[0]) if sh == 0 else (0.5, 1.0, R0 * rbr, Ap[1])
      while True:
        vt += 1; ga = [rng.gammavariate(agd, 1.0) for _ in range(4)]; ta = sum(ga)
        u = [iqa * math.log(g / ta) for g in ga]; um = 0.25 * sum(u); u = [x - um for x in u]
        if all(abs(u[0] * v[0] + u[1] * v[1] + u[2] * v[2] + u[3] * v[3]) <= h for v, h in zip(VH, VHH)): vo += 1; break
      rqa = (Ra * ((r0s + (r1s - r0s) * rng.random()) ** 0.125)) ** qa; gb = [rng.gammavariate(agd, 1.0) for _ in range(4)]; tb = sum(gb)
      rqb = (Rb * (max(rng.random(), 1e-16) ** 0.125)) ** qa; lnf = lng = lng0 = nfp = 0.0
      for m in range(4):
        rm, tm = (rqa * ga[m] / ta) ** iqa, tp * rng.random(); az = complex(rm * math.cos(tm), rm * math.sin(tm))
        bm, tb2 = (rqb * gb[m] / tb) ** iqa, tp * rng.random(); bz = complex(bm * math.cos(tb2), bm * math.sin(tb2))
        r6 = bm**6; r12 = r6 * r6; nfp += math.log(A18 + A12 * r6 + A6 * r12 + r12 * r6) * il2
        a3, b3 = az * az * az, bz * bz * bz; gz, g0z = a3 - X_BLS * b3, -X_BLS * b3
        fz = a3 * a3 * a3 - a3 * a3 * b3 + W2 * a3 * b3 * b3 - W * b3 * b3 * b3
        lng += math.log(gz.real * gz.real + gz.imag * gz.imag + 1e-300); lng0 += math.log(g0z.real * g0z.real + g0z.imag * g0z.imag + 1e-300)
        lnf += math.log(fz.real * fz.real + fz.imag * fz.imag + 1e-300)
      nf, ng, ng0 = lnf * il2, lng * il2, lng0 * il2
      pts.append((sh, nf, ng, ng0, nfp, DK.rho((nf + af) / Bf) * DK.rho((ng + ag) / Bg)))
    return pts, vo / float(vt)
  fpts, fvf = batch(nfit, 20261127); v2 = [[[0.0] * FB for _ in range(GB)] for _ in range(2)]; s3 = [[[[0.0] * NB for _ in range(FB)] for _ in range(GB)] for _ in range(2)]; tsf = sum(p[5] for p in fpts)
  for sh, nf, ng, ng0, nfp, wt in fpts:
    ig, icf, iff = max(0, min(GB - 1, int(math.floor((ng0 - G0) / GS)))), max(0, min(FB - 1, int(math.floor((nfp - F0) / FS)))), max(0, min(NB - 1, int(math.floor((nf - N0) / NS))))
    v2[sh][ig][icf] += 1.0 / nfit; s3[sh][ig][icf][iff] += wt / tsf
  tfm, lp = max(0, int(math.floor(656.0 / b["qf"])) - 127), pol["lam"]
  Tg, Tf = [[[128] * FB for _ in range(GB)] for _ in range(2)], [[[125] * FB for _ in range(GB)] for _ in range(2)]
  for sh in range(2):
    rgs, rfs, sgs, sfs = rgL[sh], rfL[sh], sgL[sh], sfL[sh]
    for ig in range(GB):
      ngu = G0 + (ig + 1) * GS + 1.0; tgm = max(0, min(128, int(math.floor(ngu / b["qg"])) - 127)); vw1, rf1 = sum(v2[sh][ig]), [0.0] * 126
      if vw1 > 0.0:
        for tf in range(tfm, 126): rf1[tf] = sum(s3[sh][ig][icf][iff] * rfs[iff][tf] for icf in range(FB) for iff in range(NB) if s3[sh][ig][icf][iff] > 0.0)
      for icf in range(FB):
        vw = v2[sh][ig][icf]
        if vw < 5.0 / nfit:
          if vw1 <= 0.0: continue
          vwu, rfw = vw1, rf1
        else:
          vwu, rfw = vw, [0.0] * 126; nzf = [(iff, s3[sh][ig][icf][iff]) for iff in range(NB) if s3[sh][ig][icf][iff] > 0.0]
          for tf in range(tfm, 126): rfw[tf] = sum(sv * rfs[iff][tf] for iff, sv in nzf)
        bo, btg, btf = -1e300, 128, 125
        for tg in range(tgm, 129):
          rgv, sgt = rgs[ig][tg], sgs[tg]
          for tf in range(tfm, 126):
            o = rgv * rfw[tf] - lp * vwu * sgt * sfs[tf]
            if o > bo: bo, btg, btf = o, tg, tf
        Tg[sh][ig][icf], Tf[sh][ig][icf] = btg, btf
        assert (128 - btg) + int(math.floor(ngu / b["qg"])) <= 255 and (128 - btf) + int(math.floor(656.0 / b["qf"])) <= 255
  epts, fve = batch(neval, 20261128); swr = swk = swk2 = sgs_ = s1s_ = nfm = ngm = 0.0
  for sh, nf, ng, ng0, nfp, wt in epts:
    if nf > nfm: nfm = nf
    if ng > ngm: ngm = ng
    ig, icf = max(0, min(GB - 1, int(math.floor((ng0 - G0) / GS)))), max(0, min(FB - 1, int(math.floor((nfp - F0) / FS))))
    tg, tf = Tg[sh][ig][icf], Tf[sh][ig][icf]
    wk = wt * rgL[sh][max(0, min(GB - 1, int(math.floor((ng - G0) / GS))))][tg] * rfL[sh][max(0, min(NB - 1, int(math.floor((nf - N0) / NS))))][tf]
    swr += wt; swk += wk; swk2 += wk * wk; sp = sgL[sh][tg]; sgs_ += sp; s1s_ += sp * sfL[sh][tf]
  retL, l2sg, l2s1 = swk / swr, math.log2(sgs_ / neval), math.log2(s1s_ / neval); l2sf = l2s1 - l2sg; mwk = swk / neval
  seb = (math.sqrt(max(1e-300, (swk2 / neval - mwk * mwk) / neval)) / mwk) * il2
  Rrl = (2.0 ** mc["Vk4"]) * mwk * (0.97843406 * (2.0 ** REF["cush"])); VkL, fvL = mc["Vk4"] + mc["sl"] * (mc["Rreq"] - math.log2(Rrl)), 0.5 * (fvf + fve)
  assert nfm <= mc["dnf"] and ngm <= mc["dng"] and abs(fvL - mc["fv"]) < 0.005 and abs(retL - pol["ret"]) < 0.03 and abs(l2sg - pol["sg"]) < 0.50 and abs(l2s1 - pol["s1"]) < 0.60 and abs(VkL - mc["Vk4"]) < 3.0 * (mc["sl"] * seb)
  rc_, L4L, g4L, hq = random.Random(20260505), 0.0, 0.0, 0.5 * qa
  for R, wsh in ((R0, w0), (R0 * rbr, w1)):
    Rq, iL, ssp, sgp, w248 = R**qa, 0.0, 0.0, 0.0, 0.0; ns, tl = max(200, int(round(nch * wsh))), int(math.ceil(R * 1.15)) + 20
    for _ in range(ns):
      gs = [rc_.gammavariate(agd, 1.0) for _ in range(4)]; gsm = sum(gs); rq = Rq * (max(rc_.random(), 1e-16) ** (qa / 8.0))
      zr, zi = [], []
      for m in range(4): rm, th = (rq * gs[m] / gsm) ** iqa, tp * rc_.random(); zr.append(rm * math.cos(th)); zi.append(rm * math.sin(th))
      ch = [tc for tc in range(-tl, tl + 1, 16) if sum(((zr[m] + tc) ** 2 + zi[m] ** 2) ** hq for m in range(4)) <= Rq]
      tlo, thi = max(-tl, (ch[0] - 32) if ch else -256), min(tl, (ch[-1] + 32) if ch else 256); tmn, tmx, ic = 10000, -10000, 0
      for t in range(tlo, thi + 1):
        if sum(((zr[m] + t) ** 2 + zi[m] ** 2) ** hq for m in range(4)) <= Rq:
          if t < tmn: tmn = t
          if t > tmx: tmx = t
          ic += 1
      if ic: sp_ = tmx - tmn + 1; w = 1.0 / ic; iL += w; ssp += w * sp_; sgp += w * (sp_ - ic); w248 += w * ((sp_ + 247) // 248)
    L4L += wsh * ((ssp / iL) / (w248 / iL)); g4L += wsh * (sgp / ssp)
  assert abs(L4L - mc["L4"]) < 6.0 and 0.0 <= g4L < 0.006
  return {"ok": True, "n_fit": nfit, "n_eval": neval, "n_chords": nch, "fvor_live": round(fvL, 5), "ret_live": round(retL, 4), "ret_4m": round(pol["ret"], 4), "l2_sg": round(l2sg, 3), "l2_sf": round(l2sf, 3), "l2_s1": round(l2s1, 3), "Vk_oos": round(VkL, 4), "Vk_4m": round(mc["Vk4"], 4), "se_Vk": round(mc["sl"] * seb, 4), "L4_live": round(L4L, 2), "gap4_live": round(g4L, 5)}

# 5. Live Bit-Exact 3-Prime 2^54-Smooth Radix-4 u2 Shoup61 NTT & 2D Mul
def _ip64(n: int) -> bool:
  d, s = n - 1, 0
  while d % 2 == 0: d //= 2; s += 1
  for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
    x = pow(a, d, n)
    if x in (1, n - 1): continue
    for _ in range(s - 1):
      x = (x * x) % n
      if x == n - 1: break
    else: return False
  return True

def _shoup(a: int, w: int, wb: int, p: int) -> int: return (a * w - ((a * wb) >> 64) * p) & 0xFFFFFFFFFFFFFFFF

def _ntt_r4(vec: list[int], p: int, g: int) -> list[int]:
  N, a, p2, p4, l4, t = len(vec), [x % p for x in vec], p << 1, p << 2, 0, len(vec)
  while t > 1: t >>= 2; l4 += 1
  for i in range(N):
    rv, x = 0, i
    for _ in range(l4): rv = (rv << 2) | (x & 3); x >>= 2
    if i < rv: a[i], a[rv] = a[rv], a[i]
  J, Jp, m = pow(g, (p - 1) // 4, p), (pow(g, (p - 1) // 4, p) << 64) // p, 1
  while m < N:
    wm, w1, W = pow(g, (p - 1) // (4 * m), p), 1, []
    for _ in range(m): w2 = (w1 * w1) % p; w3 = (w2 * w1) % p; W.append((w1, (w1 << 64) // p, w2, (w2 << 64) // p, w3, (w3 << 64) // p)); w1 = (w1 * wm) % p
    for k in range(0, N, 4 * m):
      for j in range(m):
        w1, w1p, w2, w2p, w3, w3p = W[j]; x0 = a[k + j] - (p2 if a[k + j] >= p2 else 0)
        x1, x2, x3 = _shoup(a[k + m + j], w1, w1p, p), _shoup(a[k + 2 * m + j], w2, w2p, p), _shoup(a[k + 3 * m + j], w3, w3p, p)
        t0, t1 = (x0 + x2) - (p2 if x0 + x2 >= p2 else 0), (x0 + p2 - x2) - (p2 if x0 + p2 - x2 >= p2 else 0)
        t2, t3 = (x1 + x3) - (p2 if x1 + x3 >= p2 else 0), _shoup(x1 + p2 - x3, J, Jp, p)
        a[k + j], a[k + 2 * m + j], a[k + m + j], a[k + 3 * m + j] = t0 + t2, t0 + p2 - t2, t1 + t3, t1 + p2 - t3
        assert a[k + j] < p4 and a[k + m + j] < p4 and a[k + 2 * m + j] < p4 and a[k + 3 * m + j] < p4
    m <<= 2
  return [x % p for x in a]

def _ntt1(vec: list[int], p: int, g: int, inv: bool = False) -> list[int]:
  n, a, p2, j = len(vec), [x % p for x in vec], 2 * p, 0
  for i in range(1, n):
    b = n >> 1
    while j & b: j ^= b; b >>= 1
    j ^= b
    if i < j: a[i], a[j] = a[j], a[i]
  L = 2
  while L <= n:
    hf, wl = L >> 1, pow(g, (p - 1) // L, p)
    if inv: wl = pow(wl, p - 2, p)
    ws, wp, cur = [1] * hf, [0] * hf, 1
    for k in range(hf): ws[k], wp[k] = cur, (cur << 64) // p; cur = (cur * wl) % p
    for i in range(0, n, L):
      for k in range(hf): u = a[i + k] - (p2 if a[i + k] >= p2 else 0); t = _shoup(a[i + k + hf], ws[k], wp[k], p); a[i + k], a[i + k + hf] = u + t, u + p2 - t
    L <<= 1
  return [(x * pow(n, p - 2, p)) % p for x in a] if inv else [x % p for x in a]

def _mul2d(A: int, B: int, pgs: list[tuple[int, int]]) -> int:
  m58 = (1 << 58) - 1
  def t58(x):
    r = []
    while x: r.append(x & m58); x >>= 58
    return r or [0]
  la, lb = t58(A), t58(B); k = max(1, (max(len(la), len(lb)) - 1).bit_length()); k1 = (k + 1) // 2; D1, D2 = 1 << k1, 1 << (k - k1); N1, N2 = 2 * D1, 2 * D2
  la += [0] * (D1 * D2 - len(la)); lb += [0] * (D1 * D2 - len(lb)); cm = []
  for p, g in pgs:
    gA = [_ntt1(la[r * D1 : (r + 1) * D1] + [0] * D1, p, g) for r in range(D2)] + [[0] * N1 for _ in range(D2)]
    gB = [_ntt1(lb[r * D1 : (r + 1) * D1] + [0] * D1, p, g) for r in range(D2)] + [[0] * N1 for _ in range(D2)]
    gC = [[0] * N1 for _ in range(N2)]
    for c in range(N1):
      cA, cB = _ntt1([gA[r][c] for r in range(N2)], p, g), _ntt1([gB[r][c] for r in range(N2)], p, g)
      cC = _ntt1([(cA[r] * cB[r]) % p for r in range(N2)], p, g, True)
      for r in range(N2): gC[r][c] = cC[r]
    cm.append([_ntt1(gC[r], p, g, True) for r in range(N2)])
  p1, p2, p3 = pgs[0][0], pgs[1][0], pgs[2][0]; i12, i123, p12 = pow(p1 % p2, p2 - 2, p2), pow((p1 * p2) % p3, p3 - 2, p3), p1 * p2; acc = [0] * (D1 * (N2 + 1) + 4)
  for r in range(N2):
    bs = r * D1
    for c in range(N1): a1, a2, a3 = cm[0][r][c], cm[1][r][c], cm[2][r][c]; x12 = a1 + p1 * (((a2 - a1) * i12) % p2); acc[bs + c] += x12 + p12 * (((a3 - x12) * i123) % p3)
  res = cy = 0
  for i, v in enumerate(acc): s = v + cy; res |= (s & m58) << (58 * i); cy = s >> 58
  return res | (cy << (58 * len(acc)))

def check_ntt() -> dict:
  s = REF["s1b"]; p1, p2, p3, g1, g2, g3 = s["p1"], s["p2"], s["p3"], s["g1"], s["g2"], s["g3"]
  assert all(_ip64(p) and (p - 1) % (1 << 54) == 0 and 8 * p < (1 << 64) and (1 << 58) < p for p in (p1, p2, p3))
  assert all(pow(g1, (p1 - 1) // q, p1) != 1 for q in (2, 127)) and all(pow(g2, (p2 - 1) // q, p2) != 1 for q in (2, 3, 19)) and all(pow(g3, (p3 - 1) // q, p3) != 1 for q in (2, 3))
  l2P = math.log2(p1 * p2 * p3); h1, h2 = l2P - 181.0, l2P - (math.log2(math.ceil((2.0**67.19) / 58.0)) + 116.0); assert h1 > 1.57 and h2 > 5.24
  rng, pgs = random.Random(0x5400612D), [(p1, g1), (p2, g2), (p3, g3)]
  for p, g in pgs:
    N = 64; vec = [rng.randrange(p) for _ in range(N)]; got = _ntt_r4(vec, p, g); wN = pow(g, (p - 1) // N, p)
    for k in range(N):
      wk, acc, wj = pow(wN, k, p), 0, 1
      for j in range(N): acc = (acc + vec[j] * wj) % p; wj = (wj * wk) % p
      assert got[k] == acc
  sz = [512, 1024, 2048, 4096]
  for nb in sz: A, B = rng.getrandbits(nb) | (1 << (nb - 1)), rng.getrandbits(nb) | (1 << (nb - 1)); assert _mul2d(A, B, pgs) == A * B
  return {"ok": True, "r4_bitexact": True, "log2_P": round(l2P, 6), "hr_1d": round(h1, 4), "hr_2d": round(h2, 4), "mul_2d": sz}

# 6. End-to-End 32-Col Tier-2 W4, Tier-3 W16 Zero-AND, 128B SWAR & 576-Bit Horner+SM4 Audit
def _horner576(A: list[int], sm: list[int], v: list[int]) -> int:
  r256, r2_256, x2h, m64 = R_BLS << 256, (2 * R_BLS) << 256, 0xAC45A4010001A402, (1 << 64) - 1
  T = A[0] + (A[2] << 32) + ((A[2] * x2h) << 64) + r256 - ABS_X * (A[1] + (A[3] << 32) + ((A[3] * x2h) << 64)) + sum(sm[k] * v[k] for k in range(4))
  assert T >= 0; t8 = T >> 512; assert t8 in (0, 1)
  if t8: T -= r2_256
  assert 0 <= T < (1 << 512); Th = T >> 256; t4, t5, t6, t7 = Th & m64, (Th >> 64) & m64, (Th >> 128) & m64, (Th >> 192) & m64
  M = [0x42737A020C0D6393, 0x65043EB4BE4BAD71, 0x38B5DCB707E08ED3, 0x355094EDFEDE377C]
  ph = ((t4 * M[3]) + (t5 * M[2] + ((t5 * M[3]) << 64)) + (t6 * M[1] + ((t6 * M[2]) << 64) + ((t6 * M[3]) << 128)) + (t7 * M[0] + ((t7 * M[1]) << 64) + ((t7 * M[2]) << 128) + ((t7 * M[3]) << 192))) >> 64
  s = (T - ((Th << 1) + ph) * R_BLS) & ((1 << 320) - 1); assert 0 <= s < 6 * R_BLS; s2 = s - (s >> 255) * R_BLS; assert 0 <= s2 < 2 * R_BLS
  return s2 - R_BLS if s2 >= R_BLS else s2

def check_kernels_e2e() -> dict:
  pm, bi, st, nh = 31, 5, 13, 200; cc, bc = [0] * 32, bi
  for _ in range(nh): bc = (bc + st) % pm; cc[bc] += 1
  sw4, bl4, cur = [], [0] * 8, [0] * 32
  for j in range(8):
    n0, n1, n2, n3 = cc[4 * j : 4 * j + 4]; tot = n0 + n1 + n2 + n3; bl4[j] = tot
    if not tot: continue
    sw4.append(n0 | (n1 << 16) | (n2 << 32) | (n3 << 48))
    for c in range(4): cur[4 * j + c] = len(sw4); sw4.extend([0] * cc[4 * j + c])
  bc = bi
  for aid in range(nh): bc = (bc + st) % pm; sw4[cur[bc]] = SLAB_BASE + aid; cur[bc] += 1
  ow4, pos = [[] for _ in range(32)], 0
  for j in range(8):
    if not bl4[j]: continue
    H = sw4[pos]; pos += 1
    for c in range(4):
      nc = (H >> (16 * c)) & 0xFFFF
      if nc: ow4[4 * j + c].extend(sw4[pos : pos + nc]); pos += nc
  assert pos == len(sw4); rw4, bc = [[] for _ in range(32)], bi
  for aid in range(nh): bc = (bc + st) % pm; rw4[bc].append(SLAB_BASE + aid)
  assert ow4 == rw4
  Va = 1024; sl, rf = [bytearray([120] * Va) for _ in range(32)], [bytearray([120] * Va) for _ in range(32)]; s16, bl16 = [], [257, 256]
  for j in range(2):
    for k in range(bl16[j]):
      lc, off = (k * 7 + j * 3) & 15, (k * 97 + j * 13) % Va; p64 = SLAB_BASE + off; assert (p64 >> 57) == 0
      ent = (lc << 60) | p64; assert (ent >> 57) == (lc << 3); s16.append(ent); rf[j * 16 + lc][off] = (rf[j * 16 + lc][off] + 12) & 0xFF
  rx16, pos = [[] for _ in range(32)], 0
  for j in range(2):
    for _ in range(bl16[j]): e = s16[pos]; pos += 1; rx16[j * 16 + ((e >> 57) >> 3)].append(e & 0x0FFFFFFFFFFFFFFF)
  sv16 = 0
  for c in range(32):
    for ptr in rx16[c]: aid = ptr - SLAB_BASE; sl[c][aid] = (sl[c][aid] + 12) & 0xFF
    assert sl[c] == rf[c]; gi = []
    for blk in range(0, Va, 128):
      ck = sl[c][blk : blk + 128]; m = 0
      for i in range(0, 128, 8): m |= int.from_bytes(ck[i : i + 8], "little")
      if m & 0x8080808080808080:
        for kw in range(16):
          bt = int.from_bytes(ck[kw * 8 : kw * 8 + 8], "little") & 0x8080808080808080
          while bt: gi.append(blk + (kw << 3) + (((bt & -bt).bit_length() - 1) >> 3)); bt &= bt - 1
    assert gi == [idx for idx in range(Va) if rf[c][idx] & 0x80]; sv16 += len(gi)
  rng = random.Random(0x576B1750)
  for t in range(60):
    A, sm, vv = ([(1 << 304) - 1, 0, (1 << 304) - 1, 0], [R_BLS - 1] * 4, [R_BLS - 1] * 4) if t == 0 else (([0, (1 << 314) - 1, 0, (1 << 314) - 1], [0] * 4, [0] * 4) if t == 1 else (([(1 << 300) - 1] * 4, [R_BLS - 1] * 4, [R_BLS - 1] * 4) if t == 2 else ([rng.getrandbits(314) for _ in range(4)], [rng.randrange(R_BLS) for _ in range(4)], [rng.randrange(R_BLS) for _ in range(4)])))
    assert _horner576(A, sm, vv) == (A[0] + X_BLS * A[1] + X_BLS * X_BLS * A[2] + (X_BLS**3) * A[3] + sum(sm[k] * vv[k] for k in range(4))) % R_BLS
  Vv, r60 = [rng.randrange(R_BLS) for _ in range(256)], R_BLS << 60
  for _ in range(30):
    Ap, Am, rr = [0, 0, 0, 0], [0, 0, 0, 0], 0
    for _ in range(120):
      col, k, sg = rng.randrange(256), rng.randrange(4), (1 if rng.random() < 0.5 else -1); cm = 1 if rng.random() < 0.8143 else rng.randrange(2, 1 << 30)
      if sg > 0: Ap[k] += cm * Vv[col]
      else: Am[k] += cm * Vv[col]
      rr = (rr + sg * cm * pow(X_BLS, k, R_BLS) * Vv[col]) % R_BLS
    As = [Ap[k] + r60 - Am[k] for k in range(4)]; assert all(0 <= Ap[k] < (1 << 293) and 0 <= Am[k] < (1 << 293) and 0 <= As[k] < (1 << 315) for k in range(4))
    sm, vv = [rng.randrange(R_BLS) for _ in range(4)], [rng.randrange(R_BLS) for _ in range(4)]
    assert _horner576(As, sm, vv) == (rr + sum(sm[j] * vv[j] for j in range(4))) % R_BLS
  return {"ok": True, "w4_hits": nh, "w16_hits": sum(bl16), "swar_surv": sv16, "horner576": "60/60", "signed_spmv": "30/30"}

# 7. Faithful 4-Tuple Z[zeta_12] STRICT Monomial-Pivot SGE + Exact BLS12-381 F_r Rank (64 Live Seeds)
MONO = ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (-1, 0, 1, 0), (0, -1, 0, 1), (-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1), (1, 0, -1, 0), (0, 1, 0, -1))
MMAP = {v: i for i, v in enumerate(MONO)}

def _zmulm(a, k: int):
  k %= 12; ng = k >= 6
  if ng: k -= 6
  c0, c1, c2, c3 = a
  for _ in range(k): c0, c1, c2, c3 = -c3, c0, c1 + c3, c2
  return (-c0, -c1, -c2, -c3) if ng else (c0, c1, c2, c3)

def _z12_sge_rank(U: int, eta: float, sd: int) -> dict:
  rng = random.Random(sd); Uf, lBf, lBg = int(round(0.4392 * U)), 70.80 * math.log(2.0), 71.16 * math.log(2.0)
  Ug, rf1, rg1 = U - Uf, DK.rho(6.65), DK.rho(8.31)
  ws = [min(0.45, (1.173 / U) / ((j + 0.5) / Uf) * (DK.rho(7.65 - max(0.02, 1.0 + math.log((j + 0.5) / Uf) / lBf)) / rf1) + (0.35 / (2.0 + 2.5 * j) if j < 40 else 0.0)) for j in range(Uf)] + [min(0.45, (1.196 / U) / ((j + 0.5) / Ug) * (DK.rho(9.31 - max(0.02, 1.0 + math.log((j + 0.5) / Ug) / lBg)) / rg1) + (0.40 / (2.0 + 2.5 * j) if j < 40 else 0.0)) for j in range(Ug)]
  mw, cdf, acc = sum(ws), [], 0.0
  for w in ws: acc += w / mw; cdf.append(acc)
  rows, Lp = [], math.exp(-mw)
  for _ in range(int(round(eta * U))):
    kp, pa = 0, 1.0
    while pa > Lp: kp += 1; pa *= rng.random()
    d = {}
    for _ in range(max(4, kp - 1)):
      c, m = min(U - 1, bisect.bisect_left(cdf, rng.random())), MONO[rng.randrange(12)]
      if c in d:
        a = d[c]; nv = (a[0] + m[0], a[1] + m[1], a[2] + m[2], a[3] + m[3])
        if nv == (0, 0, 0, 0): del d[c]
        else: d[c] = nv
      else: d[c] = m
    if d: rows.append(d)
  cr = {}
  for ri, d in enumerate(rows):
    for c in d: cr.setdefault(c, set()).add(ri)
  def kr(ri):
    d = rows[ri]
    if d is None: return
    for c in d:
      s = cr.get(c)
      if s:
        s.discard(ri)
        if not s: cr.pop(c, None)
    rows[ri] = None
  def psg():
    q = [c for c, s in cr.items() if len(s) == 1]
    while q:
      c = q.pop(); s = cr.get(c)
      if not s or len(s) != 1: continue
      ri = next(iter(s)); d = rows[ri]
      if d is None: continue
      nb = list(d.keys()); kr(ri)
      for nc in nb:
        ns = cr.get(nc)
        if ns and len(ns) == 1: q.append(nc)
  psg(); ar = [i for i, d in enumerate(rows) if d]
  while len(ar) - len(cr) > 12: kr(max(ar, key=lambda i: len(rows[i]))); psg(); ar = [i for i, d in enumerate(rows) if d]
  while True:
    did = False
    for deg in (2, 3, 4):
      cd = None
      for c, s in list(cr.items()):
        if len(s) != deg: continue
        rl = list(s); mi = [MMAP.get(rows[r][c], -1) for r in rl]
        if any(k < 0 for k in mi): continue
        od = sorted(range(deg), key=lambda idx: len(rows[rl[idx]])); rp, kp = rl[od[0]], mi[od[0]]; lp = len(rows[rp])
        if deg >= 3 and any(lp + len(rows[rl[o]]) > 180 for o in od[1:]): continue
        cd = (rp, kp, [(rl[o], mi[o]) for o in od[1:]]); break
      if cd is not None:
        rp, kp, tgs = cd; dp = rows[rp]
        for rt, kt in tgs:
          dt, sh = rows[rt], (kt - kp) % 12
          for col, v0 in dp.items():
            sv = _zmulm(v0, sh)
            if col in dt:
              a = dt[col]; nv = (a[0] - sv[0], a[1] - sv[1], a[2] - sv[2], a[3] - sv[3])
              if nv == (0, 0, 0, 0):
                del dt[col]; cr[col].discard(rt)
                if not cr[col]: cr.pop(col, None)
              else: dt[col] = nv
            else: dt[col] = (-sv[0], -sv[1], -sv[2], -sv[3]); cr.setdefault(col, set()).add(rt)
          if not dt: kr(rt)
        kr(rp); psg(); did = True; break
    if not did: break
  al = [d for d in rows if d]; cols = sorted({c for d in al for c in d}); mx = ov = nz = nm = 0
  for d in al:
    for v in d.values():
      nz += 1
      if v in MMAP: nm += 1
      ma = max(abs(v[0]), abs(v[1]), abs(v[2]), abs(v[3]))
      if ma > mx: mx = ma
      if ma >= (1 << 30): ov += 1
  z1 = X_BLS % R_BLS; z2 = (z1 * z1) % R_BLS; z3 = (z2 * z1) % R_BLS; cm, piv, rk = {c: i for i, c in enumerate(cols)}, {}, 0
  for d in al:
    cur = {cm[c]: val for c, (c0, c1, c2, c3) in d.items() if (val := (c0 + c1 * z1 + c2 * z2 + c3 * z3) % R_BLS)}
    while cur:
      pc = min(cur)
      if pc in piv:
        fac = cur[pc]
        for c, pv in piv[pc].items():
          nv = (cur.get(c, 0) - fac * pv) % R_BLS
          if nv: cur[c] = nv
          else: cur.pop(c, None)
      else: iv = pow(cur[pc], -1, R_BLS); piv[pc] = {c: (v * iv) % R_BLS for c, v in cur.items()}; rk += 1; break
  return {"rows": len(al), "cols": len(cols), "mf": nm / max(1, nz), "mx": mx, "ov": ov, "null": min(len(al), len(cols)) - rk}

def check_live_sge(full: bool = False) -> dict:
  ns, fl, mxc, tov, mfs = 64, 0, 0, 0, 0.0
  for s in range(ns):
    r = _z12_sge_rank(512, 0.88, 7000 + s)
    if r["rows"] < r["cols"] or r["null"] != 0: fl += 1
    mxc = max(mxc, r["mx"]); tov += r["ov"]; mfs += r["mf"]
  mm, cp1, cp2 = mfs / ns, 1.0 - (0.05 ** (1.0 / ns)), 1.0 - (0.025 ** (1.0 / ns))
  assert fl == 0 and tov == 0 and mxc < (1 << 30) and mm > 0.75 and cp2 < REF["p_retry"]
  out = {"ok": True, "ring": "Z[zeta_12] -> BLS12-381 F_r", "seeds": ns, "fails": fl, "cp95_1s": round(cp1, 5), "cp95_2s": round(cp2, 5), "p_retry": REF["p_retry"], "max_abs": mxc, "ov30": tov, "mean_mono": round(mm, 4)}
  if full:
    ex = {}
    for eta in (0.85, 0.88, 0.90):
      ef = em = 0
      for s in range(4):
        r1 = _z12_sge_rank(1024, eta, 9100 + int(eta * 1000) * 10 + s)
        if r1["rows"] < r1["cols"] or r1["null"] != 0: ef += 1
        em = max(em, r1["mx"])
      assert ef == 0; ex[f"U1024_eta_{eta:.2f}"] = {"seeds": 4, "fails": ef, "max_abs": em}
    out["u1024"] = ex
  u24, la = REF["sge"][-1], REF["la"]; assert la["nl"] > u24[1] + 2.0 * u24[2] and la["sw"] > u24[3]
  out["la_hr_x"] = round((la["nl"] * la["sw"]) / (u24[1] * u24[3]), 3); return out

# 8. Unified Stage 4 Separate-(Bf, Bg) Descent, 12D LLL & Full Dusart-Exact Ledger
def _mul_r4(Bb: float) -> float:
  s, oe = REF["s1b"], REF["oe"]; L1, L2 = (oe * (2.0**60) * 8.0) / 4.5, (oe * (2.0**60) * 8.0) / 9.0
  if Bb <= 1024.0: nl = max(1.0, Bb / 64.0); return 12.0 * nl * nl * s["saf"]
  if Bb <= float(1 << 59):
    N = float(1 << int(math.ceil(math.log2(max(2.0, 2.0 * math.ceil(Bb / 64.0)))))); assert N <= (1 << 54)
    return (3.0 * (1.5 * N * math.log2(N) * s["bfly"] + N * (s["pw"] + s["sc"])) + ((3.0 * N * 160.0) if (Bb * 2.0 > L1) else 0.0) + N * s["c1"]) * s["saf"]
  L = int(math.ceil(Bb / 58.0)); k = max(1, (L - 1).bit_length()); k1 = (k + 1) // 2; k2 = k - k1; N2 = float(1 << (k1 + 1)) * float(1 << (k2 + 1))
  return (3.0 * (N2 * (1.0 * (k1 + 1) + 1.5 * (k2 + 1)) * s["bfly"] + N2 * (s["pw"] + s["sc"])) + ((9.0 * N2 * 160.0) if (Bb * 2.0 > L2) else 0.0) + N2 * s["c2"]) * s["saf"]

def _rem_r4(nlv: float, lb: float) -> float:
  if nlv <= 1.0: return 40000.0
  dm, kp = REF["s1b"]["dm"], max(1, int(4096.0 // lb)); ns, sb, nll = max(1.0, nlv / float(kp)), float(kp) * lb, lb / 64.0
  isb = sum(12.0 * (j * nll) * nll * 1.15 for j in range(1, kp)); isd = float(kp) * 24.0 * (sb / 64.0) * nll * 1.15
  h, ip, id_, cn, cb = int(math.ceil(math.log2(ns))), 0.0, 0.0, ns, sb
  for _ in range(h): pr = cn * 0.5; mir = _mul_r4(cb); ip += pr * mir; id_ += pr * dm * mir; cn, cb = pr, cb * 2.0
  return ns * (isb + isd) + ip + 3.0 * _mul_r4(ns * sb) + id_ + nlv * max(25000.0, 40000.0 * (lb / 636.0))

def _opt_ecm(qb: float, Nb: float, _ec: dict = {}, _ep: dict = {}):
  qk, key = round(qb, 2), (round(qb, 2), int(math.ceil(Nb / 64.0)))
  if key in _ec: return _ec[key]
  if qk not in _ep:
    pl, qe = [], max(1.0, qb - 3.5)
    for lb4 in range(32, 249):
      B1, l1 = 2.0 ** (lb4 * 0.25), lb4 * 0.25; p1 = DK.rho(qe / l1)
      for m in (10.0, 25.0, 60.0, 120.0, 250.0):
        B2 = B1 * m; dt = (math.log2(B2) - l1) / 120.0
        p2 = sum((DK.rho((qe - (l1 + (i + 0.5) * dt)) / l1) * (dt / (l1 + (i + 0.5) * dt))) if (l1 + (i + 0.5) * dt) < qe else (dt / (l1 + (i + 0.5) * dt)) for i in range(120))
        pr = min(1.0, p1 + p2)
        if pr > 1e-40: pl.append((B1, B2, pr))
    _ep[qk] = pl
  L = max(4, int(math.ceil(Nb / 64.0))); mir = 212.0 if L == 4 else (13.25 * L * L + 80.0 if L <= 16 else min(13.25 * L * L + 80.0, 2.5 * _mul_r4(L * 64.0))); best = (1e300, None)
  for B1, B2, pr in _ep[qk]:
    c1 = 1.15 * (40.0 * mir + (B1 / math.log(2.0)) * (11.0 * mir + 93.0 * L) + 2.5 * max(0.0, B2 / math.log(B2) - B1 / math.log(B1)) * (mir + 40.0 * L))
    if c1 / pr < best[0]: best = (c1 / pr, (B1, B2, pr, c1))
  _ec[key] = best; return best

def _dnode(qb: float, side: str, Bf: float, Bg: float, cf: float, cg: float):
  Bfl, Bgl, best = max(Bf, cf), max(Bg, cg), (1e300, 0.0, 0.0)
  for deg in (1, 2):
    for s2 in range(45):
      s = s2 * 0.5; nf = (4.5 if deg == 1 else 3.0) * qb + 72.0 * s + (111.0 if deg == 1 else 55.0)
      ng = (1.5 if deg == 1 else 1.0) * qb + 8.0 * (max(3.0 * s, LOG2_ABS_X - 3.0 * s) if deg == 1 else max(3.0 * s, LOG2_ABS_X, 2.0 * LOG2_ABS_X - 3.0 * s)) + (45.0 if deg == 1 else 35.0)
      nfr, ngr = (max(1.0, nf - qb), ng) if side == "f" else (nf, max(1.0, ng - qb)); l2T, dim = 20.0, 8.0 * (deg + 1)
      for _ in range(4):
        nft, ngt = nfr + (72.0 / dim) * l2T, ngr + (24.0 / dim) * l2T; pj = DK.rho(nft / Bfl) * DK.rho(ngt / Bgl)
        if pj <= 1e-250: l2T = 300.0; break
        l2T = -math.log2(pj) + 1.0
      if l2T <= 120.0:
        tu = ((2.0**l2T) * 2000.0 + 3.0 * (_opt_ecm(Bfl, nft)[0] + _opt_ecm(Bgl, ngt)[0]) + 1.0e7) / C_COST
        if tu < best[0]: best = (tu, nft, ngt)
  return best

def _s4_unified() -> dict:
  Bf, Bg, eta = REF["Bf"], REF["Bg"], REF["eta"]
  def br(Nr, Qh, Ql, M=300):
    dn, dx = DK.rho(Nr / Qh), (Qh - Ql) / M
    return sum((DK.rho((Nr - (Ql + (m + 0.5) * dx)) / Qh) / dn) * (dx / (Ql + (m + 0.5) * dx)) for m in range(M))
  Q = [275.0, 155.0, 130.0, 106.0, 86.0]; b0 = br(2572.0, Q[0], Q[1]); nd, bs = [b0], [round(b0, 3)]; ut1 = 5.0 * _dnode(550.0, "g", Bf, Bg, 275.0, 275.0)[0]
  for k in range(1, 5):
    cf, cg = Q[k], Q[k]; qfn, qgn = (Q[k + 1], Q[k + 1]) if k < 4 else (Bf, Bg)
    uf, nff, ngf = _dnode(Q[k - 1], "f", Bf, Bg, cf, cg); ug, nfg, ngg = _dnode(Q[k - 1], "g", Bf, Bg, cf, cg)
    bk = max(br(nff, cf, qfn) + br(ngf, cg, qgn), br(nfg, cf, qfn) + br(ngg, cg, qgn)); ut1 += nd[-1] * max(uf, ug); bs.append(round(bk, 3)); nd.append(nd[-1] * bk)
  ut1 += nd[-1] * max(_dnode(Q[4], "f", Bf, Bg, Bf, Bg)[0], _dnode(Q[4], "g", Bf, Bg, Bf, Bg)[0])
  ub1 = ((1.0 / (0.918 * DK.rho(2572.0 / 275.0))) * (2.5 * _opt_ecm(275.0, 2572.0)[0])) / C_COST; mf = (1.04094 / eta) ** 2
  return {"bs": bs, "nodes": [round(x, 1) for x in nd] + [0.0], "s4": math.log2(2.0 * (ub1 + ut1) * mf), "ecm": math.log2(2.0 * (ut1 + 0.25 * ub1) * mf)}

def _lll(B: list[list[int]], delta: float = 0.75) -> list[list[int]]:
  n, B = len(B), [r[:] for r in B]
  def gso():
    bo, mu, B2 = [[0.0] * len(B[0]) for _ in range(n)], [[0.0] * n for _ in range(n)], [0.0] * n
    for i in range(n):
      v = [float(x) for x in B[i]]
      for j in range(i):
        mu[i][j] = sum(B[i][t] * bo[j][t] for t in range(len(v))) / B2[j] if B2[j] > 1e-18 else 0.0
        for t in range(len(v)): v[t] -= mu[i][j] * bo[j][t]
      bo[i] = v; B2[i] = sum(x * x for x in v)
    return mu, B2
  mu, B2 = gso(); k = 1
  while k < n:
    for j in range(k - 1, -1, -1):
      if abs(mu[k][j]) > 0.5:
        q = int(round(mu[k][j]))
        for t in range(len(B[0])): B[k][t] -= q * B[j][t]
        mu, B2 = gso()
    if B2[k] >= (delta - mu[k][k - 1] ** 2) * B2[k - 1]: k += 1
    else: B[k], B[k - 1] = B[k - 1], B[k]; mu, B2 = gso(); k = max(k - 1, 1)
  return B

def check_stage4(full: bool = False) -> dict:
  s4 = _s4_unified(); assert s4["s4"] < REF["cr"]["both_individual_logs_exp"] and s4["ecm"] < REF["cr"]["ecm_split_expected"]
  pt, rl = 1009, random.Random(12012); T = [[rl.randrange(pt) for _ in range(6)] for _ in range(6)]
  B12 = [[(1 if i == j else 0) for j in range(6)] + T[i] for i in range(6)] + [[0] * 6 + [(pt if i == j else 0) for j in range(6)] for i in range(6)]
  sl2 = math.sqrt(sum(x * x for x in _lll(B12)[0])); assert sl2 < 2.0 * math.sqrt(12.0 * pt)
  lll = 33.0 - 6.0 * math.log2(12.0) + 4.5 * math.log2(5.0) + 6.0 * math.log2(3.0) + 6.0 * math.log2(P_BLS) + 4.0 * math.log2(ABS_X); assert lll < 2572.0
  out = {"ok": True, "lll_ub": round(lll, 4), "term": {"Bf": REF["Bf"], "Bg": REF["Bg"], "out_f": 0, "out_g": 0}, "bs": s4["bs"], "nodes": s4["nodes"], "s4_l2": round(s4["s4"], 4), "ecm_l2": round(s4["ecm"], 4), "toy_l2": round(sl2, 2)}
  if not full: return out
  xd = lambda X, Z, A24, N: (((X + Z) ** 2 * (X - Z) ** 2) % N, ((((X + Z) ** 2 - (X - Z) ** 2) % N) * (((X - Z) ** 2 + A24 * ((X + Z) ** 2 - (X - Z) ** 2)) % N)) % N)
  xa = lambda X1, Z1, X2, Z2, Xd, Zd, N: ((Zd * (((X1 - Z1) * (X2 + Z2) + (X1 + Z1) * (X2 - Z2)) ** 2)) % N, (Xd * (((X1 - Z1) * (X2 + Z2) - (X1 + Z1) * (X2 - Z2)) ** 2)) % N)
  def xm(k, X, Z, A24, N):
    r0x, r0z, r1x, r1z = 1, 0, X % N, Z % N
    for bt in bin(k)[2:]:
      if bt == "0": r1x, r1z = xa(r0x, r0z, r1x, r1z, X, Z, N); r0x, r0z = xd(r0x, r0z, A24, N)
      else: r0x, r0z = xa(r0x, r0z, r1x, r1z, X, Z, N); r1x, r1z = xd(r1x, r1z, A24, N)
    return r0x, r0z
  nx, dx = 5, -7; t16 = [p for p in _primes(200) if p >= 101]
  for p in t16:
    x = (nx * pow(dx, -1, p)) % p; d = ((2 * x * x - 1) * pow(x, -4, p)) % p; A = (2 * (1 + d) * pow(1 - d, -1, p)) % p
    assert (1 + sum(1 if (u**3 + A * u * u + u) % p == 0 else (2 if pow((u**3 + A * u * u + u) % p, (p - 1) // 2, p) == 1 else 0) for u in range(p))) % 16 == 0
  N1 = 50023 * 50033; x1 = (nx * pow(dx, -1, N1)) % N1; A24 = pow(1 - ((2 * x1 * x1 - 1) * pow(x1, -4, N1)) % N1, -1, N1); X, Z = 3, 1
  for p in _primes(30):
    pk = p
    while pk * p <= 30: pk *= p
    X, Z = xm(pk, X, Z, A24, N1)
  assert math.gcd(Z, N1) == 1; bb = [xm(i, X, Z, A24, N1) for i in range(1, 16) if math.gcd(i, 30) == 1]; DX, DZ = xm(30, X, Z, A24, N1); gp, gc, w1 = (1, 0), (DX, DZ), 1
  for k in range(1, 13):
    for bx, bz in bb:
      g = math.gcd((gc[0] * bz - bx * gc[1]) % N1, N1)
      if 1 < g < N1: w1 = g; break
    if w1 > 1: break
    gp, gc = gc, (xd(DX, DZ, A24, N1) if k == 1 else xa(gc[0], gc[1], DX, DZ, gp[0], gp[1], N1))
  X2, Z2 = xm(2, 4, 1, 2, 35); DX2, DZ2 = xm(4, X2, Z2, 2, 35); w2 = math.gcd((DX2 * Z2 - X2 * DZ2) % 35, 35)
  assert w1 == 50023 and math.gcd(Z2, 35) == 1 and w2 == 5; out.update({"t16": f"{len(t16)}/{len(t16)}", "w1": w1, "w2": [w2, 35 // w2]}); return out

def _tx64(lam: float) -> float:
  if lam <= 0.0: return 0.0
  if lam > 45.0: return (69.0 / 64.0) * lam + 24.0
  ei, pk = 0.0, math.exp(-lam)
  for k in range(1, 160):
    pk *= lam / k; q, r = divmod(k, 64)
    ei += pk * (10.0 + q * 69.0 + (35.0 if r & 32 else 0) + (19.0 if r & 16 else 0) + (11.0 if r & 8 else 0) + (7.0 if r & 4 else 0) + (3.0 if r & 2 else 0) + (3.0 if r & 1 else 0))
    if k > lam + 25.0 and pk < 1e-14: break
  return ei

def _w4r(lam: float) -> float:
  if lam <= 0.0: return 0.0
  if lam > 40.0: return 13.5 + (13.0 / 8.0) * lam + 9.5
  ei, pk = 0.0, math.exp(-lam)
  for k in range(1, 140):
    pk *= lam / k; ik = 6.5 if k == 1 else (9.5 if k == 2 else 13.5 + (k // 8) * 13.0 + (7.0 if k & 4 else 0) + (3.75 if k & 2 else 0) + (3.0 if k & 1 else 0)); ei += pk * ik
    if k > lam + 25.0 and pk < 1e-14: break
  return ei

def _dm16(lb: float) -> float:
  if lb <= 0.0: return 0.0
  if lb > 35.0: return 4.0 + 6.50 * lb + 2.25
  ei, pk = 0.0, math.exp(-lb)
  for k in range(1, 120):
    pk *= lb / k; q, r = divmod(k, 4); ei += pk * (4.0 + 26.0 * q + (1.0 if q > 0 else 0.0) + 8.0 * r)
    if k > lb + 25.0 and pk < 1e-14: break
  return ei

def _integ(plo, pw1, pw4, phi, Lw, dok, rmx, wb3=16):
  wpc = 1.0 / (V_A_SH * Lw); ih = ip = nok = nb = n1 = n2 = n3 = nh = 0.0
  if pw1 > plo:
    dx = (pw1 - plo) / 100.0
    for i in range(100):
      x = plo + (i + 0.5) * dx; p = 2.0**x; d = (p / x) * dx; dk = d * dok; nok += dk; nb += dk * (p + 2.0); n1 += d * V_A_SH
      ip += dk * (44.0 + Lw * 7.0) * wpc; ih += dk * sum(pr * Lw * (_tx64(r * V_A_SH / p) + 1.5625 * (r * V_A_SH / p)) for r, pr in rmx) * wpc
  if pw4 > pw1:
    nb4, dx = Lw / 4.0 + 1.0, (pw4 - pw1) / 60.0
    for i in range(60):
      x = pw1 + (i + 0.5) * dx; p = 2.0**x; d = (p / x) * dx; dk = d * dok; nbp = dk * (p / 4.0 + 2.0); nok += dk; nb += nbp; n2 += d * V_A_SH; eip = eih = 0.0
      for r, pr in rmx:
        lm = r * V_A_SH / p; pb4 = 1.0 - math.exp(-4.0 * lm); pc = 1.0 - math.exp(-Lw * lm); nh += pr * nbp * pb4
        eip += pr * (44.0 + pc * nb4 * 5.0 + nb4 * pb4 * 12.0); eih += pr * Lw * (_w4r(lm) + 1.5625 * lm)
      ip += dk * eip * wpc; ih += dk * eih * wpc
  if phi > pw4:
    nb3, dx = Lw / float(wb3) + 1.0, (phi - pw4) / 60.0
    for i in range(60):
      x = pw4 + (i + 0.5) * dx; p = 2.0**x; d = (p / x) * dx; dk = d * dok; nok += dk; nb += dk * (p / float(wb3) + 2.0); n3 += d * V_A_SH; eip = eih = 0.0
      for r, pr in rmx:
        lb = float(wb3) * r * V_A_SH / p; pc = 1.0 - math.exp(-nb3 * lb); la = r * V_A_SH * Lw / p
        eip += pr * (44.0 + pc * nb3 * 5.0); eih += pr * (nb3 * _dm16(lb) + la * 1.5625)
      ip += dk * eip * wpc; ih += dk * eih * wpc
  return {"ih": ih, "ip": ip, "nok": nok, "nb": nb, "n1": n1, "n2": n2, "n3": n3, "nh": nh}

def _tprob(cof, Pt, Bb, km=8):
  if cof <= Pt: return 1.0, 0.0, 1.0
  om = lambda v: 0.0 if v < 1.0 else (1.0 / v if v <= 2.0 else ((1.0 + math.log(v - 1.0)) / v if v <= 3.0 else 0.57))
  dy, md, ma = (cof - Pt) / 200.0, DK.rho(cof / Pt), 0.0; mt = md
  for i in range(200):
    y = Pt + (i + 0.5) * dy; w = DK.rho((cof - y) / Pt) * (om(y / Pt) / Pt) * dy; mt += w
    if y <= Bb: md += w
    elif any(k * Pt <= y <= k * Bb for k in range(2, km + 1)): ma += w
  return md / max(mt, 1e-30), ma / max(mt, 1e-30), (md + ma) / max(mt, 1e-30)

def check_ledger() -> dict:
  Bf, Bg, eta, pr = REF["Bf"], REF["Bg"], REF["eta"], REF["p_retry"]; b, mc, pol, rep, la = REF["b"], REF["mc"], REF["pol"], REF["rep"], REF["la"]
  dokg, dokf, dng, dnf = 1.0 / 3.0, 23.0 / 81.0, 1.0 / 24.0, 2131.0 / 34992.0; w0 = 1.0 / (1.0 + REF["rbr"] ** 8); w1 = 1.0 - w0
  fmx = [(3, 39.0 / 46.0), (6, 6.0 / 46.0), (9, 1.0 / 46.0)]
  sg0, sf0 = _integ(b["pglo"], b["pw1"], b["pw4"], b["pghi0"], mc["L4"], dokg, [(3, 1.0)]), _integ(b["pflo"], b["pw1"], b["pw4"], b["pfhi0"], mc["L4"], dokf, fmx)
  sg1, sf1 = _integ(b["pglo"], b["pw1"], b["pw4"], b["pghi1"], mc["L4"], dokg, [(3, 1.0)]), _integ(b["pflo"], b["pw1"], b["pw4"], b["pfhi1"], mc["L4"], dokf, fmx)
  ih = w0 * (sg0["ih"] + sf0["ih"]) + w1 * (sg1["ih"] + sf1["ih"]); ip = w0 * (sg0["ip"] + sf0["ip"]) + w1 * (sg1["ip"] + sf1["ip"])
  def rsh(sg, sf):
    n12, n3, nb, nh = sg["n1"] + sf["n1"] + sg["n2"] + sf["n2"], sg["n3"] + sf["n3"], sg["nb"] + sf["nb"], sg["nh"] + sf["nh"]
    return (n12 * rep["t12"] + n3 * rep["t3"] + nb * rep["bkt"] + nh * rep["hdr"]) / C_COST, math.log2(8.0 * (n12 + nh + nb + n3) + 48.0 * (sg["nok"] + sf["nok"]))
  rep0, mem0 = rsh(sg0, sf0); rep1, mem1 = rsh(sg1, sf1); repu = rep0 + rep1
  l2Vk = mc["Vk4"] + mc["cpad"] + mc["sh"] + 3.0 * mc["se"] + math.log2(1.0 / (1.0 - (mc["g4"] + 0.00010)))
  assert abs(l2Vk - mc["Vkch"]) < 1e-9; Vt, LtM, lf, lg = 2.0**l2Vk, 0.0, REF["lf"], REF["lg"]
  def ps1b(pfm, pgm):
    nonlocal LtM
    s1 = 2.0 ** pol["s1"]; c1, c2, mg = pfm + (Bf - pfm) / 3.0, pfm + 2.0 * (Bf - pfm) / 3.0, 0.5 * (pgm + Bg)
    _, af1, _ = _tprob(lf, pfm, Bf); _, af2, _ = _tprob(lf, c1, Bf); _, af3, _ = _tprob(lf, c2, Bf); _, ag1, _ = _tprob(lg, pgm, Bg); _, ag2, _ = _tprob(lg, mg, Bg)
    rf = min(0.95, math.sqrt(pol["ret"] / 0.80)); rg = pol["ret"] / rf
    pfs = min(1.0, (DK.rho((mc["mnf"] + REF["af"]) / Bf) * rf) / (2.0 ** pol["sf"])); pgs = min(1.0, (DK.rho((mc["mng"] + REF["ag"]) / Bg) * rg) / (2.0 ** pol["sg"]))
    cr = s1 * 1700.0 + s1 * pfs * 1700.0
    for lo, hi, pin, lb, dn in [(1.0, pfm, s1, 656.0, dnf), (pfm, c1, s1 * af1, lf, dnf), (c1, c2, s1 * af2, lf, dnf), (c2, Bf, s1 * af3, lf, dnf), (1.0, pgm, s1 * pfs, 688.0, dng), (pgm, mg, s1 * pfs * ag1, lg, dng), (mg, Bg, s1 * pfs * ag2, lg, dng)]:
      Lt = max(1.0, (ppb_dus(lo, min(hi, 36.0)) if lo < 36.0 else 0.0) + (dn * ppb_dus(max(lo, 36.0), hi) if hi > 36.0 else 0.0))
      LtM = max(LtM, Lt); Mt = Vt * pin; Mb = min(Mt, max(1.0, Lt / lb)); nbb, LX = max(1.0, Mt / Mb), Mb * lb; ir = _rem_r4(Mb, lb)
      if LX < Lt: ir += max(0.0, Lt / LX - 1.0) * 2.5 * _mul_r4(LX)
      cr += (nbb * ir) / Vt
    return cr + s1 * pfs * pgs * (110.0 * 256.0 * 150.0 * 4.0 + 0.05 * 5.0e7 * 2.0)
  cs1b = w0 * ps1b(b["pfhi0"], b["pghi0"]) + w1 * ps1b(b["pfhi1"], b["pghi1"])
  sll = 0.171875 + 0.234375 + 150.0 * (2.0 ** pol["sg"]); slc = sll + 0.002; ipcl, ipcc = ih + ip + sll + cs1b, ih + ip + slc + cs1b
  c_l, c_c = (Vt * ipcl) / C_COST, (Vt * ipcc) / C_COST; U = U_ub(Bf, Bg); Nl, Nz = U * la["nl"], U * la["sw"]; u24 = REF["sge"][-1]
  assert la["nl"] > u24[1] + 2.0 * u24[2] and la["sw"] > u24[3] and (la["nl"] * la["sw"]) > 2.15 * (u24[1] * u24[3])
  bsp = 0.8143 * 10.0 + 0.1857 * 25.0 + (Nl * 8.0) / Nz; assert bsp < la["x4"] and (56.0 + 872.0) < la["rv"]
  Mb, Nb = la["m"], la["n"]; Lsq, dmk = Nl / Mb + Nl / Nb + 64.0, Nl / Nb + 32.0
  ibw = 1.02 * (Nb * Lsq * (Nz * la["x4"] + Nl * la["rv"]) * la["ov"] + dmk * (Nz * la["x1"] + Nl * la["rv"]) * la["ov"])
  imk = 1.02 * dmk * Nb * Nl * 283.0 * la["ov"]; ilg = 20.0 * ((Mb + Nb) ** 3) * Lsq * math.log2(Lsq) * 236.0 * la["ov"]
  lau = (ibw + imk + ilg + 0.008 * (ibw + imk) + 30.0 * U * 283.0 * la["ov"]) / C_COST; sgeu = ((12500.0 * U + 25.0 * U * math.log2(U) + 100.0 * Nz) * 1.15) / C_COST
  snl, snc = pr * (((1.0 - eta) / eta) * c_l + repu + lau), pr * (((1.0 - eta) / eta) * c_c + repu + lau); s4 = _s4_unified()
  lr = {"relation_collection": math.log2(c_l), "linear_algebra": math.log2(lau), "replay_a_tables": math.log2(repu), "snell_resume_tail": math.log2(snl), "both_individual_logs_exp": s4["s4"], "ecm_split_expected": s4["ecm"], "meter_safety_buffer": 93.4, "band_split_and_trees": 90.0, "canonical_a_domain": 84.0, "schirokauer_maps": 80.0, "final_log_backsub_and_gcd": 72.0, "matrix_generator_sge": math.log2(sgeu), "stage0_setup_and_iso": 45.0, "final_verification_Q": 40.0}
  cr = REF["cr"]; assert len(lr) == len(cr) == 14 and all(cr[k] > lr[k] for k in lr) and cr["relation_collection"] >= math.log2(c_c) and cr["snell_resume_tail"] >= math.log2(snc)
  wol, woc = sum(2.0**v for v in lr.values()), sum(2.0**v for v in cr.values()); lol, loc = math.log2(wol), math.log2(woc)
  pb = math.ceil((2.0**mem0) / (12.0 * (2.0**60))); gsB = 248.0 * V_A_SH * (2.0 ** pol["sg"]) * 8.0
  s1p = math.log2((2.0**mem0) + 248.0 * (V_A_SH + 262144.0) + gsB); ntB, cb, nn, lv = 0.0, 4080.0, LtM / 4080.0, 0
  while nn >= 1.0:
    if cb <= (1 << 59): ntB += nn * (float(1 << int(math.ceil(math.log2(max(2.0, 2.0 * math.ceil(cb / 64.0)))))) * 24.0)
    else: L2 = int(math.ceil(cb / 58.0)); k = max(1, (L2 - 1).bit_length()); k1 = (k + 1) // 2; ntB += nn * (float(1 << (k1 + 1)) * float(1 << (k - k1 + 1)) * 24.0)
    cb *= 2.0; nn *= 0.5; lv += 1
  s1bB = ntB + (9.0 + lv) * LtM / 8.0; s1b_oe = s1bB / (2048.0 * (2.0**60)); s2r = math.log2((eta * U) * 18.0 * 12.0)
  lmB, lvB = Nz * (0.8143 * 8.0 + 0.1857 * 12.0) + Nl * 8.0, Nl * 32.0; la_oe = (lmB / la["p"] + lvB) / (2.0**60)
  s3r = math.log2(lmB + (Nb + 6.0 + la["p"]) * lvB + (Mb * Nb * Lsq + 4.0 * Nb * Nb * dmk + 2.0 * Mb * (Mb + Nb) * Lsq) * 32.0)
  assert V_A_SH / (2.0**60) < 0.06 and s1b_oe < REF["oe"] < 16.0 and la_oe < 4.0 < 16.0 and max(s1p, math.log2(s1bB), s2r, s3r) < REF["ram_ub"] < 120.0
  mu, ds = eta * U, (2.0 ** (3.0 * mc["se"])) - 1.0; chl, chh = 0.5 * mu * ds * ds, (mu * ds * ds) / (2.0 + ds); assert chl > 1.9e13 and chh > 1.9e13
  pm = woc / (2.0**124.0); pf = pm + 2.0**(-48) + 2.0**(-100) + 2.0**(-132) + 2.0**(-200); po = 1.0 - pf; Kg, t64 = REF["Kg"], 1 << 64
  assert pf < 2.0e-5 and po > 0.99998 and Fraction(Kg, t64) * Fraction(99998, 100000) >= Fraction(39, 100) and Fraction(Kg - 1, t64) * Fraction(99998, 100000) < Fraction(39, 100)
  pc = (Kg / float(t64)) * po; num, den = REF["num"], REF["den"]; lgc, lgl = loc + math.log2(Kg / float(t64)), lol + math.log2(Kg / float(t64))
  assert pc > 0.39 and math.gcd(num, den) == 1 and lgc < float(Fraction(num, den))
  return {
    "ok": True, "claim": f"{num}/{den}", "U_l2": round(math.log2(U), 6), "U_ub": float(f"{U:.8e}"),
    "c_s1b": round(cs1b, 6), "slab_live": round(sll, 6), "slab_ch": round(slc, 6), "ipc_live": round(ipcl, 6), "ipc_ch": round(ipcc, 6),
    "blended_spmv": round(bsp, 4), "la_mat_eib": round(lmB / (2.0**60), 4), "la_owner_eib": round(la_oe, 4),
    "opened_live": round(lol, 6), "opened_ch": round(loc, 6), "gated_live": round(lgl, 6), "gated_ch": round(lgc, 6), "slack_b": round(float(Fraction(num, den)) - lgc, 6),
    "chernoff_H1": [float(f"{chl:.4e}"), float(f"{chh:.4e}")], "p_meter": float(f"{pm:.7e}"), "p_open": round(po, 8), "p_cold": round(pc, 8),
    "ledger_live": {k: round(v, 6) for k, v in lr.items()},
    "ram_l2_B": {"S1_sh0": round(s1p, 4), "S1_sh1": round(mem1, 4), "P_banks": pb, "S1b": round(math.log2(s1bB), 4), "S1b_eib": round(s1b_oe, 4), "S2": round(s2r, 4), "S3": round(s3r, 4)},
  }

def main() -> None:
  t0 = time.perf_counter(); raw = sys.stdin.read().strip(); inp = json.loads(raw) if raw else {}
  if not isinstance(inp, dict) or (set(inp.keys()) - {"part"}): raise ValueError("Input must be a JSON object with optional key 'part'")
  part = inp.get("part", "core")
  if part not in ("core", "diagnostics", "all"): raise ValueError(f"Unsupported part: {part}")
  out: dict = {"package_id": REF["pkg"], "config_tag": REF["tag"], "part": part, "reads_external_files": 0}
  if part in ("core", "all"):
    out["stage0"] = check_stage0()
    out["chebotarev"] = check_chebotarev(26000)
    out["murphy_norms"] = check_murphy_norms()
    out["live_panjer_mc_chords"] = check_live_panjer_mc_chords(10000 if part == "all" else 5000, 10000 if part == "all" else 5000, 2000 if part == "all" else 1000)
    out["ntt_r4_shoup61"] = check_ntt()
    out["kernels_e2e"] = check_kernels_e2e()
    out["live_strict_z12_sge"] = check_live_sge(full=(part == "all"))
    out["stage4_unified"] = check_stage4(full=(part == "all"))
    out["ledger_ram_budget"] = check_ledger()
  elif part == "diagnostics":
    out["live_panjer_mc_chords"] = check_live_panjer_mc_chords(10000, 10000, 2000)
    out["live_strict_z12_sge"] = check_live_sge(full=True)
    out["stage4_unified"] = check_stage4(full=True)
  out["checks_ok"] = all(v.get("ok", True) for v in out.values() if isinstance(v, dict))
  out["elapsed_seconds"] = round(time.perf_counter() - t0, 3)
  print(json.dumps(out, indent=2))
  if not out["checks_ok"]: sys.exit(1)

if __name__ == "__main__":
  main()
