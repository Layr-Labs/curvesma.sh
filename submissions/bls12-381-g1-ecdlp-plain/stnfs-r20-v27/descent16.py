"""stnfs descent model (DESCENT.md):"""
import math, functools
from decimal import Decimal, localcontext
LN2 = math.log(2); UNIT = 20000.0
KMAX = 64
RSTEP = 1.0 / 1024
def _rho_coefficients(kmax=KMAX, nt=600):
 with localcontext() as ctx:
  ctx.prec = 200
  zero, one = Decimal(0), Decimal(1)
  coefs = [None, [one] + [zero] * (nt - 1)]
  for k in range(2, kmax + 1):
   previous, c = coefs[k - 1], [zero] * nt
   for m in range(nt - 1):
    c[m + 1] = (previous[m] + m * c[m]) / (k * (m + 1))
   c[0] = previous[0] - sum(c[1:])
   coefs.append(c)
  out = [None]
  for c in coefs[1:]:
   big = max(abs(x) for x in c)
   while len(c) > 1 and abs(c[-1]) < big * Decimal('1e-25'):
    c = c[:-1]
   out.append([float(x) for x in c])
  return out
RHO_COEFS = _rho_coefficients()
def rho_series(u):
 if u <= 1:
  return 1.0
 k = math.ceil(u)
 if k > KMAX:
  return 0.0
 x, value = k - u, 0.0
 for c in reversed(RHO_COEFS[k]):
  value = value * x + c
 return max(value, 0.0)
RHO_TAB = [rho_series(i * RSTEP) for i in range(int(KMAX / RSTEP) + 2)]
def rho_lo(u):
 """lower bracket of rho(u): the table value at the next grid point."""
 if u <= 1:
  return 1.0
 i = math.ceil(u / RSTEP)
 return RHO_TAB[i] if i < len(RHO_TAB) else 0.0
def rho_hi(u):
 """upper bracket: the previous grid point;"""
 if u <= 1:
  return 1.0
 i = math.floor(u / RSTEP)
 if i + 1 >= len(RHO_TAB):
  return 1.0 / math.factorial(math.floor(u))
 return RHO_TAB[i]
def kappa(u, t):
 return rho_series(u - t) / rho_series(u) if rho_series(u) > 0 else 0.0
def pois2(lam):
 return math.exp(-lam) * (1 + lam + lam * lam / 2)
# H5c purge-loss domination (red line n; DESCENT.md).
H5C_MUC=1.2;H5C_C0=1.01;H5C_TRI=2.0;H5C_BBM=5.0;H5C_RHO=1e-2;H5C_PHI=0.3;H5C_MUM=0.99
def _e(x): return math.exp(-x)
def h5c_r1(lam,d):
 if lam<=0 or d<=0: return 0.0
 x=lam*d
 if x>50: return _e(lam*(1-d)) if d<1 else 1.0
 return _e(lam)*(math.exp(x)-1-x-x*x/2)
def h5c_vi(lam,d):
 if lam<=0 or d<=0: return 0.0
 x=lam*d
 if x>50: return lam*_e(lam*(1-d)) if d<1 else 1.0
 return lam*_e(lam)*(math.exp(x)-1-x)
def h5c_o3(lam,ds):
 return max(_e(lam*(1-ds))-_e(lam)*(1+lam*ds),0.0) if lam>0 else 0.0
def h5c_extra(lam,dpp,bbar,dall,bun):
 e=h5c_r1(lam,dpp)+bbar*h5c_vi(lam,dpp)+H5C_RHO*h5c_r1(lam,dall)+h5c_r1(lam,H5C_MUM*bun)
 return min(max(e,0.0),max(1-pois2(lam),0.0))
def _gup(f,lo,hi,c2):
 return max(f(lo+(hi-lo)*i/8) for i in range(9))+(c2*(hi-lo)**2/512 if hi>lo else 0.0)
BANDS = [40 + i / 50.0 for i in range(1450)]   # bands [2^40, 2^69)
def h5c_tab(lr,um):
 out=[]
 for b in BANDS:
  a=kappa(um,b/69.0);c=kappa(um,(b+0.02)/69.0)
  klo,khi=min(a,c),max(a,c)
  lo=2.0**(lr-b-0.02)*2*klo;hi=2.0**(lr-b)*2*khi
  if lo>hi: lo,hi=hi,lo
  out.append((2*khi*math.log((b+0.02)/b),lo,hi))
 return out
def h5c_masses(tab,phi):
 dl=be=o2=0.0
 for w,lo,hi in tab:
  dl+=w*_e(lo)
  be+=w*(hi*_e(hi) if hi<=1 else (lo*_e(lo) if lo>=1 else 0.36787944117144233))
  o2+=w*(hi*hi*_e(hi)/2 if hi<=2 else (lo*lo*_e(lo)/2 if lo>=2 else 0.2706705664732254))
 def o3s(ds):
  return sum(w*_gup(lambda l:h5c_o3(l,ds),lo,hi,_e(lo*(1-ds))+_e(lo)*(1+lo*ds)) for w,lo,hi in tab)
 mu0=H5C_MUC*(be*(1-dl)+o3s(dl))
 muc=max(mu0,H5C_MUC*(be*0.68+o3s(0.32)))
 if muc>=1: raise ValueError('h5c!qb')
 ds=dl/(1-muc)*(1+1e-9)+1e-12
 if ds>0.32: raise ValueError('h5c!ds')
 o3=o3s(ds)
 casc=ds-dl;dpp=dl+H5C_C0*casc;dall=dl+casc+phi
 sw=sum(w for w,lo,hi in tab)
 r1m=sum(w*_gup(lambda l:h5c_vi(l,dpp)/l*dpp,lo,hi,dpp*(1+dpp)**2*_e(lo*(1-dpp))) for w,lo,hi in tab)
 vim=sum(w*_gup(lambda l:((l*dpp)+1)*_e(l*(1-dpp))-(2*(l*dpp)+1)*_e(l),lo,hi,_e(lo*(1-dpp))*(((hi*dpp)+1)+2*dpp*((hi*dpp)+2)+dpp*dpp*((hi*dpp)+3))) for w,lo,hi in tab)
 clm=H5C_RHO*(sw if dall>=1 else sum(w*_gup(lambda l:h5c_vi(l,dall)/l*dall,lo,hi,dall*(1+dall)**2*_e(lo*(1-dall))) for w,lo,hi in tab))
 mgm=H5C_MUM*r1m/dpp
 pr=dl*dl/2;tri=H5C_TRI*12*dl**3*o2**2;b=pr+tri+r1m+clm
 for _ in range(100):
  b2=pr+tri+o2*b*b+r1m+b*vim+clm+mgm*b
  if abs(b2-b)<1e-15: b=b2;break
  b=b2
 b=b*(1+1e-9)+1e-12
 if pr+tri+o2*b*b+r1m+b*vim+clm+mgm*b>b: raise ValueError('h5c!bbar')
 if not (dpp<1 and H5C_MUM*b<1 and H5C_MUM*b<=dpp): raise ValueError('h5c!prem')
 return dict(delta=dl,beta=be,o2=o2,ds=ds,o3=o3,muc=muc,casc=casc,dpp=dpp,r1mass=r1m,vimass=vim,bbar=min(b*H5C_BBM,1.0),bun=b,dall=dall)
_H5CM={};_H5CQ={}
def h5c_m(lr,uf,ug):
 k=(lr,uf,ug)
 if k not in _H5CM: _H5CM[k]=h5c_masses(h5c_tab(lr,uf)+h5c_tab(lr,ug),H5C_PHI)
 return _H5CM[k]
def h5c_mp(lr,uf,ug,um):
 k=(lr,uf,ug,um)
 if k not in _H5CQ:
  M=h5c_m(lr,uf,ug)
  dpp,bbar,dall,bun=M['dpp'],M['bbar'],M['dall'],M['bun']
  if dall>=1: raise ValueError('h5c!dall')
  mp={}
  for b,(w,lo,hi) in zip(BANDS,h5c_tab(lr,um)):
   c2=(2*dpp+1)*_e(lo*(1-dpp))+bbar*_e(lo*(1-dpp))*((2+hi)+2*dpp*(1+hi)+hi*dpp*dpp)+H5C_RHO*(2*dall+1)*_e(lo*(1-dall))+(2*H5C_MUM*bun+1)*_e(lo*(1-H5C_MUM*bun))
   mp[b]=min(1.0,pois2(lo)+_gup(lambda l:h5c_extra(l,dpp,bbar,dall,bun),lo,hi,c2))
  _H5CQ[k]=(M,mp)
 return _H5CQ[k]
def band(log2_rels, u_main, b):
 """band [2^b, 2^(b+.02)): m_I^+ upper (DESCENT.md H5c)"""
 if log2_rels!=P.LOG2_RELS or u_main not in (P.MAIN_UF,P.MAIN_UG): raise ValueError('h5c!params')
 return h5c_mp(log2_rels,P.MAIN_UF,P.MAIN_UG,u_main)[1][b], math.log((b + 0.02) / b)
@functools.lru_cache(None)
def live(log2_rels, u_main):
 """bands with nonzero m (exact zeros contribute nothing): (b, m, d)"""
 return [(b,) + band(log2_rels, u_main, b) for b in BANDS if band(log2_rels, u_main, b)[0] > 0]
@functools.lru_cache(None)
def miss_side(u4, log2_rels, u_main):
 """mu_hi(u) >= E S_u(M) = sum_I w_I m_I, w_I = kappa(u, t_I)/N(I)"""
 return sum(kappa(u4 / 4.0, (b + 0.02) / 69.0) * d * m for b, m, d in live(log2_rels, u_main))
def missing(Pm, nf, ng):
 """missing-log ideals per relation, charged at (1 + MISS_EPS) mu_hi"""
 f = lambda m, um: miss_side(math.ceil(4 * m / 69.0), Pm.LOG2_RELS, um)
 return (1 + Pm.MISS_EPS) * ((f(nf, Pm.MAIN_UF) if nf else 0.0) + (f(ng, Pm.MAIN_UG) if ng else 0.0))
def miss_tail(Pm, side, m):
 """log2 of Chebyshev's P(S_u(M) > (1 + MISS_EPS) mu_hi(u)) <= Var /..."""
 u4 = math.ceil(4 * m / 69.0); um = Pm.MAIN_UF if side == 'f' else Pm.MAIN_UG
 return math.log2(miss_pair(u4, u4, Pm.LOG2_RELS, um)) - 2 * math.log2(Pm.MISS_EPS * miss_side(u4, Pm.LOG2_RELS, um))
def p_none(Pm, nf, ng):
 """P(no missing ideal | M), Bonferroni lower bound per side (DESCENT.md)"""
 out = 1.0
 for m, um in ((nf, Pm.MAIN_UF), (ng, Pm.MAIN_UG)):
  if m:
   u4 = math.ceil(4 * m / 69.0); u = u4 / 4.0; S = (1 + Pm.MISS_EPS) * miss_side(u4, Pm.LOG2_RELS, um)
   r = kappa(u, 2.0) / kappa(u, 1.0) ** 2; c = kappa(u, (Pm.MISS_LO + 0.02) / 69.0) * 2.0 ** -Pm.MISS_LO
   out *= max(0.0, 1 - S + r * (S * S - c * S) / 2 - S ** 3 / 6)
 return out
@functools.lru_cache(None)
def miss_pair(u4i, u4j, log2_rels, u_main):
 """sum over sigma-orbits C of w_C(i) w_C(j) m_C"""
 return 2 * sum(kappa(u4i / 4.0, (b + 0.02) / 69.0) * kappa(u4j / 4.0, (b + 0.02) / 69.0) * 2.0 ** -b * d * m
            for b, m, d in live(log2_rels, u_main))
KA_MAX = 16
def _conv(a, b):
 """out[k] = sum_i a[i] b[k - i] over the nonzero supports"""
 from operator import mul
 A = len(a); out = [0.0] * (A + len(b) - 1)
 while A and a[A - 1] == 0.0: A -= 1
 B0, B1 = 0, len(b)
 while B0 < B1 and b[B0] == 0.0: B0 += 1
 while B1 > B0 and b[B1 - 1] == 0.0: B1 -= 1
 if not A or B0 >= B1: return out
 bb = b[B0:B1][::-1]; L = B1 - B0
 for k in range(B0, A + B1 - 1):
  kk = k - B0; i0 = max(0, kk - L + 1); i1 = min(kk, A - 1)
  out[k] = sum(map(mul, a[i0:i1 + 1], bb[L - 1 - kk + i0:L - kk + i1]))
 return out
@functools.lru_cache(None)
def window_rows(n, base, iA, mesh, kmax=KA_MAX):
 """rows[(jA, jB)]: measure/(jA! jB!) of jA A-type, jB B-type interval..."""
 h = (n - base) / mesh
 wA = [math.log1p(h / (base + i * h)) if i < iA else 0.0 for i in range(mesh)]
 wB = [math.log1p(h / (base + i * h)) if i >= iA else 0.0 for i in range(mesh)]
 PA = [[1.0]]
 for j in range(1, kmax + 1):
  PA.append([x / j for x in _conv(PA[-1], wA)])
 PB = [[1.0], wB]
 rows = {}
 for jA in range(kmax + 1):
  for jB in range(2):
   rows[(jA, jB)] = PA[jA] if jB == 0 else _conv(PA[jA], PB[1])
 return rows, h
@functools.lru_cache(None)
def ss_cum2(m, base, n, iA, mesh, alpha, trunc, kmax=KA_MAX):
 """lower sums L[jA][jB] (cumulative in jA, exact in jB) at integer m."""
 rows, h = window_rows(n, base, iA, mesh, kmax)
 out = [[0.0, 0.0] for _ in range(kmax + 1)]
 for (jA, jB), row in rows.items():
  j = jA + jB
  imax = len(row) - 1
  ifeas = min(imax, int(math.floor((m - j * h - j * base) / h))) if h > 0 else imax
  if ifeas < 0:
   continue
  imin = 0 if not trunc else max(0, ifeas - int(trunc * base / h) - 1)
  tab, ntab, ceil = RHO_TAB, len(RHO_TAB), math.ceil
  k0 = (m - j * base + alpha) / (base * RSTEP); hk = h / (base * RSTEP)
  s = 0.0
  for i in range(imin, ifeas + 1):
   weight = row[i]
   if weight:
    k = ceil(k0 - i * hk)
    s += weight * (1.0 if k <= 0 else (tab[k] if k < ntab else 0.0))
  out[jA][jB] = s
 cum = [[0.0, 0.0] for _ in range(kmax + 1)]
 for jB in range(2):
  acc = 0.0
  for jA in range(kmax + 1):
   acc += out[jA][jB]
   cum[jA][jB] = acc
 return tuple(tuple(x) for x in cum)
def ss_lower2(m, base, n, iA, kA, kB, mesh, alpha, trunc=3.0):
 """P(all primes <= 2^n, <= kA in (2^base, 2^nA], <= kB in (2^nA, 2^n])...))"""
 m = math.ceil(m)
 if n <= base:
  return rho_lo((m + alpha) / base)
 cum = ss_cum2(m, base, n, iA, mesh, alpha, trunc, kA)
 return cum[kA][0] + (cum[kA][1] if kB else 0.0)
def surv_upper(m, base, T):
 """P(product of the prime factors above 2^base of an m-bit integer <=...)"""
 if T <= 0:
  return 0.0
 if T >= m:
  return 1.0
 return min(1.0, max(1.0, T / base) * rho_hi((math.floor(m) - T) / base))
MONT5 = 367.76               # measured five-limb Montgomery multiplication (ecm12.c)
def ecm_instr(fbits, nbits, success):
 """ecm16 schedule instructions (KERNELS.md)"""
 import ecm16
 return ecm16.ecm_instr(fbits, nbits, success)
def log2li(l):
 x = l * LN2
 return l - math.log2(x) + math.log2(1 + 1 / x + 2 / x / x)
CELLS = {4:13315520,6:1573359464,8:47982017060,96:409193506988359583551012,192:1675659512501707896389090868,
  224:10654623198171176305341831664,256:52897551248440914796084637156,288:217401149762189752879411886484,
  320:769748406668707279529221545272,384:6863095027560050329435265332068}
class P:
 """the parameter set (proof.md / DESCENT.md)."""
 LB = 69.0; LX = 69.0; PMIN = 8.0; SLACK = 60.0
 ALPHA_F = 5.0; ALPHA_G = 1.0; LOG2U6 = 6 * 63.714
 QBINS = [(1 / 6, -26.4 + 23, 14.5 + 6), (1 / 6, -11.2 + 23, 14.5 + 6), (1 / 6, -0.1 + 23, 14.5 + 6),
   (1 / 6, 10.5 + 23, 14.5 + 6), (1 / 6, 23.2 + 23, 14.5 + 6), (1 / 6, 86.6 + 23, 14.5 + 6)]
 BOOT_BINS = [(0.5, 1.4 + 4), (0.4, 6.8 + 4), (0.1, 12.7 + 4)]
 QBINS_LO = [(-140.05 - 23, -54.56 - 6), (-26.4 - 23, -54.56 - 6), (-11.2 - 23, -54.56 - 6),
    (-0.1 - 23, -54.56 - 6), (10.5 - 23, -54.56 - 6), (23.2 - 23, -54.56 - 6)]
 LOG2_RELS = 63.97; MAIN_UF = 738.0 / 69; MAIN_UG = 497.9 / 69   # 2-SE lower bounds (audit)
 THIN_MAX = 90               # nodes with q <= 2^95 discard
 MISS_LO = 40
 UDUP = 0.25                 # distinct relations per raw relation (unit duplicates)
 # pr over 24 cf-quantile sub-cells and 10 boot deciles
 QBINS_SUB = [[-26.47,-0.76],[-15.7,1.36],[-8.64,2.41],[-3.29,4.1],[1.16,5.45],[5.07,4.44],[8.45,5.05],[11.67,5.09],[14.69,6.45],[17.52,5.98],[20.28,6.25],[22.96,6.93],[25.61,7.17],[28.26,8.98],[30.95,9.36],[33.69,9.07],[36.56,9.25],[39.6,10.43],[42.84,10.22],[46.5,11.73],[50.67,11.47],[55.75,12.05],[63.02,13.19],[109.6,20.5]]
 BOOT_SUB = [-13.24,-9.46,-7.11,-5.11,-3.15,-1.19,0.59,2.53,5.2,16.7]
 PRCAP_NODE = 4.0; PRCAP_BOOT = 1.5
 MISS_EPS = 2.0 ** -8        # Chebyshev margin on the realized missing set...
 MIN_RELS = 16.0; PSI_MAX = 0.8; PHI_FLOOR = 0.02
 TABLE_SLACK = 2.0 ** -16    # relative margin on the boot failure for the realized...
 NB_R1 = 0.5; NB_R2 = 1.0    # negative-binomial shapes
 BOUND_CELL = 2400.0         # certified per-cell norm bounds, both sides (84 forms)
 HITS_PER_CELL = 2 * math.log(69.0 / 8.0)
 SIEVE_HIT = 1100.0; SCAN = 10.0; ZERO = 2.0   # separate f and g score arrays
 PAIR_SETUP = 2.0 ** 25.1; LLL_NODE = 2.0 ** 36; NORMS = 2.0 ** 26; TDIV = 2.0 ** 24.5; IDEAL = 2.0 ** 34; PRIME = 2.0 ** 24
 ECM_SUCCESS = 0.999
 RADII = (4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256, 288, 320, 384)
 MESH = 384; TRUNC = 0.0
 LOG2P6 = 6 * 380.99700
 LMAX = 225; BIN = 5
def pairs_count(Pm):
 return 2.0 * 2 ** log2li(Pm.LX)
def node_cost(Pm, l_lo, l_hi, side, n, iA, kA, kB, R, thin=None):
 """own cost (instructions) and statistics of one node plan;"""
 pairs = pairs_count(Pm)
 cfmax = max(b[1] for b in Pm.QBINS); cgmax = max(b[2] for b in Pm.QBINS)
 S = math.log2(CELLS[R]); c = (l_hi + S) / 12
 h = (n - Pm.LX) / Pm.MESH if n > Pm.LX else 0.0
 nA = Pm.LX + iA * h if n > Pm.LX else Pm.LX
 T = kA * nA + kB * n + Pm.SLACK
 c_lo = (l_lo + S) / 12
 pr = 0.0; sv = 0.0
 for (w, cf, cg), (cf_lo, cg_lo) in zip(Pm.QBINS, Pm.QBINS_LO):
  nf = 72 * c + cf - (l_lo if side == 'f' else 0)
  ng = 12 * c + Pm.LOG2U6 + cg - (l_lo if side == 'g' else 0)
  pf = ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC)
  pg = ss_lower2(ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
  pr += w * pf * pg
  nf_lo = 72 * c_lo + cf_lo - (l_hi if side == 'f' else 0)
  ng_lo = 12 * c_lo + Pm.LOG2U6 + cg_lo - (l_hi if side == 'g' else 0)
  sv += w * surv_upper(nf_lo, Pm.LX, T) * surv_upper(ng_lo, Pm.LX, T)
 pri = 0.0
 for cf, cg in Pm.QBINS_SUB:
  nf = 72 * c + cf - (l_lo if side == 'f' else 0)
  ng = 12 * c + Pm.LOG2U6 + cg - (l_lo if side == 'g' else 0)
  pri += ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC) * ss_lower2(ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
 pr = min(pri / len(Pm.QBINS_SUB), Pm.PRCAP_NODE * pr)
 if pr <= 0:
  return None
 nf_max = 72 * c + cfmax - (l_lo if side == 'f' else 0)
 ng_max = 12 * c + Pm.LOG2U6 + cgmax - (l_lo if side == 'g' else 0)
 det = 1.0 - (1 - Pm.ECM_SUCCESS) * ((nf_max + ng_max) / 20.0 + 2 * (kA + kB))
 miss = missing(Pm, nf_max, ng_max)
 if thin is None or l_lo <= Pm.THIN_MAX:
  thin = thin or l_lo <= Pm.THIN_MAX
 mu = CELLS[R] * pr * Pm.UDUP * det * (p_none(Pm, nf_max, ng_max) if thin else 1.0)
 if mu < Pm.MIN_RELS:
  return None
 surv_cost = (Pm.NORMS + Pm.TDIV + ecm_instr(Pm.LX, nf_max, Pm.ECM_SUCCESS) + ecm_instr(Pm.LX, ng_max, Pm.ECM_SUCCESS)
    + 2 * (ecm_instr(nA, T, Pm.ECM_SUCCESS) + Pm.PRIME) + Pm.IDEAL)
 survivors = CELLS[R] * sv + CELLS[R] * pr
 sieve = CELLS[R] * (Pm.HITS_PER_CELL * Pm.SIEVE_HIT + Pm.SCAN + Pm.ZERO + Pm.BOUND_CELL) + pairs * Pm.PAIR_SETUP
 trials = CELLS[R] * surv_cost
 if trials < sieve + survivors * surv_cost:
  own = Pm.LLL_NODE + trials; how = 'ecm'; survivors = CELLS[R]
  ret, loss_f, loss_g = 1.0, 0.0, 0.0
 else:
  own = Pm.LLL_NODE + sieve + survivors * surv_cost; how = 'sieve'
  import retain16
  ret, loss_f, loss_g = retain16.retention(nf_max, ng_max, n, iA, kA, kB, Pm.SLACK)
  mu *= ret                 # relation-conditioned retention (DESCENT.md)
  if mu < Pm.MIN_RELS:
   return None
 return dict(own=own, n=n, iA=iA, nA=nA, kA=kA, kB=kB, R=R, S=S, mu=mu, how=how, miss=miss, ret=ret, loss_f=loss_f, loss_g=loss_g,
    miss_children=0.0 if thin else miss, thin=bool(thin), pr_log2=math.log2(pr),
    surv_log2=math.log2(survivors), own_log2=math.log2(own), nf=nf_max, ng=ng_max, det=det, T=T,
    sieve_log2=math.log2(sieve), survcost_log2=math.log2(survivors * surv_cost))
def bins_of(Pm):
 bins = [(int(Pm.LB), int(Pm.LB) + Pm.BIN + 1)]
 l = int(Pm.LB) + Pm.BIN + 1
 while l < Pm.LMAX:
  bins.append((l, l + Pm.BIN)); l += Pm.BIN
 return bins
def miss_nodes(Pm, PLAN):
 """missing-log ideal Q (2^MISS_LO <= N(Q) < 2^69)"""
 out = []
 for lo in range(Pm.MISS_LO, int(Pm.LB), 5):
  hi = min(lo + 5, Pm.LB)
  for side in ('f', 'g'):
   n, iA, kA, kB, R, thin = PLAN[(int(Pm.LB), side)]
   d = node_cost(Pm, lo, hi, side, float(n), iA, kA, kB, R, True)
   if d is None or kA or kB:
    raise ValueError('infeasible missing-log node %s %s' % (lo, side))
   out.append(dict(lo=lo, side=side, W=d['own'], phi=phi_nb(d['mu'], Pm, 0), mu=d['mu'], how=d['how'], nf=d['nf'], ng=d['ng']))
 return out
def phi_nb(mu, Pm, window, psi=0.0):
 """P(no relation whose children all succeed)"""
 r = Pm.NB_R1 if window <= 1 else Pm.NB_R2
 return max((1.0 + mu * (1.0 - psi) / r) ** (-r), Pm.PHI_FLOOR)
def evaluate(Pm, PLAN, BOOT):
 """evaluate the shipped plan"""
 shared = ('LX', 'MESH', 'TRUNC', 'ALPHA_F', 'ALPHA_G')
 if any(getattr(Pm, name) != getattr(P, name) for name in shared):
  raise ValueError('descent and retention parameters must agree')
 import retain16
 retain16.loss.cache_clear()
 bins = bins_of(Pm)
 W = {}; PHI = {}; DET = {}
 def wmax(side, n):
  return max([W[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])
 def phimax(side, n):
  return max([PHI[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])
 MISS = miss_nodes(Pm, PLAN)
 def miss_terms():
  return max(x['phi'] for x in MISS), max(x['W'] for x in MISS)
 for b in bins:
  for side in ('f', 'g'):
   n, iA, kA, kB, R, thin = PLAN[(b[0], side)]
   d = node_cost(Pm, b[0], b[1], side, float(n), iA, kA, kB, R, bool(thin))
   if d is None:
    raise ValueError('infeasible plan entry %s %s' % (b, side))
   nA = d['nA']
   pm, wm = miss_terms()
   psi = kA * (phimax('f', nA) + phimax('g', nA)) + kB * (phimax('f', n) + phimax('g', n)) + d['miss_children'] * pm
   if psi > Pm.PSI_MAX:
    raise ValueError('psi too large at %s %s' % (b, side))
   retry = 1.0 / (1.0 - psi)
   sub = kA * (wmax('f', nA) + wmax('g', nA)) + kB * (wmax('f', n) + wmax('g', n)) + d['miss_children'] * wm
   W[(b, side)] = d['own'] + retry * sub
   window = (kA if nA > Pm.LX + 1e-9 else 0) + (kB if n > nA + 1e-9 else 0)
   PHI[(b, side)] = phi_nb(d['mu'], Pm, window, psi)
   DET[(b, side)] = dict(d, psi=psi, retry=retry, sub_log2=math.log2(sub) if sub > 0 else None,
        W_log2=math.log2(W[(b, side)]), phi=PHI[(b, side)])
 n1, iA, kA, Rb = BOOT
 if n1 > Pm.LMAX:
  raise ValueError('boot window above the last bin')
 Sb = math.log2(CELLS[Rb]); cb = (Pm.LOG2P6 + Sb) / 12
 h = (n1 - Pm.LX) / Pm.MESH; nA = Pm.LX + iA * h
 m_max = 12 * cb + Pm.LOG2U6 + max(x[1] for x in Pm.BOOT_BINS)
 pr = sum(w * ss_lower2(12 * cb + Pm.LOG2U6 + cg, Pm.LX, float(n1), iA, kA, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
   for w, cg in Pm.BOOT_BINS)
 pri = sum(ss_lower2(12 * cb + Pm.LOG2U6 + cg, Pm.LX, float(n1), iA, kA, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC) for cg in Pm.BOOT_SUB)
 pr = min(pri / len(Pm.BOOT_SUB), Pm.PRCAP_BOOT * pr)
 det = 1.0 - (1 - Pm.ECM_SUCCESS) * (m_max / 20.0 + kA + 1)
 import retain16
 trials = Pm.MIN_RELS / (pr * det * Pm.UDUP * retain16.counted_fraction())
 balls = 1.0 + trials / CELLS[Rb]
 trial = Pm.NORMS + ecm_instr(nA, m_max, Pm.ECM_SUCCESS) + Pm.PRIME + Pm.IDEAL
 pm, wm = miss_terms(); mb = missing(Pm, 0, m_max)
 psi = kA * phimax('g', nA) + phimax('g', float(n1)) + mb * pm
 if psi > Pm.PSI_MAX:
  raise ValueError('boot psi too large')
 retry = 1.0 / (1.0 - psi)
 sub = kA * wmax('g', nA) + wmax('g', float(n1)) + mb * wm
 one_log = Pm.LLL_NODE * balls + trials * trial + retry * sub
 k = int(Pm.MIN_RELS)
 phi_boot = psi ** k * (1 + Pm.TABLE_SLACK)   # all k kept relations fail
 shared_log2 = math.log2(k * (k - 1) / 2) - math.log2(balls)   # two kept relations from one lattice (same exponent)
 boot = dict(n1=n1, nA=nA, kA=kA, R=Rb, S=Sb, m=m_max, pr_log2=math.log2(pr), trials_log2=math.log2(trials),
    trial_log2=math.log2(trial), balls=balls, miss=mb, own_log2=math.log2(trials * trial), sub_log2=math.log2(sub),
    psi=psi, retry=retry, phi=phi_boot, shared_log2=shared_log2, one_log_log2=math.log2(one_log))
 return dict(W=W, PHI=PHI, DET=DET, bins=bins, boot=boot, miss=MISS, two_log_units=math.log2(2 * one_log / UNIT),
    one_log_instr=one_log)
_PLAN = {
 69:((69,48,0,0,256,1),(69,48,0,0,288,1)),
 75:((75,42,1,0,224,1),(75,48,1,0,224,1)),
 80:((75,48,1,0,288,1),(72,48,2,0,192,1)),
 85:((78,26,2,0,192,1),(75,47,2,0,192,1)),
 90:((75,46,2,0,224,1),(75,45,2,0,288,1)),
 95:((79,22,1,1,256,1),(80,25,1,1,256,0)),
 100:((84,13,1,1,288,1),(80,16,2,1,320,0)),
 105:((85,18,2,1,256,1),(85,17,3,1,256,0)),
 110:((109,25,2,0,256,0),(93,22,2,1,288,1)),
 115:((92,22,2,1,320,1),(95,25,3,1,224,1)),
 120:((95,28,2,1,288,1),(96,46,3,0,288,0)),
 125:((123,23,3,0,288,1),(121,23,4,0,288,1)),
 130:((99,24,3,1,320,1),(100,24,4,1,320,0)),
 135:((100,29,3,1,320,1),(110,30,3,1,320,1)),
 140:((99,36,3,1,320,0),(109,36,3,1,320,0)),
 145:((111,29,3,1,320,1),(110,35,4,1,320,1)),
 150:((110,36,3,1,288,0),(110,42,4,1,320,0)),
 155:((110,29,4,1,288,0),(125,30,4,1,320,1)),
 160:((115,27,4,1,320,0),(130,28,4,1,320,0)),
 165:((125,26,4,1,384,1),(128,32,5,1,320,1)),
 170:((125,30,4,1,384,1),(135,36,4,1,384,0)),
 175:((136,30,4,1,288,1),(150,30,4,1,384,0)),
 180:((130,36,4,1,320,0),(128,35,6,1,320,0)),
 185:((140,30,4,1,384,0),(136,32,6,1,384,1)),
 190:((145,30,4,1,384,0),(160,24,6,1,320,0)),
 195:((150,30,4,1,384,0),(159,39,4,1,384,0)),
 200:((150,36,4,1,320,0),(146,35,6,1,320,0)),
 205:((165,30,4,1,384,1),(158,30,6,1,384,0)),
 210:((160,36,4,1,320,0),(165,29,6,1,384,0)),
 215:((165,36,4,1,320,0),(170,36,6,1,224,0)),
 220:((170,36,4,1,384,1),(171,36,6,1,224,0)),
}
PLAN16 = {(b, s): (v[i][0], 8 * v[i][1], *v[i][2:])
     for b, v in _PLAN.items() for i, s in enumerate('fg')}
BOOT16 = (215, 368, 16, 4)
def report(res):
 """compact per-bin table for the experiment report."""
 rows = []
 for b in res['bins']:
  for side in ('f', 'g'):
   d = res['DET'][(b, side)]
   rows.append(dict(bin=b[0], side=side, n=int(d['n']), nA=round(d['nA'], 1), kA=d['kA'], kB=d['kB'], R=d['R'], ret=round(d['ret'], 4),
       loss_f=d['loss_f'], loss_g=d['loss_g'], how=d['how'], mu=round(d['mu'], 2), own_log2=round(d['own_log2'], 3),
       W_log2=round(d['W_log2'], 3), psi=round(d['psi'], 4), phi=round(d['phi'], 6), nf=round(d['nf'], 1), ng=round(d['ng'], 1), T=round(d['T'], 1)))
 return rows
