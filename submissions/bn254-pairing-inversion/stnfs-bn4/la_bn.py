"""bn1 linear algebra, Schirokauer structure and filtering table (pure Python; BN254 only).
Kernel prices are the accepted BLS12-381 package's straight-line 4-limb kernels (instruction counts do not depend on the modulus value)."""
import math, json, base64, lzma, os
from math import log2, ceil, log
from fractions import Fraction as F
U = 4965661367192848881
P = 36*U**4 + 36*U**3 + 24*U**2 + 6*U + 1
R = 36*U**4 + 36*U**3 + 18*U**2 + 6*U + 1
H = [1, -1, -2, 5, -2, -1, 1]; W = [-2, 2, 2, -3, 0, 1]; FC = [36, 36, 24, 6, 1]
LN2 = log(2); UNIT = 20000.
LB = 56.; NBAR = 2.
W_ROW = 401.; C_ENTRY_K = 8.38; C_CHI = 200.; C_RED9 = 270.; C_FULL = 97.; C_PROJ = 60.; N_SM = 7
OWNERS = 1; BW_M, BW_N = 1024, 40
REPLAY12 = 219                   # filtsim at 2^12 orbits, seed 1: rows N after filtering (la_bn.replay)
NU = 0.03544                     # filtering ratio charged: max over 2^18, 2^19 of mean + 2 SE, rounded up (filtbn.b85)
FILT_EXCESS = 1.095; FILT_MEANS = (420., 473.)
FFT_A, FFT_B, FFT_XFER, FFT_LOCAL = 8.5, 1200., 504., 2**26
def log2_li(x): t = x*LN2; return x - log2(t) + log2(1 + 1/t + 2/t/t)
def log2_fb(lb=LB): return log2(NBAR/2*2**log2_li(lb) + 2**31) - 1        # tau-orbits U'
def fft_F(N): return FFT_A*N*log2(max(N, 2))**2 + FFT_B + (FFT_XFER*N if N > FFT_LOCAL else 0.)
def row_tail(send=15, combine=62, redc=65, part=267, owners=OWNERS):
    return (owners-1)*(send+18) + (owners-1)*(combine+2) + redc + 6 + 10*owners + 40 + owners*ceil((part+6)/5)
C_ROW = float(row_tail())
def generator(log2_N, m=None, n=None, exact=False):
    m = m or BW_M; n = n or BW_N; k = m + n; L = (2**log2_N + N_SM)*(1/m + 1/n) + 2; J = ceil(log2(L))
    if exact: return sum(2**j*(m+k)*k*k*(fft_F(18*l+18) + 640*l) for j in range(J) for l in [-(-L//2**j)]) + 2**J*m*k*k*(C_FULL + C_RED9)
    return log2(4*k**3*log2(L)*(fft_F(18*L) + 640*L)/UNIT)
def linear_algebra(nu, m=None, n=None, lb=LB):
    m = m or BW_M; n = n or BW_N; log2_N = log2_fb(lb) + log2(nu); c_prod = W_ROW*C_ENTRY_K + C_ROW + C_CHI
    per_n2 = c_prod*(n/m + 1) + (1/m + 1/n)*(c_prod + n*C_FULL + C_RED9) + C_PROJ*(n/m + 1); per_n1 = (m + n)*C_RED9
    la = log2((2**(2*log2_N)*per_n2 + 2**log2_N*per_n1)/UNIT); g = generator(log2_N, m, n)
    return dict(log2_N=log2_N, c_prod=c_prod, log2_la=la, log2_gen=g, log2_both=log2(2**la + 2**g), m=m, n=n)
def bw_grid(nu):
    """(m, n): m a power of two 2^5..2^13, n a multiple of 5 (the Krylov kernel advances five vectors) from 10 to 200, both >= N_SM."""
    best = None
    for e in range(5, 14):
        for n in range(10, 201, 5):
            r = linear_algebra(nu, 2**e, n)
            if best is None or r['log2_both'] < best['log2_both']: best = r
    return best
# ---- kernel bound arithmetic at the BN254 r (STORAGE: _row_red, _chi_mul, _row_redc) ----
def kernel_bounds():
    chi = pow(P, 3, R); chi_abs = min(chi, R - chi); mu = 2**327//R; K = -(-2**323//R)
    Dmax = 2**323 + K*R                                   # D = pos - neg + K r, pos, neg < 2^323
    slack = 1 + F(2**248, R) + F(Dmax, 2**327)           # D/r - q < slack (one Barrett step)
    mass = 2**323//(2*R)                                  # coefficient mass a (pos or neg) sum may carry with terms < 2r
    return dict(r_bits=R.bit_length(), r_limbs=-(-R.bit_length()//64), chi_sq=(chi*chi + 1) % R == 0, chi_bits=chi_abs.bit_length(),
                mu_bits=mu.bit_length(), K_r_ge=K*R >= 2**323, D_lt=Dmax < 2**324 + R, barrett_slack=float(slack), barrett_ok=slack < 2,
                q_bits=(Dmax//R).bit_length(), prod_bits=(Dmax >> 248).bit_length() + mu.bit_length(), mass_log2=log2(mass),
                row_mass_ok=mass > 255*2*W_ROW*64, mont256_ok=R % 2 == 1 and 2*R < 2**256, mont320_ok=R < 2**255,
                kron_ok=(BW_M + BW_N)*2**50*R*R < 2**576)
# ---- Schirokauer structure ----
def _trim(a):
    while a and a[-1] == 0: a = a[:-1]
    return a
def _pmod(a, b, q):
    a = a[:]; inv = pow(b[-1], -1, q)
    while len(a) >= len(b):
        c = a[-1]*inv % q; d = len(a) - len(b)
        if c:
            for i, x in enumerate(b): a[d+i] = (a[d+i] - c*x) % q
        a.pop()
    return _trim(a)
def _res(a, b, q):
    a = _trim([x % q for x in a]); b = _trim([x % q for x in b]); s = 1
    while b and len(b) > 1:
        da, db = len(a)-1, len(b)-1; rem = _pmod(a, b, q)
        if not rem: return 0
        s = s*(-1)**(da*db)*pow(b[-1], da-(len(rem)-1), q) % q; a, b = b, rem
    return s*pow(b[0], len(a)-1, q) % q if b else 0
def _pgcd(a, b, q):
    a, b = _trim([x % q for x in a]), _trim([x % q for x in b])
    while b: a, b = b, _pmod(a, b, q)
    return a
def _wpow(k):
    r = [1]
    for _ in range(k):
        c = [0]*(len(r) + 5)
        for i, x in enumerate(r):
            for j, y in enumerate(W): c[i+j] += x*y
        for i in range(len(c)-1, 5, -1):
            for j in range(6): c[i-6+j] -= c[i]*H[j]
        r = c[:6]
    return r + [0]*(6 - len(r))
F_TERMS = [(8 - 2*k, [FC[k]*x for x in _wpow(k)]) for k in range(5)]      # f = sum FC[k] w^k x^(8-2k)
G_TERMS = [(2, [1, 0, 0, 0, 0, 0]), (0, [-U*x for x in W])]               # g = x^2 - u w
def absolute_poly(terms, deg, q=R):
    """Res_y(h(y), poly(x, y)) mod q by interpolation at 6 deg + 1 points (Newton form)."""
    n = 6*deg + 1; xs = list(range(n)); ys = [_res(H, [sum(c[i]*pow(k, d, q) for d, c in terms) % q for i in range(6)], q) for k in xs]
    for j in range(1, n):
        for i in range(n-1, j-1, -1): ys[i] = (ys[i] - ys[i-1])*pow(xs[i] - xs[i-j], -1, q) % q
    out = [ys[-1]]
    for i in range(n-2, -1, -1): out = [(-xs[i]*out[0] + ys[i]) % q] + [(out[k-1] - xs[i]*out[k]) % q for k in range(1, len(out))] + [out[-1]]
    return out
def unramified():
    out = {}
    for name, terms, deg in ('f', F_TERMS, 8), ('g', G_TERMS, 2):
        a = absolute_poly(terms, deg); d = [i*c % R for i, c in enumerate(a)][1:]
        out['deg_'+name] = len(_trim(a)) - 1; out['lead_'+name] = a[-1]
        out['sqf_'+name] = a[-1] != 0 and len(a) == 6*deg + 1 and len(_pgcd(a, d, R)) == 1
    out['lead_ok'] = out['lead_f'] == pow(36, 6, R) and out['lead_g'] == 1
    out['unramified'] = out['sqf_f'] and out['sqf_g']; return out
def signatures():
    """Sturm counts at the real roots of h (isolated to 1e-30 in exact rationals)."""
    def sturm(p):
        s = [_trim([F(c) for c in p])]; s.append([i*c for i, c in enumerate(s[0])][1:])
        while len(s[-1]) > 1:
            a, b = s[-2][:], s[-1]
            while len(a) >= len(b):
                q = a[-1]/b[-1]; d = len(a) - len(b); a = _trim([x - q*b[i-d] if i >= d else x for i, x in enumerate(a)][:-1])
            if not a: break
            s.append([-x for x in a])
        return s
    def ch(s, x): v = [t for t in (sum(c*x**i for i, c in enumerate(q)) for q in s) if t]; return sum((a > 0) != (b > 0) for a, b in zip(v, v[1:]))
    real = lambda p, B=F(10)**60: (lambda s: ch(s, -B) - ch(s, B))(sturm(p)); sh = sturm(H); ys = []
    def iso(a, b):
        n = ch(sh, a) - ch(sh, b)
        if n == 1 and b - a < F(1, 10**30): ys.append((a + b)/2)
        elif n: iso(a, (a + b)/2); iso((a + b)/2, b)
    iso(F(-10), F(10)); rf = rfp = rg = 0; wv = []
    for y in ys:
        w = sum(F(c)*y**i for i, c in enumerate(W)); wv.append(float(w))
        fx = [0]*9
        for k in range(5): fx[8 - 2*k] = FC[k]*w**k
        rf += real(fx); rfp += real(fx[0::2]); rg += real([-U*w, 0, 1])
    rk = lambda n, r1: r1 + (n - r1)//2 - 1; r1h = len(ys)
    ranks = dict(f=rk(48, rf), fp=rk(24, rfp), g=rk(12, rg), h=rk(6, r1h)); odd = ranks['f'] - ranks['fp'] + ranks['g'] - ranks['h']
    return dict(h=(r1h, (6 - r1h)//2), f=(rf, 24 - rf//2), fp=(rfp, 12 - rfp//2), g=(rg, 6 - rg//2), w_real=wv, ranks=ranks, odd=odd, sm=odd//2,
                sm_f=(ranks['f'] - ranks['fp'])//2, sm_g=(ranks['g'] - ranks['h'])//2)
def schirokauer_cost(log2_rels=51.):
    """Bound per relation: (phi^eps - 1)/r mod r in (Z/r^2)[y,x]/(h, f/36) (rank 48) and /(h, g) (rank 12), left-to-right powering,
    eight-limb coefficients (r^2 < 2^508). Products are schoolbook 8 x 8 limbs at <= 512 instructions, Barrett <= 1,024."""
    limbs = -(-(2*R.bit_length())//64); assert limbs == 8
    eps_f, eps_g = 48*R.bit_length(), 12*R.bit_length()
    sq_f = 48*48*512                      # all coefficient products of a squaring (no symmetry used)
    red_f = 7*4*36*512                    # 7 high x-terms, 4 lower terms of the monic f/36 (even powers), 6 x 6 coefficient products each
    red_h = 2**17                         # reductions by the monic h: small-integer multiples of eight-limb words (allowance)
    bar_f = 3*48*1024                     # three Barrett passes over the 48 coefficients
    phi_f = 2*48*6*8*40 + 48*1024         # times a - b x: 6 x 6 products of 64-bit by eight-limb words (40 per limb), two terms; one Barrett pass
    bit_f = sq_f + red_f + red_h + bar_f + phi_f
    sq_g = 12*12*512; red_g = 1*6*6*512 + 2**15; bar_g = 3*12*1024; phi_g = 2*12*6*8*40 + 12*1024; bit_g = sq_g + red_g + bar_g + phi_g
    per = eps_f*bit_f + eps_g*bit_g
    return dict(eps_bits=(eps_f, eps_g), bit_f_log2=log2(bit_f), bit_g_log2=log2(bit_g), per_rel_log2=log2(per), log2_rels=log2_rels,
                log2_units=log2(per) + log2_rels - log2(UNIT))
# ---- filtering table ----
def filt_table():
    here = os.path.dirname(os.path.abspath(__file__))
    d = json.loads(lzma.decompress(base64.b85decode(open(os.path.join(here, 'filtbn.b85')).read().strip())))
    return {int(k): tuple(map(tuple, v)) for k, v in d['FILT_RUNS'].items()}, d
def filt_stats(tab):
    out = {}
    for k, vs in sorted(tab.items()):
        x = [v[1]/2**k for v in vs]; n = len(x); m = sum(x)/n; sd = (sum((y - m)**2 for y in x)/(n - 1))**.5 if n > 1 else float('nan')
        out[k] = dict(n=n, mean=m, se=sd/n**.5, sd_sqrtU=sd*2**(k/2), maxc=max(v[2] for v in vs))
    return out
def charged_ratio(st):
    ks = sorted(st)[-2:]; return max(st[k]['mean'] + 2*st[k]['se'] for k in ks), ks
def replay(seed=1, k=12):
    import filtsim
    Us, res = filtsim.run(LB, FILT_MEANS[0]/LB, FILT_MEANS[1]/LB, FILT_EXCESS, seed, k, kmax=64, targets=(400,), verbose=False, c4=True); return res[400][1]
if __name__ == '__main__':
    import sys
    print(json.dumps(dict(kernel=kernel_bounds(), unram=unramified(), sig=signatures(), sm=schirokauer_cost(), row_tail=C_ROW, log2_Uprime=log2_fb()), default=str, indent=1))
