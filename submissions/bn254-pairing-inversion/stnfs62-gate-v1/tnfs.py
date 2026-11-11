"""Shared constants, exact arithmetic and the cost model (Python standard library only)."""
import json
import math
from pathlib import Path

U = 0x44E992B44A6909F1
P = 36 * U**4 + 36 * U**3 + 24 * U**2 + 6 * U + 1
R_ORDER = 36 * U**4 + 36 * U**3 + 18 * U**2 + 6 * U + 1
H = [-1, -2, 0, 0, 1, 0, 1]  # h(y) = y^6 + y^4 - 2y - 1, low degree first
LN2 = math.log(2.0)
C_UNIT = 20000.0


# ---------- Z[y]/h and F_q[y]/h arithmetic (coefficient lists of length 6) ----------
def hmul(a, b, mod=None):
    t = [0] * 11
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                t[i + j] += ai * bj
    for k in range(10, 5, -1):  # y^6 = -y^4 + 2y + 1
        c = t[k]
        if c:
            t[k] = 0
            t[k - 2] -= c
            t[k - 5] += 2 * c
            t[k - 6] += c
    t = t[:6]
    return [c % mod for c in t] if mod else t


def hpow(a, e, mod):
    r = [1, 0, 0, 0, 0, 0]
    while e:
        if e & 1:
            r = hmul(r, a, mod)
        a = hmul(a, a, mod)
        e >>= 1
    return r


def det_int(m):
    """Exact integer determinant (Bareiss)."""
    m = [row[:] for row in m]
    n = len(m)
    sign, prev = 1, 1
    for k in range(n - 1):
        if m[k][k] == 0:
            for i in range(k + 1, n):
                if m[i][k]:
                    m[k], m[i] = m[i], m[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                m[i][j] = (m[i][j] * m[k][k] - m[i][k] * m[k][j]) // prev
        prev = m[k][k]
    return sign * m[n - 1][n - 1]


def hnorm(e):
    """Norm of e in Z[y]/h down to Z: determinant of multiplication by e."""
    cols = []
    cur = e[:]
    for _ in range(6):
        cols.append(cur)
        cur = hmul(cur, [0, 1, 0, 0, 0, 0])
    return det_int([[cols[c][r] for c in range(6)] for r in range(6)])


Y1 = [0, 1, 0, 0, 0, 0]
Y2 = [0, 0, 1, 0, 0, 0]
Y3 = [0, 0, 0, 1, 0, 0]
Y4 = [0, 0, 0, 0, 1, 0]


def hscale(a, c):
    return [c * v for v in a]


def hadd(*xs):
    return [sum(v) for v in zip(*xs)]


def norm_f(a, b):
    """|N_{K_h/Q}(sum_i f_i a^i b^(8-i))| for f = 36x^8 + 36y x^6 + 24y^2 x^4 + 6y^3 x^2 + y^4."""
    a2, b2 = hmul(a, a), hmul(b, b)
    a4, b4 = hmul(a2, a2), hmul(b2, b2)
    a6, b6 = hmul(a4, a2), hmul(b4, b2)
    a8, b8 = hmul(a4, a4), hmul(b4, b4)
    s = hadd(hscale(a8, 36),
             hscale(hmul(Y1, hmul(a6, b2)), 36),
             hscale(hmul(Y2, hmul(a4, b4)), 24),
             hscale(hmul(Y3, hmul(a2, b6)), 6),
             hmul(Y4, b8))
    return abs(hnorm(s))


def norm_g(a, b):
    """|N_{K_h/Q}(a^2 - u y b^2)| for g = x^2 - u y."""
    a2, b2 = hmul(a, a), hmul(b, b)
    return abs(hnorm(hadd(a2, hscale(hmul(Y1, b2), -U))))


def log2_int(n):
    if n <= 0:
        return float("-inf")
    k = n.bit_length()
    if k <= 900:
        return math.log2(n)
    return math.log2(n >> (k - 900)) + (k - 900)


def small_primes(limit):
    s = bytearray([1]) * (limit + 1)
    s[0:2] = b"\0\0"
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(limit + 1) if s[i]]


PRIMES_1000 = small_primes(1000)
C1000 = sum(math.log(q) / (q - 1) for q in PRIMES_1000)  # nats: mean log of the 1000-smooth part of a random integer


def smooth_part_bits(n):
    """log2 of the part of n made of primes below 1000."""
    s = 0.0
    for q in PRIMES_1000:
        while n % q == 0:
            n //= q
            s += math.log2(q)
    return s


# ---------- Dickman rho ----------
class Dickman:
    def __init__(self, step=512, umax=40):
        self.st = step
        n = umax * step + 1
        rho = [1.0] * n
        run = float(step - 1)  # sum of rho[j], j = i-step+1 .. i-1, for i = step+1
        for i in range(step + 1, n):
            u = i / step
            # u rho(u) = int_{u-1}^{u} rho(t) dt  (trapezoid, implicit in rho[i])
            rho[i] = (rho[i - step] / 2 + run) / step / (u - 0.5 / step)
            run += rho[i] - rho[i - step + 1]
        self.rho = rho
        self.umax = umax

    def __call__(self, u):
        if u <= 1:
            return 1.0
        if u >= self.umax:
            return 0.0
        x = u * self.st
        k = int(x)
        fr = x - k
        return self.rho[k] * (1 - fr) + self.rho[k + 1] * fr


def ecm_mults(bits):
    """Modular multiplications of the ECM schedule charged for primes below 2^bits: 8 x L[1/2, sqrt 2]."""
    lp = bits * LN2
    return 8.0 * math.exp(math.sqrt(2.0 * lp * math.log(lp)))


def modmul_instr(bits):
    k = math.ceil(bits / 64)
    return 16 * k * k + 64


def fft_mul_instr(bits):
    """H-FFT: instructions charged for one product of two integers of at most `bits` bits (bits >= 2^16)."""
    return 8.0 * bits * math.log2(bits)


def batch_test_instr(bound_bits, number_bits):
    """Instructions charged per tested number by the product/remainder-tree smoothness test for primes below 2^bound_bits."""
    pbits = 1.4428 * 2.0**bound_bits            # bits of the product of all primes below the bound
    k = pbits / number_bits                     # numbers per batch
    products = bound_bits + 6 * math.log2(k) + 16
    per_number = products * fft_mul_instr(pbits) / k
    return per_number + (math.log2(number_bits) + 4) * modmul_instr(number_bits) + number_bits * number_bits, pbits / number_bits * per_number


BATCH_MAX_BOUND = 66


def mean_se(v):
    n = len(v)
    m = sum(v) / n
    var = sum((x - m) ** 2 for x in v) / (n - 1)
    return m, math.sqrt(var / n)


def l2(x):
    return math.log2(x)


def model(data, verbose=False):
    """Recompute every row of the ledger from the stored production-size samples."""
    D = Dickman()
    cfg = data["config"]
    b, A, slack = cfg["log2_B"], cfg["A"], cfg["slack_bits"]
    lnB = b * LN2
    out = {}
    # ---- main phase ----
    main = data["main"]
    n_all = main["draws"]
    # primitive draws: [nf, ng, m, Xf, Xg, sf, sg], stored in tenths of a bit
    rows = [[r[0] / 10, r[1] / 10, r[2], r[3] / 10, r[4] / 10, r[5] / 10, r[6] / 10] for r in main["rows"]]
    sf_m, sf_se = mean_se([r[5] for r in rows])
    sg_m, sg_se = mean_se([r[6] for r in rows])
    alpha_f = C1000 - (sf_m - 3 * sf_se) * LN2
    alpha_g = C1000 - (sg_m - 3 * sg_se) * LN2
    af, ag = max(alpha_f, 0.0), max(alpha_g, 0.0)
    nf_sorted = sorted(r[0] for r in rows)
    ng_sorted = sorted(r[1] for r in rows)
    nf_min = nf_sorted[len(rows) // 200]
    ng_min = ng_sorted[len(rows) // 200]
    ys, yraw = [], []
    for nf, ng, m, xf, xg, _, _ in rows:
        pr = D((nf * LN2 + af) / lnB) * D((ng * LN2 + ag) / lnB)
        ok = xf <= slack and xg <= slack and nf >= nf_min and ng >= ng_min
        ys.append(pr / m if ok else 0.0)
        yraw.append(pr)
    ys += [0.0] * (n_all - len(rows))
    yraw += [0.0] * (n_all - len(rows))
    y_m, y_se = mean_se(ys)
    yr_m, yr_se = mean_se(yraw)
    S = 2 * A + 1
    V = S**12 / 2.0
    n_ideals = 2.2 * 2.0**b / lnB
    dim = n_ideals + 29 + 64
    rel_low = V * (y_m - 3 * y_se) * 0.98
    rel_all_high = V * (yr_m + 3 * yr_se)
    hits = sum(main["small_ideal_inverse_sums"]) + 2 * 1.05 * math.log(b / 12.0)
    tf, tg = nf_min - slack - 1, ng_min - slack - 1
    p1f = (nf_sorted[-1] + 8 - tf) * LN2 * D(tf / b)
    p1g = (ng_sorted[-1] + 8 - tg) * LN2 * D(tg / b)
    per_cell = 32 * hits + 96 * n_ideals / S**6 + 8 + min(p1f, p1g) * 2**17
    collection = V * per_cell
    factor_rel = rel_all_high * (64 * n_ideals + 2**30 + 2**12)
    schirokauer = dim * 2.0**38
    factor_base = 2.0**b * 2**24 / lnB + 2.0**b * 64
    la = 22400.0 * dim * dim
    out.update(alpha_f=alpha_f, alpha_g=alpha_g, nf_min=nf_min, ng_min=ng_min,
               log2_cells=l2(V), log2_yield=l2(y_m), yield_rel_se=y_se / y_m,
               log2_relations_lower=l2(rel_low), log2_dim=l2(dim), relations_over_dim=rel_low / dim,
               hits_per_cell=hits, instr_per_cell=per_cell, retention=sum(1 for r in rows if r[3] <= slack and r[4] <= slack) / len(rows),
               mean_inverse_multiplicity=sum(1.0 / r[2] for r in rows) / len(rows), primitive_fraction=len(rows) / n_all,
               mean_nf=sum(r[0] for r in rows) / len(rows), mean_ng=sum(r[1] for r in rows) / len(rows))
    # ---- occurrence of factor-base ideals among the dim relations used ----
    u_f = (out["mean_nf"] * LN2 + af) / lnB
    u_g = (out["mean_ng"] * LN2 + ag) / lnB

    def kappa(u, t):
        return D(u - t) / D(u) if u - t > 0 else 1.0 / D(u)

    def missing(norm_bits, bound_bits, side_u):
        """Expected number of prime ideals of norm in (2^(b-12), 2^b] of a bound-smooth norm that never occurred."""
        uu = norm_bits / bound_bits
        tot = 0.0
        steps = 240
        for i in range(steps):
            lo = b - 12 + 12 * i / steps
            hi = lo + 12 / steps
            if hi > bound_bits:
                break
            dens = kappa(uu, hi / bound_bits) * math.log(hi / lo)
            lam = dim * kappa(side_u, lo / b) / 2.0**hi
            tot += dens * math.exp(-lam)
        return tot

    def children(norm_bits, bound_bits, lo_bits, hi_bits):
        """Expected number of prime factors with norm in (2^lo_bits, 2^hi_bits] of a bound-smooth norm (Dickman local density)."""
        hi_bits = min(hi_bits, bound_bits)
        if hi_bits <= lo_bits:
            return 0.0
        uu = norm_bits / bound_bits
        tot, steps = 0.0, 100
        for i in range(steps):
            lo = lo_bits + (hi_bits - lo_bits) * i / steps
            hi = lo_bits + (hi_bits - lo_bits) * (i + 1) / steps
            tot += kappa(uu, hi / bound_bits) * math.log(hi / lo)
        return tot

    def child_work(mf, mg, bound_bits, W):
        """Expected work of all children of one relation: each child is charged at its own bin."""
        tot, kids, lo = 0.0, 0.0, b
        for e in edges:
            if lo >= bound_bits:
                break
            k = children(mf, bound_bits, lo, e) + children(mg, bound_bits, lo, e)
            kids += k
            tot += k * W[e]
            lo = e
        return tot, kids

    # ---- descent nodes ----
    edges = data["descent"]["edges"]  # top edges of the special-q bins, in bits
    bins = data["descent"]["bins"]
    W = {}
    node_rows = []
    lower = b
    for n in edges:
        best = {}
        for side in (0, 1):
            cand = None
            for Rk, cells in sorted(bins[f"{side}:{n}"].items(), key=lambda kv: int(kv[0])):
                R = int(Rk)
                mf = mg = 0.0
                inv = []
                for ideal_cells in cells:  # one list of sampled ball cells per sampled special-q ideal
                    ps = []
                    for c in ideal_cells:
                        if not c:  # non-primitive draw: not counted
                            ps.append(0.0)
                            continue
                        nfb, ngb = c[0] / 10, c[1] / 10
                        if side == 0:
                            nfb -= n - 1  # the special-q ideal (norm >= 2^(n-1)) is removed from its side
                        else:
                            ngb -= n - 1
                        nfb = max(nfb, 1.0) + af / LN2
                        ngb = max(ngb, 1.0) + ag / LN2
                        ps.append(D(nfb / lower) * D(ngb / lower))
                        mf, mg = max(mf, nfb), max(mg, ngb)
                    pi = sum(ps) / len(ps)
                    inv.append(1.0 / pi if pi > 0 else float("inf"))
                if max(inv) == float("inf") or (2 * R + 1) ** 12 / max(inv) < 64:
                    continue
                spread = math.log2(max(inv) / min(inv))
                p = 1.0 / (1.5 * sum(inv) / len(inv))  # 1/p charged = 1.5 x mean over sampled ideals of 1/p(ideal)
                miss = missing(mf, lower, u_f) + missing(mg, lower, u_g)
                if miss >= 0.5:
                    continue
                p_use = p * (1 - miss) * 0.98 * 0.99 ** (2 * math.ceil(max(mf, mg) / b))
                if lower <= BATCH_MAX_BOUND:
                    tf_, batch_f = batch_test_instr(lower, mf)
                    tg_, batch_g = batch_test_instr(lower, mg)
                    trial = 2**13 + 2**17 + tf_ + tg_
                    extra = batch_f + batch_g + 2 * math.ceil(max(mf, mg) / b) * ecm_mults(lower) * modmul_instr(max(mf, mg))
                else:
                    trial = 2**13 + 2**17 + ecm_mults(lower) * (modmul_instr(mf) + modmul_instr(mg))
                    extra = 0.0
                cw, kids = child_work(mf, mg, lower, W)
                w = 2.0**44 + extra + trial / p_use + cw
                if cand is None or w < cand[0]:
                    cand = (w, R, p_use, kids, miss, mf, mg, spread)
            if cand is None:
                raise ValueError(f"no admissible radius for bin {side}:{n}")
            best[side] = cand
            node_rows.append(dict(side="fg"[side], top_bits=n, child_bound_bits=lower, radius=cand[1],
                                  log2_p=l2(cand[2]), children=cand[3], missing=cand[4], ideal_spread_bits=cand[7],
                                  log2_work=l2(cand[0])))
        W[n] = max(best[0][0], best[1][0])
        lower = n
    # ---- boot ----
    boot = {"rows": [[r[0] / 10, r[1] / 10, r[2]] for r in data["boot"]["rows"]]}
    n1 = edges[-1]
    bs_m, bs_se = mean_se([r[1] for r in boot["rows"]])
    alpha_b = max(C1000 - (bs_m - 3 * bs_se) * LN2, 0.0)
    pb = [D((r[0] * LN2 + alpha_b) / (n1 * LN2)) if r[2] else 0.0 for r in boot["rows"]]
    pb_m, pb_se = mean_se(pb)
    nb_max = max(r[0] for r in boot["rows"])
    miss_b = missing(nb_max, n1, u_g)
    cw_b, kids_b = child_work(nb_max, 0.0, n1, W)
    p_boot = (pb_m - 2 * pb_se) * (1 - miss_b) * 0.98 * 0.99 ** math.ceil(nb_max / b)
    boot_trial = 2.0**44 + 2.0**24 + 2**17 + ecm_mults(n1) * modmul_instr(nb_max + 8)
    one_log = boot_trial / p_boot + cw_b
    two_logs = 8 * 2 * one_log  # charged at eight times the model's expectation
    wrapper = 2.0**70
    fixed = collection + factor_rel + schirokauer + factor_base + la + wrapper
    run = fixed + two_logs
    gated = bool(cfg.get("gate"))
    gate = 13 / 32 if gated else 1.0
    total_units = (1 if gated else 0) + gate * run / C_UNIT
    out.update(nodes=node_rows, log2_boot_p=l2(p_boot), boot_children=kids_b, boot_missing=miss_b, alpha_boot=alpha_b,
               boot_norm_mean=sum(r[0] for r in boot["rows"]) / len(boot["rows"]), boot_norm_max=nb_max,
               ledger_log2_instructions=dict(
                   relation_collection=l2(collection), relation_factoring=l2(factor_rel), schirokauer_maps=l2(schirokauer),
                   factor_base=l2(factor_base), linear_algebra=l2(la), two_individual_logarithms_expected=l2(two_logs),
                   wrapper=l2(wrapper), run_total_expected=l2(run)),
               log2_run_units=l2(run / C_UNIT), descent_cap_log2_instructions=l2(32 * two_logs),
               gate="13/32" if gated else "none", post_gate_success=0.96, success=gate * 0.96,
               log2_expected_total_units=l2(total_units))
    return out


def load_data(base):
    """data_main.json (config, region sample), data_desc_f/g.json (special-q ball samples), data_boot.json."""
    base = Path(base)
    data = json.loads((base / "data_main.json").read_text())
    boot = json.loads((base / "data_boot.json").read_text())
    bins = {}
    for name in ("data_desc_f.json", "data_desc_g.json"):
        bins.update(json.loads((base / name).read_text()))
    data["descent"] = {"edges": boot["edges"], "bins": bins, "check_lattices": boot["check_lattices"]}
    data["boot"] = boot["boot"]
    return data
