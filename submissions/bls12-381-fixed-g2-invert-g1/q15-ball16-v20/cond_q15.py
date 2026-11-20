"""Descent conditioning evidence (DESCENT.md [A], [B]): recompute every reported statistic from the
shipped records (cond_q15.b85: lzma, base85, JSON) and replay sample records from scratch.

[A] production (f_9, g_3): encountered children vs sampled ideals of equal side and bit size.
[B] same field Q(zeta_15), toy pair f_3 = X^3 - (2 + w), g_3 = X^3 + 2 w^2, real factor base of all
    G_12-orbit ideals of norm <= 2^14 kept by a 2-core filter; one-level descents.
"""
import json, math
from pathlib import Path

HERE = Path(__file__).resolve().parent
LX_B = 14
B_FB = 1 << LX_B
SMALL_NS_MAX = 128


def _prod(xs):
    xs = list(xs) or [1]
    while len(xs) > 1:
        xs = [xs[i] * xs[i + 1] if i + 1 < len(xs) else xs[i] for i in range(0, len(xs), 2)]
    return xs[0]


def _strip(n, p1, p2):
    for P in (p1, p2):
        while n > 1:
            g = math.gcd(n, P)
            if g == 1:
                break
            n //= g
    return n


def _w_roots(p):
    for g in range(2, p):
        w = pow(g, (p - 1) // 3, p)
        if w != 1:
            w2 = (w * w) % p
            return (w, w2) if w < w2 else (w2, w)
    raise RuntimeError(p)


def _ks(x, y):
    x = sorted(x); y = sorted(y); i = j = 0; d = 0.0; n, m = len(x), len(y)
    while i < n and j < m:
        if x[i] <= y[j]:
            i += 1
        else:
            j += 1
        d = max(d, abs(i / n - j / m))
    ne = n * m / (n + m); lam = (math.sqrt(ne) + 0.12 + 0.11 / math.sqrt(ne)) * d
    p = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 101))
    return d, max(0.0, min(1.0, p))


class Toy:
    """The [B] tower: tables, 2-core sets, ball geometry, one ball."""

    def __init__(self, Z, CM, ret_f, ret_g):
        self.Z, self.CM = Z, CM
        self.w = Z["y_pow"](5); w2 = Z["y_pow"](10)
        self.c_f = Z["add"]([2] + [0] * 7, self.w); self.c_g = Z["scale"](2, w2)
        self.sig = [Z["y_pow"]((7 * j) % 15) for j in range(8)]
        self.split, self.small = {}, {}
        for p in CM.primes_upto(B_FB):
            if p % 15 == 1:
                ws = _w_roots(p)
                self.split[p] = [(oi, w, pow((2 + w) % p, (p - 1) // 3, p) == 1, pow((-2 * w * w) % p, (p - 1) // 3, p) == 1) for oi, w in enumerate(ws)]
            elif p <= SMALL_NS_MAX:
                self.small[p] = _w_roots(p) if p % 3 == 1 else None
        self.sps = sorted(self.small); self.pps = sorted(self.split)
        self.P1 = _prod(self.sps + [p for p in self.pps if p <= 256]); self.P2 = _prod([p for p in self.pps if p > 256])
        self.ret_f = {tuple(x) for x in ret_f}; self.ret_g = {tuple(x) for x in ret_g}
        self.als = [[1, 0, 0, 0, 0, 0, 0, 0]]
        for j in range(1, 8):
            for s in (1, -1):
                v = [1, 0, 0, 0, 0, 0, 0, 0]; v[j] = s; self.als.append(v)
        for i in range(1, 8):
            for j in range(i + 1, 8):
                v = [1, 0, 0, 0, 0, 0, 0, 0]; v[i] = 1; v[j] = 1; self.als.append(v)
        self.bts = []
        for i in range(8):
            for s in (1, -1):
                v = [0] * 8; v[i] = s; self.bts.append(v)
        for i in range(8):
            for j in range(i + 1, 8):
                for si in (1, -1):
                    for sj in (1, -1):
                        v = [0] * 8; v[i] = si; v[j] = sj; self.bts.append(v)
        assert len(self.als) == 36 and len(self.bts) == 128
        self.units = [Z["y_pow"](k) for k in range(15)]

    def _sigma(self, v):
        out = [0] * 8
        for j, c in enumerate(v):
            if c:
                for i in range(8):
                    out[i] += c * self.sig[j][i]
        return out

    def _nzw(self, el):
        m = self.Z["mul"]; s1 = self._sigma(el); s2 = self._sigma(s1); s3 = self._sigma(s2)
        r = m(m(el, s1), m(s2, s3)); assert r[1] == r[2] == r[3] == r[4] == r[6] == r[7] == 0
        return r[0], r[5]

    def _cols(self, el_f, el_g, nf, ng, q, side):
        out = []
        for sd, rem, el, full in (("f", nf, el_f, nf * q if side == "f" else nf), ("g", ng, el_g, ng * q if side == "g" else ng)):
            c0, c1 = self._nzw(el); assert c0 * c0 - c0 * c1 + c1 * c1 == full
            cols = set()
            for p in self.sps:
                if rem % p == 0:
                    while rem % p == 0:
                        rem //= p
                    ws = self.small[p]
                    if ws is None:
                        cols.add((p, 0))
                    else:
                        h = [oi for oi, w in enumerate(ws) if (c0 + c1 * w) % p == 0]; assert h
                        cols.update((p, oi) for oi in h)
            hits = []
            for p in self.pps:
                if rem == 1:
                    break
                if rem % p == 0:
                    while rem % p == 0:
                        rem //= p
                    hits.append(p)
            if rem != 1:
                return None
            for p in hits:
                h = 0
                for oi, w, rf, rg in self.split[p]:
                    if (c0 + c1 * w) % p == 0:
                        if not (rf if sd == "f" else rg):
                            return None
                        cols.add((p, oi)); h += 1
                if h == 0:
                    return None
            out.append(cols)
        return out

    def canon(self, a, b):
        """Representative of (a, b) up to the torsion units +-y^k (k < 15)."""
        best = None
        for u in self.units:
            ua = self.Z["mul"](u, a); ub = self.Z["mul"](u, b)
            for s in (1, -1):
                t = tuple(s * x for x in ua + ub)
                if best is None or t < best:
                    best = t
        return best

    def ball(self, q, side, ry, rt):
        """Distinct usable relations M (whole ball; one per torsion-unit class), M9 (first 9 of 36 alpha rows), survivors,
        model means."""
        Z, CM = self.Z, self.CM
        B = [[0] * 16 for _ in range(16)]; B[0][0] = q; ryp = [pow(ry, i, q) for i in range(8)]
        for i in range(1, 8):
            B[i][0] = (-ryp[i]) % q; B[i][i] = 1
        for i in range(8):
            B[8 + i][0] = (rt * ryp[i]) % q; B[8 + i][8 + i] = 1
        R = CM.lll_16d_int(B, 99, 100)[0]
        a1, b1 = R[0][:8], R[0][8:]
        a2, b2 = next((v[:8], v[8:]) for v in R[1:] if abs(Z["norm"](Z["sub"](Z["mul"](a1, v[8:]), Z["mul"](b1, v[:8])))) == q)
        m = m4 = sv = 0; nfl, ngl, nfq, ngq = [], [], [], []; seen = set()
        for ia, al in enumerate(self.als):
            a1a = Z["mul"](al, a1); b1a = Z["mul"](al, b1)
            for bt in self.bts:
                a = Z["add"](a1a, Z["mul"](bt, a2)); b = Z["add"](b1a, Z["mul"](bt, b2))
                a3 = Z["pow"](a, 3); b3 = Z["pow"](b, 3)
                el_g = Z["add"](a3, Z["mul"](self.c_g, b3)); el_f = Z["sub"](a3, Z["mul"](self.c_f, b3))
                nf = abs(Z["norm"](el_f)); ng = abs(Z["norm"](el_g))
                if side == 0:
                    assert nf % q == 0; nf //= q
                else:
                    assert ng % q == 0; ng //= q
                if nf == 0 or ng == 0:
                    continue
                nfl.append(math.log2(nf)); ngl.append(math.log2(ng))
                if ia < 9:
                    nfq.append(nfl[-1]); ngq.append(ngl[-1])
                if _strip(nf, self.P1, self.P2) == 1 and _strip(ng, self.P1, self.P2) == 1:
                    sv += 1
                    c = self._cols(el_f, el_g, nf, ng, q, "f" if side == 0 else "g")
                    if c is not None and c[0] <= self.ret_f and c[1] <= self.ret_g:
                        k = self.canon(a, b)
                        if k not in seen:
                            seen.add(k); m += 1; m4 += ia < 9
        rf = len([1 for p, o in self.ret_f if p % 15 == 1]) / sum(1 for p in self.split for x in self.split[p] if x[2])
        rg = len([1 for p, o in self.ret_g if p % 15 == 1]) / sum(1 for p in self.split for x in self.split[p] if x[3])

        def mod(fl, gl):
            fl = sorted(fl); gl = sorted(gl); n = len(fl); pr = 0.0
            for s in range(6):
                i = min(n - 1, int((s + 0.5) * n / 6)); ef = fl[i] + 3.0; eg = gl[i] + 3.0
                pr += CM.rho_db(ef / LX_B, ef) * CM.rho_db(eg / LX_B, eg)
            return n * pr / 6.0 * rf * rg
        # 4230 of the 4608 cells are distinct up to torsion units (structural; the first 9 rows are all distinct)
        return m, m4, sv, round(mod(nfl, ngl) * 4230 / 4608, 4), round(mod(nfq, ngq), 4)


def stats_b(recs):
    """Failure = no usable relation. Model: negative binomial r = 1 at the model mean (no fitted input)."""
    out = {}
    for lab, mi, ui in (("child", 1, 3), ("ctrl", 5, 7), ("child_q", 2, 4), ("ctrl_q", 6, 8)):
        ph = lambda u: 1.0 / (1.0 + u)
        ms = [r[mi] for r in recs]; n = len(ms); mean = sum(ms) / n
        byp = {}
        for r in recs:
            byp.setdefault(r[0], []).append((r[mi] == 0, ph(r[ui])))
        obs = sum(f for g in byp.values() for f, _ in g); mod = sum(p for g in byp.values() for _, p in g)
        co = cm = npair = 0
        for g in byp.values():
            for i, (fi, _) in enumerate(g):
                if fi:
                    for j, (fj, pj) in enumerate(g):
                        if i != j:
                            npair += 1; co += fj; cm += pj
        out[lab] = dict(n=n, mean_M=round(mean, 2), mean_mu=round(sum(r[ui] for r in recs) / n, 2), fail=obs,
                        fail_model=round(mod, 1), after_sib_fail=f"{co}/{npair}", after_sib_model=round(cm, 1))
        assert obs <= mod and co <= cm
    return out


def stats_a(enc, ctrl):
    out = {}
    for k, i in (("cf", 1), ("cg", 2)):
        x = [e[i] for e in enc]; y = [c[i] for c in ctrl]
        mx, my = sum(x) / len(x), sum(y) / len(y)
        se = math.sqrt(sum((v - mx) ** 2 for v in x) / len(x) ** 2 + sum((v - my) ** 2 for v in y) / len(y) ** 2)
        d, p = _ks(x, y)
        out[k] = dict(diff=round(mx - my, 2), se=round(se, 2), ks_p=round(p, 2))
    out["n"] = [len(enc), len(ctrl)]
    out["max_cf"] = max(e[3] for e in enc + ctrl); out["max_cg"] = max(e[4] for e in enc + ctrl)
    return out


def replay_a(rec, Z, CM, U, s8):
    """Recompute one [A] record: lattice of the degree-1 ideal (q, y -> ry, X -> rt), 16 fixed points."""
    import random
    side, q, ry, rt = rec[:4]
    B = [[0] * 16 for _ in range(16)]; B[0][0] = q; ryp = [pow(ry, i, q) for i in range(8)]
    for i in range(1, 8):
        B[i][0] = (-ryp[i]) % q; B[i][i] = 1
    for i in range(8):
        B[8 + i][0] = (rt * ryp[i]) % q; B[8 + i][8 + i] = 1
    R = CM.lll_16d_int(B, 99, 100)[0]
    a1, b1 = R[0][:8], R[0][8:]
    a2, b2 = next((v[:8], v[8:]) for v in R[1:] if abs(Z["norm"](Z["sub"](Z["mul"](a1, v[8:]), Z["mul"](b1, v[:8])))) == q)
    rp = random.Random(20261120); pts = []
    while len(pts) < 16:
        al = [rp.randint(-4, 4) for _ in range(8)]; bt = [rp.randint(-4, 4) for _ in range(8)]
        if any(al) and sum(x * x for x in al) <= 64 and sum(x * x for x in bt) <= 64:
            pts.append((al, bt))
    w = Z["y_pow"](5); wp1 = Z["add"](w, [1] + [0] * 7); w2 = Z["y_pow"](10); lu8 = 8.0 * math.log2(abs(U))
    c = (math.log2(q) + s8) / 16.0; cf = []; cg = []
    for al, bt in pts:
        a = Z["add"](Z["mul"](al, a1), Z["mul"](bt, a2)); b = Z["add"](Z["mul"](al, b1), Z["mul"](bt, b2))
        a3, b3 = Z["pow"](a, 3), Z["pow"](b, 3); a6, b6 = Z["mul"](a3, a3), Z["mul"](b3, b3)
        ng = abs(Z["norm"](Z["sub"](a3, Z["scale"](U, Z["mul"](w2, b3)))))
        fe = Z["sub"](Z["add"](Z["add"](Z["mul"](a6, a3), Z["mul"](wp1, Z["mul"](a6, b3))), Z["mul"](a3, b6)), Z["mul"](w, Z["mul"](b6, b3)))
        nf = abs(Z["norm"](fe))
        assert (ng if side == "g" else nf) % q == 0
        cf.append(math.log2(nf) - 72.0 * c); cg.append(math.log2(ng) - 24.0 * c - lu8)
    return round(sum(cf) / 16, 3), round(sum(cg) / 16, 3)


def verify(Z, CM, U):
    import lzma, base64
    D = json.loads(lzma.decompress(base64.b85decode("".join((HERE / "cond_q15.b85").read_text().split()))))
    A, Bc = D["A"], D["B"]
    sa = stats_a(A["enc"], A["ctrl"]); sb = stats_b(Bc["records"])
    assert sa == A["summary"] and sb == Bc["summary"], "summary mismatch"
    na = 0
    for rec in A["live"]:
        assert replay_a(rec, Z, CM, U, A["s8"]) == (rec[4], rec[5]); na += 1
    toy = Toy(Z, CM, Bc["ret_f"], Bc["ret_g"]); nb = 0
    for q, side, ry, rt, m, m4, sv, mu, mu4 in Bc["live"]:
        assert toy.ball(q, side, ry, rt) == (m, m4, sv, mu, mu4); nb += 1
    return dict(A=sa, B=sb, replayed_A=na, replayed_B_balls=nb, checks_ok=True)


def unit_orbit(Z, x, R2):
    """m(x) = #{units e of O_H : |e x|^2 <= R2} / 2 (exact): Fincke-Pohst over e in Z^8 on the form
    |e alpha|^2 + |e beta|^2, then N(e) = +-1."""
    al, be = x[:8], x[8:]
    Ma = [Z["mul"](Z["y_pow"](i), al) for i in range(8)]; Mb = [Z["mul"](Z["y_pow"](i), be) for i in range(8)]
    G = [[sum(Ma[i][k] * Ma[j][k] + Mb[i][k] * Mb[j][k] for k in range(8)) for j in range(8)] for i in range(8)]
    L = [[0.0] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(i + 1):
            s = G[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            L[i][j] = math.sqrt(s) if i == j else s / L[j][j]
    e = [0] * 8; cnt = [0]
    def rec(i, rem):
        if i < 0:
            if any(e) and abs(Z["norm"](e)) == 1: cnt[0] += 1
            return
        c = sum(L[j][i] * e[j] for j in range(i + 1, 8)); r = math.sqrt(max(rem, 0.0)) / L[i][i]; ce = -c / L[i][i]
        for t in range(math.ceil(ce - r - 1e-9), math.floor(ce + r + 1e-9) + 1):
            d = (L[i][i] * t + c) ** 2
            if d <= rem + 1e-6:
                e[i] = t; rec(i - 1, rem - d)
        e[i] = 0
    rec(7, float(R2))
    return cnt[0] / 2


def verify_udup(Z, ev, udup_rel):
    """ev: runs [R, lattice, side, qbits, n, pr_charged, E[s], E[s/m]] and the 12 heaviest cells of run (48, 0)."""
    ys = [r[7] / r[5] for r in ev["runs"]]
    assert len(ys) == 48 and min(ys) >= udup_rel
    for x, s, m in ev["top200_R48_i0"][:6]:
        assert unit_orbit(Z, x, 48 * 48) == m / 2
    return dict(runs=len(ys), yr_min=round(min(ys), 4), yr_median=round(sorted(ys)[24], 4), yr_max=round(max(ys), 4), orbits_replayed=6)
