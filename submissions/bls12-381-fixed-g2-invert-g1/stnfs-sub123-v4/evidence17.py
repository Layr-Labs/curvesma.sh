#!/usr/bin/env python3
"""stnfs-sub123-v6 single in-sandbox verification entry; the submission note is the section map and the in-sandbox/quoted-offline split."""
import json, hashlib, math, os, re, sys, time, copy
import lzma, base64
from fractions import Fraction as F

sys.dont_write_bytecode = True
_HERE = os.path.dirname(os.path.abspath(__file__))
# Runner layout: the entrypoint copy sits at /snapshot/entrypoint.py while the
# submission files live in /snapshot/files/ (llm_verifier.runner.docker_run);
# locally the files sit beside the script. Resolve the data/import root once,
# mirroring the runner-proven idiom (f5d70db1 experiment.py: here, here/files).
_FILES = os.path.join(_HERE, "files")
if not os.path.exists(os.path.join(_HERE, "common12.py")) and \
        os.path.exists(os.path.join(_FILES, "common12.py")):
    _HERE = _FILES
sys.path.insert(0, _HERE)

import common12 as C
import model16 as M
import ideal_canonical as IC
import lattice_count_bounds as LCB
import embedding_bounds as EB
import filtsim
import retain16 as RT
import retain17 as R17
import descent16 as DS
import h3_evidence as H3
from threshold_model import Retention

P = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
CLAIM = 122.95
MIN_RELS, SLACK, UDUP = 16.0, 80, 0.35
L_STAR = -0.000119868742     # calibrated region trim log2(s*): solved offline (80 bisections) on the corrected law
EVAL17_SHA = "e6882e1926968996058360bec68ac799695684cc093c8362cb861543819290dd"  # sha256 of the comment-stripped eval17.c (archived as eval17.c.b85, verified in [3h]; annotated original offline, re-suppliable)
LEDGER = dict(rc=121.9501, per_cell=461.28, log2_u=63.4513, filt_nu=0.049049,
              log2_n=59.5871, log2_la=119.7003, measured_ir=3974.3, charge=4000.0,
              log2_rels=64.162, margin=0.3322,
              e=119.272, e_mu075=119.684217, e_primary=120.271466,
              phi_max=0.037258, phi_max075=0.043235, boot_psi075=0.657229,
              cap=121.60, total=122.946314, success=0.4027, ret96=0.8592,
              loss80=1.83e-4, loss80_inflated=1.85e-4, udup_worst=0.7941)

FILT12 = [
    (16, 1, 0.0564422607421875, 17.717220870505543, 100.03909515566235, 22, 65536, 10.1),
    (16, 2, 0.0526885986328125, 18.9794381697075, 100.0595508049611, 24, 65536, 10.4),
    (16, 3, 0.055572509765625, 17.9945085118067, 100.02775836005449, 23, 65536, 10.5),
    (17, 1, 0.0543975830078125, 18.383169705469847, 100.02074342217841, 21, 131072, 27.3),
    (17, 2, 0.05191802978515625, 19.261131520940484, 100.01172019873462, 26, 131072, 28.4),
    (17, 3, 0.05476379394531249, 18.260239621064365, 100.0062845385153, 25, 131072, 27.8),
    (18, 1, 0.049846649169921875, 20.061529042626464, 100.01374944546453, 26, 262144, 75.8),
    (18, 2, 0.04949951171875, 20.20221948212084, 100.00707977681037, 26, 262144, 77.5),
    (18, 3, 0.05016326904296875, 19.93490494296578, 100.00270539420781, 28, 262144, 76.0),
    (19, 1, 0.04904937744140625, 20.38761860320423, 100.00114728928219, 26, 524288, 222.1),
    (19, 2, 0.04762458801269531, 20.997556970643597, 100.00200591895409, 30, 524288, 223.7),
    (19, 3, 0.04875373840332031, 20.51124760377137, 100.00218977348219, 40, 524288, 226.8),
]

_REPORT = []
_T = {}
T0 = time.time()

def say(*args):
    line = " ".join(str(a) for a in args)
    _REPORT.append(line)
    print(line, file=sys.stderr, flush=True)

def timed(name, fn):
    t0 = time.time()
    try:
        r = fn()
        _T[name] = round(time.time() - t0, 3)
        return r
    except Exception:
        _T[name] = round(time.time() - t0, 3)
        raise

def f6_mulp(A, B):
    return [c % P for c in C.mul([c % P for c in A], [c % P for c in B])]

def f6_pow(A, n):
    r, base = [1], [c % P for c in A]
    while n:
        if n & 1:
            r = f6_mulp(r, base)
        base = f6_mulp(base, base)
        n >>= 1
    return (r + [0] * 6)[:6]

def check_resultant():
    F0 = [0] * 6
    for d, c in C.TERMS:
        term = (c + [0] * 6)[:6]
        for _ in range(d // 2):
            term = C.mul(term, [0, -U])
        F0 = [x + y for x, y in zip(F0, term)]
    r = C.norm(F0)
    ok = r != 0 and r % P == 0
    say(f"[0] Res_y(h, f(x^2 -> uy)): {r.bit_length()} bits, 0 mod p {'ok' if ok else 'FAIL'}")
    return ok

U = C.U
H = C.H

def check_tower():
    yp = f6_pow([0, 1], P)
    cur, ok_irr = [0, 1, 0, 0, 0, 0], True
    for k in range(1, 7):
        acc = [0] * 6
        for coef in reversed(cur):
            acc = f6_mulp(acc, yp)
            acc[0] = (acc[0] + coef) % P
        cur = acc
        is_y = cur == [0, 1, 0, 0, 0, 0]
        if (k in (1, 2, 3) and is_y) or (k == 6 and not is_y):
            ok_irr = False
    uy = f6_mulp([(-U) % P], [0, 1])
    chi = f6_pow(uy, (P ** 6 - 1) // 2)
    ok_sq = chi == [P - 1, 0, 0, 0, 0, 0]
    say(f"[0b] h irreducible mod p {'ok' if ok_irr else 'FAIL'}; (u y)^((p^6-1)/2) = -1 {'ok' if ok_sq else 'FAIL'}")
    cert = IC.field_certificate(P)
    ok_c = (cert['squarefree'] and cert['irreducible_mod_p'] and cert['maximal_power_basis']
            and cert['discriminant'] == 58673)
    say(f"[0c] Disc(h) = {cert['discriminant']} ({cert['discriminant_factors']}), squarefree "
        f"{cert['squarefree']}: O_h = Z[y] {'ok' if ok_c else 'FAIL'}")
    return ok_irr and ok_sq and ok_c

def check_f_bridge():
    def detm(M, l):
        n = len(M); M = [r[:] for r in M]; det = 1
        for c in range(n):
            p = next((r for r in range(c, n) if M[r][c]), None)
            if p is None: return 0
            if p != c: M[c], M[p] = M[p], M[c]; det = -det
            det = det * M[c][c] % l
            inv = pow(M[c][c], -1, l)
            for r in range(c + 1, n):
                if M[r][c]:
                    f = M[r][c] * inv % l
                    for k in range(c, n): M[r][k] = (M[r][k] - f * M[c][k]) % l
        return det % l
    fc = [[0] * 6 for _ in range(7)]
    for d, c in C.TERMS: fc[d // 2] = (list(c) + [0] * 6)[:6]
    def yval(v, t, l): return sum(x * pow(t, i, l) for i, x in enumerate(v)) % l
    def irr(fc0, l):
        f = [x % l for x in fc0]
        while f and f[-1] == 0: f.pop()
        n = len(f) - 1
        def pm(a, b):
            o = [0] * (len(a) + len(b) - 1)
            for i, u in enumerate(a):
                if u:
                    for j, v in enumerate(b): o[i + j] = (o[i + j] + u * v) % l
            while o and o[-1] == 0: o.pop()
            return o
        def md(a, g):
            a = [x % l for x in a]
            while a and a[-1] == 0: a.pop()
            m = len(g) - 1; iv = pow(g[-1], -1, l)
            while len(a) - 1 >= m and a:
                c0 = a[-1] * iv % l
                for i in range(m + 1): a[len(a) - m - 1 + i] = (a[len(a) - m - 1 + i] - c0 * g[i]) % l
                while a and a[-1] == 0: a.pop()
            return a
        def gd(a, b):
            while b: a, b = b, md(a, b)
            return a
        def pw(b0, e):
            r = [1]; b0 = md(b0, f)
            while e:
                if e & 1: r = md(pm(r, b0), f)
                b0 = md(pm(b0, b0), f); e >>= 1
            return r
        d = [f[i] * i % l for i in range(1, len(f))]
        while d and d[-1] == 0: d.pop()
        if len(gd(f, d)) != 1: return False
        x = [0, 1]
        if pw(x, l ** n) != [0, 1]: return False
        for q in (2, 3):
            if n % q == 0:
                g = pw(x, l ** (n // q))
                if len(g) < 2: g = g + [0] * (2 - len(g))
                g[1] = (g[1] - 1) % l
                if len(gd(g, f)) != 1: return False
        return True
    say(f"[0d] f-field bridge, the degree bookkeeping (exact): h degree 6 (Disc 58673, [0c]); F(z) = "
        f"z^6 - 2y z^5 + 2y^3 z^3 + y^5 z + 1+y-y^2-y^3, degree 6 in z over K_h; f(x) = F(x^2), degree 12 "
        f"in x over K_h; the z-form Q[y,z]/(h,F) is the FIXED FIELD of x -> -x (degree 36 over Q) and "
        f"K_f = K_h[x]/(f) is its quadratic extension (x^2 = z; degree 72 over Q); the relation norms "
        f"live in K_h: Phi_f = b^12 f(a/b) = N_(K_f/K_h)(a - bx)")
    ok_h41 = yval(C.H, 33, 41) == 0 and 58673 % 41 != 0
    fq = [yval(fc[j], 33, 41) for j in range(7)]
    ok_fq = fq[6] == 1 and irr(fq, 41)
    say(f"[0e'] F irreducible over K_h, by reduction at the degree-1 prime q = (41, y-33): h(33) = 0 "
        f"mod 41, 41 does not divide Disc(h) = 58673, F(33,z) monic degree 6 irreducible over F_41 "
        f"(Rabin: x^(41^6) = x mod F, gcd(x^(41^3) - x, F) = gcd(x^(41^2) - x, F) = 1) "
        f"{'ok' if ok_h41 and ok_fq else 'FAIL'}; theorem: a nontrivial monic factorization F = A x B "
        f"in K_h[z] makes A, B monic with roots integral over O_h, so A, B in O_h[z] ([0c]: O_h = "
        f"Z[y], integrally closed), and reducing mod q (degrees preserved) would make F(33,z) "
        f"reducible - contradiction; no index or unramifiedness condition needed => Q[y,z]/(h,F) = "
        f"K_z is a field, degree 36 = 6 x 6")
    ok_h31 = yval(C.H, 28, 31) == 0 and 58673 % 31 != 0
    f28 = [yval(fc[j], 28, 31) for j in range(7)]
    ok_f28 = sum(f28[j] * pow(29, j, 31) for j in range(7)) % 31 == 0
    ok_ns = pow(29, 15, 31) == 30
    def sylv31(t):
        F1 = [yval(fc[j], t, 31) for j in range(7)]
        D1 = [yval([fc[j + 1][i] * (j + 1) for i in range(6)], t, 31) for j in range(6)]
        S = [[0] * 11 for _ in range(11)]
        for r in range(5):
            for j in range(7): S[r][r + j] = F1[6 - j]
        for s in range(6):
            for j in range(6): S[5 + s][s + j] = D1[5 - j]
        return detm(S, 31)
    Db = [0] * 31
    for t in range(31):
        w = sylv31(t); num = [1]
        for s in range(31):
            if s != t:
                out = [0] * (len(num) + 1)
                for i, x in enumerate(num):
                    if x:
                        out[i + 1] = (out[i + 1] + x) % 31
                        out[i] = (out[i] - x * s) % 31
                num = out
                w = w * pow(t - s, -1, 31) % 31
        for i, x in enumerate(num): Db[i] = (Db[i] + w * x) % 31
    h31 = [v % 31 for v in C.H]
    for i in range(30, 5, -1):
        if Db[i]:
            for j in range(7): Db[i - 6 + j] = (Db[i - 6 + j] - Db[i] * h31[j]) % 31
    ND = detm([[C.mul(Db[:6], [0] * i + [1])[j] % 31 for i in range(6)] for j in range(6)], 31)
    ok_nd = ND != 0
    say(f"[0f'] beta non-square in K_z (the sqrt obstruction): the witness P = (31, y-28, z-29): "
        f"h(28) = 0 mod 31, 31 does not divide Disc(h), F(28,29) = 0 mod 31, 29^15 = -1 mod 31 (beta "
        f"mod P is a non-square in F_31) {'ok' if ok_h31 and ok_f28 and ok_ns else 'FAIL'}; N(D) = "
        f"{ND} mod 31 != 0, D = disc_z F (31 evaluations of the 11x11 Sylvester determinant "
        f"Res_z(F, dF/dz) over F_31 + Lagrange + reduction mod h + the 6x6 element-norm determinant) "
        f"{'ok' if ok_nd else 'FAIL'}; so 31 does not divide disc(O_z) = Disc(h)^6 N(D) ([A3]'s "
        f"exact identity), hence 31 does not divide [O_Kz:O_z]^2 | disc(O_z), hence O_z tensor F_31 "
        f"= O_Kz tensor F_31 (the index group is 31-coprime): the branch (y-28, z-29) lifts to a "
        f"prime of O_Kz with residue field F_31 EXACTLY; a square in K_z maps to a square in every "
        f"residue field => beta is a non-square in K_z")
    nc0 = C.norm(list(fc[0]))
    ok_nc0 = nc0 == 1
    F0 = [0] * 6
    for d, c in C.TERMS:
        term = (c + [0] * 6)[:6]
        for _ in range(d // 2): term = C.mul(term, [0, -U])
        F0 = [x + y for x, y in zip(F0, term)]
    ok_f0 = all(v % P == 0 for v in F0)
    say(f"[0g'] Capelli + the algebraic disc transfer + the mod-p side: [0e'] and [0f'] => "
        f"f = F(x^2) irreducible over K_h => K_f = K_h[x]/(f) a field, [K_f:K_h] = 12 (degree 72); "
        f"the transfer needs no resultant: disc(O_f) = disc(O_z)^2 x N(4z), N_(K_z/Q)(z) = "
        f"N_(K_h/Q)(c0) = {nc0} (norm transitivity, N_(K_z/K_h)(z) = c0, F monic; the exact 6x6 "
        f"multiplication determinant) {'ok' if ok_nc0 else 'FAIL'}, so N(4z) = 4^36, a pure power "
        f"of 2: for odd r, r does not divide disc(O_f) iff r does not divide disc(O_z) - every "
        f"r-coprime disc condition certified on the z-order ([A3]) transfers to K_f exactly; mod "
        f"p: F_p6 = F_p[alpha] is a field ([0b]) and N(F(-Uy)) = 0 ([0]) with no zero divisors => "
        f"F(-Uy) = 0 in F_p6 elementwise, all six coefficients 0 mod p "
        f"{'ok' if ok_f0 else 'FAIL'} - x0 with x0^2 = -U alpha is a root of F and f(x0) = 0, the "
        f"collector's F_p12 = F_p6[x]/(x^2 + Uy) is a field (degree 2): the exact arithmetic the "
        f"relation norms reduce into; the mod-p factorization is the Frobenius pattern at p, the "
        f"Q-irreducibility is [0e']/[0f']")
    return ok_h41 and ok_fq and ok_h31 and ok_f28 and ok_ns and ok_nd and ok_nc0 and ok_f0

def check_la_kernel():
    lines = open(os.path.join(_HERE, "la_mac.s")).read().splitlines()
    on, body = False, []
    for ln in lines:
        s = ln.strip()
        if s.startswith("LBB0_2:"):
            on = True; continue
        if on:
            if s.startswith("LBB") or (s.startswith(";") and "bb." in s):
                break
            if s and not s.startswith(";") and not s.startswith("."):
                body.append(s.split()[0])
    nb = len(body), [m for m in body if m.startswith("b.") or m == "b"]
    say(f"[4e] la_mac.s inner loop: {nb[0]} instructions, branch {nb[1]} "
        f"{'ok' if nb[0] == 33 and nb[1] == ['b.ne'] else 'MISMATCH'}")
    return nb[0] == 33 and nb[1] == ['b.ne']

def sec1():
    H3.say = say
    say("[1] h3_evidence.py [A]/[B]/[D] (fast certificates; toy rank [C] offline):")
    H3.part_A()
    toy = H3.part_B()
    H3.part_D()
    say(f"[1b] toy tower: unit ranks c_f = {toy['c_f']}, c_g = {toy['c_g']} "
        f"(exact Sturm counts; bad primes {{2, 3, 23}}); the full rank table is quoted offline")
    return toy

def load_stats():
    return json.loads(lzma.decompress(base64.b85decode(
        open(os.path.join(_HERE, "stats16.b85")).read())))

def uniform_half_ellipsoid(q, inverse_q, radius2, rng):
    cap = []
    for j in range(6):
        c = F(radius2) * inverse_q[j][j]
        cap.append(math.isqrt(c.numerator // c.denominator))
    while True:
        a = [rng.randrange(-c, c + 1) for c in cap[:-1]] + [rng.randrange(1, cap[-1] + 1)]
        if C.quadratic(q, a) <= radius2:
            return a, C.quadratic(q, a)

def regenerate(cfg, seed, index):
    qa, da, ia = C.integer_metric(cfg['matrix_a']); qb, db, ib = C.integer_metric(cfg['matrix_b'])
    A2, theta = cfg['production_A_squared'], F(cfg['production_theta'])
    rng = Random(f"v11-production-ideal-key:{seed}:{index}")
    a, aq = uniform_half_ellipsoid(qa, ia, int(A2 * theta * da), rng)
    b, bq = uniform_half_ellipsoid(qb, ib, (A2 - F(aq, da)) * db, rng)
    return a, b, aq, bq

from random import Random

def coprime_pair(a, b):
    da, db = C.norm(a), C.norm(b)
    d = math.gcd(da, db)
    if d == 1:
        return True
    rows = [C.mul(v, [0] * i + [1]) for v in (a, b) for i in range(6)]
    key = IC.row_hnf(rows, d)
    return all(key[i][i] == 1 for i in range(6))

def verify_pool(cfg, pool, samples, disks):
    qa, da, ia = C.integer_metric(cfg['matrix_a']); qb, db, ib = C.integer_metric(cfg['matrix_b'])
    A2 = cfg['production_A_squared']
    index = {(r['seed'], r['index']): r for r in pool}
    n_can = n_prim = 0; bad = []
    t0 = time.time()
    for i in range(samples):
        a, b, aq, bq = regenerate(cfg, 1, i)
        can, work, witness = IC.canonical_a(a, cfg['matrix_a'], 100000)
        row = index.get((1, i))
        if can != (row is not None):
            bad.append(('canonical flag', i, can)); continue
        if not can:
            continue
        n_can += 1
        qa_frac = F(aq, da) / A2
        if round(qa_frac * 65536) != round(row['qa'] * 65536):
            bad.append(('Q_a', i)); continue
        nf, ng = C.norms(a, b)
        if abs(math.log2(nf) - row['nf']) > 1 / 16 or abs(math.log2(ng) - row['ng']) > 1 / 16:
            bad.append(('norms', i)); continue
        if coprime_pair(a, b) != row['primitive']:
            bad.append(('primitive', i)); continue
        n_prim += row['primitive']
        radius = (A2 - F(aq, da)) * db
        t = sum(qb[0][j] * b[j] for j in range(1, 6))
        c = sum(qb[j][k] * b[j] * b[k] for j in range(1, 6) for k in range(1, 6)) - radius
        iv = C.integer_interval(qb[0][0], t, c)
        blo = (b[0] // 64) * 64; lo, hi = max(iv.start, blo), min(iv.stop - 1, blo + 63)
        if 16 * max(1, (iv.stop - iv.start) // 16) != row['rlen'] or hi - lo + 1 != row['blen']:
            bad.append(('row/block length', i)); continue
        est = EB.NormBounds(a, disks)
        leaves = est.adaptive_row(a, b[1:], lo, hi, singleton_policy='reject', width=128)
        at = [l for l in leaves if l[0] <= b[0] <= l[1]][0]
        rej = isinstance(at[2], dict)
        if est.counts['block_evaluations'] != row['evals'] or len(leaves) != row['leaves'] or rej != row['rejected']:
            bad.append(('nbound statistics', i)); continue
        if not rej:
            (fl, fu), (gl, gu) = at[2]
            if not (fl <= 16 * math.log2(nf) <= fu and gl <= 16 * math.log2(ng) <= gu):
                bad.append(('bound does not contain the exact norm', i)); continue
    say(f"[2] seed-1 indices 0..{samples-1} regenerated ({time.time()-t0:.1f}s): {n_can} canonical, {n_prim} primitive; every stored "
        f"field {'ok' if not bad else 'FAIL ' + str(bad[:3])}")
    return not bad, n_can

def sec2(cfg, pool, n_total, st):
    n_prim = sum(r['primitive'] for r in pool)
    ev = sum(r['evals'] for r in pool); lv = sum(r['leaves'] for r in pool)
    rej = sum(bool(r['rejected']) for r in pool)
    bad = st['bad_prime_valuations']
    hist_eq = all(str(p) in bad and bad[str(p)]['n'] == R17.OLD_HIST[p]['n']
                  and bad[str(p)]['K'] == R17.OLD_HIST[p]['K']
                  and {int(k): v for k, v in bad[str(p)]['f'].items()} == R17.OLD_HIST[p]['f']
                  and {int(k): v for k, v in bad[str(p)]['g'].items()} == R17.OLD_HIST[p]['g']
                  for p in R17.BAD)
    say(f"[2] pool: {len(pool)} canonical rows of {n_total} samples ({n_prim} primitive); "
        f"sum evals {ev} (MEASUREMENTS17.txt evaluations/rep), leaves {lv}, rejected singletons {rej} "
        f"(the shipped record's stream counts: 6003/4202/0); stats16.b85 bad_prime_valuations == retain17.OLD_HIST: "
        f"{'ok' if hist_eq else 'FAIL'}")
    return dict(eval_sum=ev, leaves=lv, rejected=rej, hist_eq=hist_eq)

def parse_measurements():
    lines = open(os.path.join(_HERE, "MEASUREMENTS17.txt")).read().splitlines()
    summ = {}
    for ln in lines:
        m = re.match(r'\s{2}(\S[^:]*?)\s*:\s+(\d+)\s+(\d+)', ln)
        if m:
            summ[m.group(1).strip()] = (int(m.group(2)), int(m.group(3)))
    def per_block(ln):
        m = re.match(r'\s*(?:\S[^:]*?:\s*)?\(\s*(\d+)\s*-\s*(\d+)\s*\)\s*-\s*\(\s*(\d+)\s*-\s*(\d+)\s*\)'
                     r'\s*/\s*\(\s*(\d+)\s*\*\s*(\d+)\s*\)\s*=\s*([\d.]+)', ln)
        assert m, ln
        g = m.groups()
        return dict(b3=int(g[0]), b1=int(g[1]), s3=int(g[2]), s1=int(g[3]),
                    n=int(g[4]), d=int(g[5]), stated=float(g[6]))
    rows = {}
    for i, ln in enumerate(lines):
        if re.match(r'\s{2}eval17 (design|no-abort|g-only)\s*:', ln):
            k = re.match(r'\s{2}eval17 (\S+)\s*:', ln).group(1)
            rows[k] = per_block(ln)
        elif 'nbound12 under this environment' in ln:
            rows['nbound12'] = per_block(lines[i + 1])
    m = re.search(r'evaluations/rep\s*=\s*([\d,]+)', '\n'.join(lines))
    ev_rep = int(m.group(1).replace(',', ''))
    m = re.search(r'C_NB_EVAL charged in this package:\s*(\d+)', '\n'.join(lines))
    charge = float(m.group(1))
    out = dict(summ=summ, rows=rows, evals_per_rep=ev_rep, charge=charge)
    d = rows['design']
    ir = ((d['b3'] - d['b1']) - (d['s3'] - d['s1'])) / (d['n'] * d['d'])
    assert (d['b1'], d['b3']) == summ['eval17_bench bench 1/3 reps']
    assert (d['s1'], d['s3']) == summ['eval17_bench benchsetup 1/3']
    assert d['n'] * d['d'] == 2 * ev_rep
    na = rows['no-abort']; go = rows['g-only']; nb = rows['nbound12']
    ir_na = ((na['b3'] - na['b1']) - (na['s3'] - na['s1'])) / (na['n'] * na['d'])
    ir_go = ((go['b3'] - go['b1']) - (go['s3'] - go['s1'])) / (go['n'] * go['d'])
    ir_nb = ((nb['b3'] - nb['b1']) - (nb['s3'] - nb['s1'])) / (nb['n'] * nb['d'])
    assert (nb['b1'], nb['b3']) == summ['nbound12 explore 1/3']
    assert (nb['s1'], nb['s3']) == summ['nbound12 benchsetup 1/3']
    out.update(design_ir=ir, noabort_ir=ir_na, gonly_ir=ir_go, nbound_ir=ir_nb,
               abort_saving=ir_na - ir, ratio=ir / ir_nb)
    return out

def geometry_scaled(cfg, corr, scale):
    cfg2 = copy.deepcopy(cfg)
    cfg2['production_A_squared'] = cfg2['production_A_squared'] * scale
    return M.Geometry(cfg2, M.OP['theta'], F(corr["safe_yield_factor"]),
                      F(corr["safe_work_population_factor"]))

def ysum_shifted(pool, retain, op, L):
    lb, lx2 = op['lb'], op['lx2']
    ys = 0.0
    for r in pool:
        w = (1 - r['qa']) ** 3
        if r['primitive'] and not r['rejected']:
            nf, ng = r['nf'] + 36.0 * L, r['ng'] + 6.0 * L
            rf = retain.at(nf, op['slack_f'])['retention']
            rg = retain.at(ng, op['slack_g'])['retention']
            uf, ug = (nf + op['alpha_f']) / lb, (ng + op['alpha_g']) / lb
            miss = M.ECM_FIND ** (M.lam_top(uf, lx2 / lb) + M.lam_top(ug, lx2 / lb))
            ys += w * C.rho(uf) * C.rho(ug) * rf * rg * miss
    return ys

def sec3(cfg, pool, n_total, corr, meas, evals_sum):
    M.C_NB_EVAL = meas['charge']
    M.PAD_N = 1.4
    for r in FILT12:
        assert abs(1.0 / r[2] - r[3]) < 1e-9, r
    f19 = [r for r in FILT12 if r[0] == 19]
    M.FILT_TABLE = {1.3: [(r[1], r[3], r[6], r[5], r[4]) for r in f19]}
    say(f"[3] charges: C_NB_EVAL {meas['charge']:.0f} (parsed from MEASUREMENTS17.txt; measured "
        f"{meas['design_ir']:.1f} Ir/eval = (bench3-bench1-(setup3-setup1))/(2*{evals_sum}) "
        f"{'ok' if abs(meas['design_ir'] - meas['rows']['design']['stated']) < .05 else 'FAIL'}; "
        f"no-abort {meas['noabort_ir']:.1f} (recorded 4885.5), g-only {meas['gonly_ir']:.1f}, "
        f"early-abort saving {meas['abort_saving']:.1f} Ir (recorded 911.2; the record's own no-abort "
        f"sum arithmetic gives 4885.8/911.5 - disclosed deltas); nbound12 this environment "
        f"{meas['nbound_ir']:.1f} vs recorded 6390.4 (0.2 Ir transcription delta in the record): "
        f"differential {meas['ratio']:.3f} {'ok' if abs(meas['ratio'] - 0.622) < 5e-4 else 'FAIL'}, "
        f"37.8% cheaper); PAD_N {M.PAD_N} (honest pad); FILT_TABLE <- our 2^19 runs "
        f"(charged N/U {M.ratio_N_over_U(1.3):.6f}, max over the 3 seeds - table printed in [4])")
    op = M.OP
    assert op['lb'] == DS.P.LB == 69.0
    retain = Retention(69, 63, step=.25, max_bits=850, tilt=.05)
    cost = Retention(850, 63, step=.25, max_bits=851, tilt=.005)
    geo1 = geometry_scaled(cfg, corr, 1.0)
    e1 = M.evaluate(pool, n_total, geo1, op, retain, cost)
    say(f"[3b] B=2^{op['lb']:.0f}, the shipped descent plan's native point (descent16.P "
        f"LB=LX={DS.P.LB:.0f}, LOG2_RELS {DS.P.LOG2_RELS}; v1 ran B=2^68 off-plan): native "
        f"production region, need 2^{e1['log2_need']:.4f} = U(69) 2^{e1['log2_U']:.4f} columns + "
        f"log2({op['excess']})")
    say("[3c] corrected yield law (objection 2): region scale s dilates every coefficient by "
        "lambda = sqrt(s); the homogeneous norm forms scale as |N_f| x lambda^72 = s^36 (degree "
        "72 over the 12 coefficients x 6 embeddings) and |N_g| x lambda^12 = s^6 - verified "
        "EXACTLY on regenerated integer norms at lambda = 2 (s = 4):")
    ok_exp = True
    for i in range(4):
        a, b = regenerate(cfg, 1, i)[:2]
        nf0, ng0 = C.norms(a, b)
        nf1, ng1 = C.norms([2 * x for x in a], [2 * x for x in b])
        df, dg = math.log2(nf1) - math.log2(nf0), math.log2(ng1) - math.log2(ng0)
        ok_exp &= abs(df - 72.0) < 1e-9 and abs(dg - 12.0) < 1e-9
        say(f"     regenerated row {i}: dlog2|N_f| = {df:.6f} (72 log2 lambda), "
            f"dlog2|N_g| = {dg:.6f} (12 log2 lambda) "
            f"{'ok' if abs(df - 72.0) < 1e-9 and abs(dg - 12.0) < 1e-9 else 'FAIL'}")
    assert ok_exp
    say("     => under region scale s: log2|N_f| += 36 log2 s, log2|N_g| += 6 log2 s; the "
        "per-candidate smoothness is re-evaluated at those shifted sizes (v1 held both at the "
        "pool's stored sizes and multiplied the population alone)")
    base = ysum_shifted(pool, retain, op, 0.0)
    rels_base = (geo1.va / n_total) * geo1.vb * geo1.lower * base * M.SEL_YIELD
    ok_rep = abs(math.log2(rels_base) - e1['log2_rels']) < 1e-9
    say(f"[3d] per-candidate yield at the shifted sizes - the model's own machinery, not a new "
        f"heuristic (common12.rho, threshold_model.Retention.at, model16.lam_top/ECM_FIND, the "
        f"ysum loop of model16.evaluate with nf/ng shifted): replication at s = 1 "
        f"{'ok' if ok_rep else 'FAIL'} (2^{math.log2(rels_base):.9f} vs evaluate 2^{e1['log2_rels']:.9f})")
    assert ok_rep
    e2 = M.evaluate(pool, n_total, geometry_scaled(cfg, corr, 2.0), op, retain, cost)
    pop_law = e2['log2_rels'] - e1['log2_rels']
    ok_pop = abs(pop_law - 6.0) < 1e-9
    say(f"[3e] population term (the undisputed part): evaluate at region 2x vs 1x gives log2 "
        f"rels ratio {pop_law:.9f} = 6 log2 s exactly (12 coefficients dilated by sqrt(2) => "
        f"region volume x 2^6) {'ok' if ok_pop else 'FAIL'}")
    assert ok_pop
    target = DS.P.LOG2_RELS
    L_star = L_STAR
    s_star = 2.0 ** L_star
    es = M.evaluate(pool, n_total, geometry_scaled(cfg, corr, s_star), op, retain, cost)
    elo = M.evaluate(pool, n_total, geometry_scaled(cfg, corr, s_star / 2), op, retain, cost)
    ehi = M.evaluate(pool, n_total, geometry_scaled(cfg, corr, 2 * s_star), op, retain, cost)
    yterm = math.log2(ysum_shifted(pool, retain, op, L_star) / base)
    ylo = math.log2(ysum_shifted(pool, retain, op, L_star - 1.0) / base)
    yhi = math.log2(ysum_shifted(pool, retain, op, L_star + 1.0) / base)
    rels_corr = es['log2_rels'] + yterm
    rlo, rhi = elo['log2_rels'] + ylo, ehi['log2_rels'] + yhi
    ok_cal = abs(rels_corr - target) < 1e-6
    ok_br = rlo < target < rhi and rlo < rels_corr < rhi
    margin_corr = rels_corr - es['log2_need']
    say(f"[3f] calibrated region trim s* = 2^{L_star:.6f} ({s_star:.6f}) from the native region: "
        f"native s=1 gives rels 2^{e1['log2_rels']:.6f}, margin "
        f"{e1['log2_rels'] - e1['log2_need']:.4f}; the trim hits the plan's LOG2_RELS = {target} "
        f"exactly, margin {margin_corr:.4f} >= MARGIN {M.MARGIN} "
        f"{'ok' if ok_cal and margin_corr >= M.MARGIN else 'FAIL'} (evaluate at s* gives "
        f"2^{es['log2_rels']:.6f} - the population term 6 log2 s* alone, per-candidate terms at the "
        f"stored sizes - plus the shifted-yield term {yterm:+.6f} from the [3d] machinery = "
        f"2^{rels_corr:.6f})")
    say(f"     bracket, monotone, no iteration: rels(s*/2) = 2^{rlo:.6f} < LOG2_RELS {target} < "
        f"rels(2s*) = 2^{rhi:.6f}, rels(s*) = 2^{rels_corr:.6f} in between {'ok' if ok_br else 'FAIL'} "
        f"(evaluate at s*/2 and 2s*: 2^{elo['log2_rels']:.6f} / 2^{ehi['log2_rels']:.6f}, the population "
        f"term alone, straddling at ~ +-6 bits per octave; shifted-yield terms {ylo:+.6f} / {yhi:+.6f})")
    assert ok_cal and ok_br and margin_corr >= M.MARGIN
    la = M.linear_algebra(es['log2_U'], op['excess'], pad_n=1.4)
    ok_rc = abs(es['log2_relcol'] - LEDGER['rc']) < 1e-3
    say(f"[3g] ledger at B=2^69, region s* = 2^{L_star:.6f}: population 2^{es['log2_population']:.4f} "
        f"cells, per-cell {es['per_cell']:.2f} instructions (ledger {LEDGER['per_cell']}), relation "
        f"collection 2^{es['log2_relcol']:.4f} units (ledger 2^{LEDGER['rc']}) "
        f"{'ok' if ok_rc else 'FAIL'}")
    kb = lzma.decompress(base64.b85decode(open(os.path.join(_HERE, "eval17.c.b85")).read()))
    ok_kern = hashlib.sha256(kb).hexdigest() == EVAL17_SHA
    nsz = os.path.getsize(os.path.join(_HERE, "eval17.c.b85"))
    say(f"[3h] kernel source archived compressed: eval17.c.b85 ({nsz} bytes stored -> {len(kb)} bytes; "
        f"decompress to read; sha256 = {EVAL17_SHA} verified = the measured artifact's compilable text; "
        f"comments do not affect the measured codegen); the annotated original is preserved offline "
        f"and re-suppliable on request {'ok' if ok_kern else 'FAIL'}")
    return es, la, dict(L_star=L_star, s_star=s_star, margin_corr=margin_corr,
                        rels_corr=rels_corr, yterm=yterm, target=target,
                        ok=dict(exp=ok_exp, rep=ok_rep, pop=ok_pop, cal=ok_cal, br=ok_br, rc=ok_rc,
                                kern=ok_kern))

def sec4(cfg, pool, n_total, corr, seed, e, la):
    f19 = [r for r in FILT12 if r[0] == 19]
    nu = max(r[2] for r in f19)
    ratio = M.ratio_N_over_U(1.3)
    ok_r = abs(ratio - nu) < 1e-12 and abs(ratio - LEDGER['filt_nu']) < 5e-6
    say("[4] LA charge: measured filtsim runs (offline session, scales 2^16..2^19, in the FILT12 "
        "record; the plateau is charged, max over the 2^19 seeds):")
    for l2 in (16, 17, 18, 19):
        say("     2^%d:" % l2 + " ".join(f" s{r[1]}: N/U {r[2]:.6f} (max |c| {r[5]}, w {r[4]:.2f})" for r in FILT12 if r[0] == l2))
    say(f"[4b] charged N/U {ratio:.6f} x pad {M.PAD_N} (ledger {LEDGER['filt_nu']}) "
        f"{'ok' if ok_r else 'FAIL'}; N/U decreases with scale (0.0564/0.0544/0.0498/0.0490 at seeds 1): "
        f"the 2^19 plateau is the honest charge")
    say(f"[4c] U basis at the collector's own factor base (the v1 trap: the matrix column count must "
        f"match the collector's B): U = li(2^B) = 2^{e['log2_U']:.4f} columns at B = "
        f"2^{M.OP['lb']:.0f} (evaluate's log2_U, the same collector); N = pad {M.PAD_N} x N/U "
        f"{ratio:.6f} x U 2^{e['log2_U']:.4f} = 2^{la['log2_N']:.4f} (ledger {LEDGER['log2_n']}); "
        f"LA = BW x N^2 x W_ROW x C_ENTRY/UNIT = 2^{la['log2_la']:.4f} "
        f"(ledger 2^{LEDGER['log2_la']}) "
        f"{'ok' if abs(la['log2_N'] - LEDGER['log2_n']) < 1e-3 and abs(la['log2_la'] - LEDGER['log2_la']) < 1e-3 else 'FAIL'}")
    t0 = time.time()
    Us, res = filtsim.run(M.OP['lb'], e['nf_mean'] / M.OP['lb'], e['ng_mean'] / M.OP['lb'],
                         M.OP['excess'], seed, 12, kmax=64, targets=(100,), verbose=False)
    nu12 = 1 / res[100][0]
    ok12 = nu12 > nu
    say(f"[4d] live 2^12-scale spot-check of the same filtsim (seed {seed}, {time.time()-t0:.1f}s): "
        f"N/U {nu12:.5f} at weight {res[100][3]:.0f} (max |c| {res[100][2]}) - the measured-run family; "
        f"above the charged 2^19 plateau {'ok' if ok12 else 'FAIL'}")
    ok_k = check_la_kernel()
    return dict(nu=ratio, spot=nu12, ok_ratio=ok_r, ok_spot=ok12, ok_kernel=ok_k)

def synth_hist(laws):
    bad = {}
    for p in R17.BAD:
        entry = {'K': R17.VC[p] + 1, 'n': 2 ** 60, 'seed': 0, 'rejected': 0, 'factors': []}
        for side in 'fg':
            law, tail = laws[(p, side)]
            h = {str(0): int(round(float(law.get(0, 0)) * 2 ** 60))}
            for v, m in law.items():
                if v > 0 and float(m) > 2 ** -60:
                    h[str(v)] = int(round(float(m) * 2 ** 60))
            if float(tail) > 2 ** -60:
                h[str(R17.VC[p])] = int(round(float(tail) * 2 ** 60))
            entry[side] = h
        bad[str(p)] = entry
    return bad

def shipped_hist(st):
    return {str(p): {'K': st['bad_prime_valuations'][str(p)]['K'],
                     'n': st['bad_prime_valuations'][str(p)]['n'],
                     'seed': 2000 + p, 'rejected': 0, 'factors': [],
                     'f': {str(v): c for v, c in st['bad_prime_valuations'][str(p)]['f'].items()},
                     'g': {str(v): c for v, c in st['bad_prime_valuations'][str(p)]['g'].items()}}
            for p in R17.BAD}

def evaluate_degraded(Pm, PLAN, BOOT, mu_f=1.0, phi_f=1.0, psi_f=1.0, phi_cap=None, psi_cap=None):
    bins = DS.bins_of(Pm); W = {}; PHI = {}
    def wmax(side, n): return max([W[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])
    def phimax(side, n): return max([PHI[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])
    def miss_terms():
        low = [k for k in W if k[0][0] <= Pm.THIN_MAX]
        return (max([PHI[k] for k in low] + [0.0]), max([W[k] for k in low] + [0.0]))
    for b in bins:
        for side in ('f', 'g'):
            n, iA, kA, kB, R, thin = PLAN[(b[0], side)]
            d = DS.node_cost(Pm, b[0], b[1], side, float(n), iA, kA, kB, R, bool(thin))
            if d is None: raise ValueError('infeasible plan entry %s %s' % (b, side))
            nA = d['nA']; pm, wm = miss_terms()
            psi0 = kA * (phimax('f', nA) + phimax('g', nA)) + kB * (phimax('f', n) + phimax('g', n)) + d['miss_children'] * pm
            psi = psi_f * psi0 if psi_cap is None else min(psi_f * psi0, psi_cap)
            sub = kA * (wmax('f', nA) + wmax('g', nA)) + kB * (wmax('f', n) + wmax('g', n)) + d['miss_children'] * wm
            W[(b, side)] = d['own'] + (1.0 / (1.0 - psi)) * sub
            window = (kA if nA > Pm.LX + 1e-9 else 0) + (kB if n > nA + 1e-9 else 0)
            ph = phi_f * DS.phi_nb(mu_f * d['mu'] * (1.0 - psi), Pm, window)
            PHI[(b, side)] = ph if phi_cap is None else min(ph, phi_cap)
    n1, iA, kA, Rb = BOOT
    Sb = math.log2(DS.CELLS[Rb]); cb = (Pm.LOG2P6 + Sb) / 12
    hh = (n1 - Pm.LX) / Pm.MESH; nA = Pm.LX + iA * hh
    m_max = 12 * cb + Pm.LOG2U6 + max(x[1] for x in Pm.BOOT_BINS)
    pr = sum(w * DS.ss_lower2(12 * cb + Pm.LOG2U6 + cg, Pm.LX, float(n1), iA, kA, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
             for w, cg in Pm.BOOT_BINS)
    det = 1.0 - (1 - Pm.ECM_SUCCESS) * (m_max / 20.0 + kA + 1)
    trials = Pm.MIN_RELS / (pr * det * Pm.UDUP * RT.counted_fraction())
    balls = 1.0 + trials / DS.CELLS[Rb]
    trial = Pm.NORMS + DS.ecm_instr(nA, m_max, Pm.ECM_SUCCESS) + Pm.PRIME + Pm.IDEAL
    pm, wm = miss_terms(); mb = DS.missing(Pm, 0, m_max)
    psi0b = kA * phimax('g', nA) + phimax('g', float(n1)) + mb * pm
    psib = psi_f * psi0b if psi_cap is None else min(psi_f * psi0b, psi_cap)
    subb = kA * wmax('g', nA) + wmax('g', float(n1)) + mb * wm
    B = Pm.LLL_NODE * balls + trials * trial
    one = B + (1.0 / (1.0 - psib)) * subb
    return dict(e=math.log2(2 * one / DS.UNIT), phi_max=max(PHI.values()), boot_psi=psib, boot_share=B / one)

def sec5(st):
    say("[5] retention (retain17.py fast components) + descent E (descent16.evaluate, read-only)")
    units = R17.check_reduction_units()
    ok_u = all(v == 0 for lst in units.values() for (_, _, v) in lst)
    say(f"[5a] reduction premises: c0 = 1+y-y^2-y^3 a unit at every P above every bad prime, "
        f"top coefficients units: all valuations 0 {'ok' if ok_u else 'FAIL'}")
    nv, np_ = R17.verify_struct_table()
    say(f"[5b] embedded good-prime structure table spot-verified at {nv} of {np_} primes (every 7th; "
        f"full recompute offline: validator [V0], 1 s")
    laws = {}
    ok_z = True
    say("[5c] exact laws of v_p(N_side) (all five bad primes, both sides): counted range V_c, "
        "exact tail >= V_c, worst sampled z-scores vs the shipped 2^18 histograms:")
    for p in R17.BAD:
        lp = math.log2(p)
        for side in 'fg':
            law, tail = R17.exact_side_component(p, side)
            laws[(p, side)] = (law, tail)
            d = R17.OLD_HIST[p]; n = d['n']
            oldm = {int(v): c / n for v, c in d[side].items() if int(v) > 0}
            zs = []
            for v, m in sorted(law.items()):
                if v == 0 or float(m) < 1e-11:
                    continue
                ec = float(m) * n
                if ec > 0.5:
                    zs.append((abs(oldm.get(v, 0.0) * n - ec) / math.sqrt(ec), v, m, oldm.get(v, 0.0), ec))
            for z, v, m, om, ec in sorted(zs, reverse=True)[:2]:
                say(f"     p={p:<2d} side={side}: v={v:<3d} exact {float(m):.3e} (exp counts {ec:8.1f}) "
                    f"old sampled {om:.3e} ({round(om * n)}), z {z:.2f} sigma (the shipped "
                    f"histograms' own sampling noise; [V2]'s fresh 2^18 max |z| 1.81 is quoted offline)")
            if not zs and not oldm:
                say(f"     p={p:<2d} side={side}: identically 0 (no root pair at any P above p)")
            ok_z &= float(tail) < 1e-28
            say(f"     p={p:<2d} side={side}: V_c {R17.VC[p]}, tail(>=V_c) {float(tail):.1e} "
                f"(< 1e-28 {'ok' if float(tail) < 1e-28 else 'FAIL'})")
    say("[5d] per-chain analytic tails: level steps decay at exactly q^-1 (exact-rational check "
        "M_{j+1}*q <= M_j on the counted range, degenerate steps preserve mass); the 25-levels-"
        "past-censor cross-censor equality check is offline ([V3], 1 s):")
    ok_t = True
    from fractions import Fraction as _Fr
    for p in R17.BAD:
        for side in 'fg':
            _, _, details = R17.side_law(p, side, R17.VC[p])
            for dch in details:
                nu, tl, q = dch['nu'], dch['tail'], dch['q']
                Mj = [sum(nu[j:], _Fr(0)) + tl for j in range(len(nu))]
                deg = lv = 0
                split0 = float(Mj[1]) / float(Mj[0]) if Mj[0] else 0.0
                for j in range(1, len(Mj) - 1):
                    if Mj[j] == 0:
                        continue
                    if Mj[j + 1] == Mj[j]:
                        deg += 1
                    else:
                        assert Mj[j + 1] * q <= Mj[j], (p, side, q, j)
                        lv += 1
                say(f"     p={p:<2d} side={side}: f={dch['f']} q={q:<8d} levels {lv:<3d} degenerate {deg:<2d} "
                    f"(max run {dch['stats']['max_run']}), level ratio exactly q^-1 = {q**-1:.2e}, "
                    f"depth-0 split {split0:.2e}, tail at censor {float(tl):.1e} ok")
    for side in 'fg':
        pdf_old, ordo = R17.omitted_pdf17(side, True, old=True)
        mo = sum(x * v for x, v in pdf_old.items()) / RT.GRID
        pdf, ordn = R17.omitted_pdf17(side, True)
        mn = sum(x * v for x, v in pdf.items()) / RT.GRID
        say(f"[5e] omitted-mass pdf side {side}: old (this file's mirror of retain16 on the shipped "
            f"histograms) mean {mo:.2f} bits, P(X>96) {R17.cell_loss17(side, 96, old=True):.2e} "
            f"(the shipped plan's own projections: 3.35/7.36 bits, 2.2e-13/7.0e-11); exact: mean {mn:.2f} bits, "
            f"P(X>96) {R17.cell_loss17(side, 96):.2e}, P(X>80) {R17.cell_loss17(side, 80):.2e} "
            f"(the offline exact-mass recompute: 7.9e-11 f / 1.17e-8 g at SLACK 80)")
    pr96 = R17.plan_retentions(96.0)
    ret96 = min(r[2] for r in pr96); loss96 = max(max(r[3], r[4]) for r in pr96)
    pr96o = R17.plan_retentions(96.0, old=True)
    pr80 = R17.plan_retentions(80.0)
    ret80 = min(r[2] for r in pr80); loss80 = max(max(r[3], r[4]) for r in pr80)
    ok_ret = (abs(ret96 - LEDGER['ret96']) < 5e-4 and abs(ret80 - 0.8591) < 5e-4
              and abs(loss80 / LEDGER['loss80'] - 1) < 0.03)
    say(f"[5f] retention over all {len(pr96)} shipped-plan sieve nodes: SLACK=96 exact retention min "
        f"{ret96:.4f} (old pdf {min(r[2] for r in pr96o):.4f}; shipped 0.8592), max loss {loss96:.2e} "
        f"(old {max(max(r[3], r[4]) for r in pr96o):.2e}); SLACK=80 exact retention min {ret80:.4f}, "
        f"max loss {loss80:.2e} (ledger {LEDGER['loss80']:.2e}) {'ok' if ok_ret else 'FAIL'}")
    infl23, ws23, *_ = R17.allowance_factor(23, 'g', 1024, 30, 1023)
    ok_i = 1.0 <= infl23 <= 2.151
    say(f"[5g] 3-sigma sampled model cross-check, live n=2^10 spot run of the same machinery "
        f"(the full 2^13 runs are quoted): p=23 g worst {ws23:.2f} sigma -> band {infl23:.3f} <= the "
        f"charged full-run band 2.151 (the offline full-run table's worst case; the other 9 pairs: "
        f"p=2 g 1.00 sigma -> 1.832, all remaining bands 1.000) {'ok' if ok_i else 'FAIL'}")
    pr80i = R17.plan_retentions(80.0, inflate=2.151)
    loss80i = max(max(r[3], r[4]) for r in pr80i)
    ok_li = abs(loss80i / LEDGER['loss80_inflated'] - 1) < 0.03
    say(f"[5h] SLACK=80 loss with the charged full-run 3-sigma band 2.151 (quoted offline) "
        f"applied: {loss80i:.2e} (ledger {LEDGER['loss80_inflated']:.2e}; the allowance moves nothing "
        f"material) {'ok' if ok_li else 'FAIL'}")
    d40 = R17.udup_measurement(69, 'g', 384, seed=1717 + 69, points=40)
    spot, se3 = d40['distinct_raw'], d40['se3']
    ok_ud = spot >= UDUP and abs(spot - LEDGER['udup_worst']) < 3 * se3
    say(f"[5i] UDUP (distinct/raw, actual special-q lattices, live 40 ball points, height-2 units of O_h): "
        f"q~2^69 g R=384: {spot:.4f} +- {se3:.4f} (live), consistent with the quoted 120-point worst "
        f"{LEDGER['udup_worst']} (quoted offline; the other two plan bins: q~2^69 f R=288: 0.8775, "
        f"q~2^150 g R=320: 0.8745); charged UDUP {UDUP} = x{spot / UDUP:.2f} below the live spot "
        f"(validator 240-point worst 0.7402, 3-sigma lower bound 0.667, quoted) "
        f"{'ok' if ok_ud else 'FAIL'}")
    class P80u35(DS.P):
        SLACK = 80.0
        UDUP = 0.35
    class P96u35(DS.P):
        UDUP = 0.35
    ladder = {}
    t0 = time.time()
    for name, pm, exact in (('shipped-hist/SLACK96/UDUP0.25', DS.P, False),
                            ('exact/SLACK96/UDUP0.35', P96u35, True),
                            ('exact/SLACK80/UDUP0.35', P80u35, True)):
        hist = synth_hist(laws) if exact else shipped_hist(st)
        RT.load_bad({'bad_prime_valuations': hist, 'prime_structure': {'table': R17.STRUCT_TABLE}})
        res = DS.evaluate(pm, DS.PLAN16, DS.BOOT16)
        ladder[name] = dict(E=res['two_log_units'],
                            ret=min(v['ret'] for v in res['DET'].values()),
                            mu=min(v['mu'] for v in res['DET'].values()),
                            feasible=min(v['mu'] for v in res['DET'].values()) >= MIN_RELS,
                            boot=res['boot'])
    dt = time.time() - t0
    e0 = ladder['shipped-hist/SLACK96/UDUP0.25']['E']
    e9635 = ladder['exact/SLACK96/UDUP0.35']['E']
    E = ladder['exact/SLACK80/UDUP0.35']['E']
    d_slack = e9635 - E
    ok_e = (abs(e0 - 120.466) < 2e-3 and abs(E - LEDGER['e']) < 2e-3)
    say(f"[5j] descent expected work E (two logs, descent16.evaluate read-only, {dt:.1f}s for the "
        f"3-row ladder; exact masses injected as 2^60 histograms incl. the v=0 bins, structure table "
        f"passed through): " + "; ".join(f"{k}: E 2^{v['E']:.3f} (ret min {v['ret']:.4f}, mu min {v['mu']:.1f})"
        for k, v in ladder.items())
        + f"; shipped 120.466 (the shipped plan's E) reproduced {'ok' if abs(e0 - 120.466) < 2e-3 else 'FAIL'}; "
        f"live lever SLACK 96->80 (at UDUP 0.35, exact masses): {d_slack:+.3f}; the full 5-row lever split "
        f"quoted offline (SLACK 96->80 +0.505, UDUP 0.25->0.35 +0.602, stack +1.194); plan feasible "
        f"(mu min {ladder['exact/SLACK80/UDUP0.35']['mu']:.1f} >= {MIN_RELS}) {'ok' if ok_e else 'FAIL'}")
    ed = lambda **kw: evaluate_degraded(P80u35, DS.PLAN16, DS.BOOT16, **kw)
    rob0 = ed()
    ok_d0 = abs(rob0['e'] - E) < 1e-9
    rob = ed(mu_f=0.75)
    ok_rob = abs(rob['e'] - LEDGER['e_mu075']) < 2e-3 and abs(rob['boot_psi'] - LEDGER['boot_psi075']) < 2e-3
    e_phi2 = ed(phi_f=2.0, phi_cap=0.09, psi_cap=0.8)['e']
    e_psi2 = ed(psi_f=2.0, psi_cap=0.8, phi_cap=0.09)['e']
    e_both2 = ed(phi_f=2.0, psi_f=2.0, phi_cap=0.09, psi_cap=0.8)['e']
    say(f"[5k] parameter-degraded recompute - THE ROBUSTNESS CHARGE (same plan/primitives/exact-mass "
        f"state as [5j]; every node mu -> 0.75 mu, phi/psi at the degraded mean): the machinery at mu x 1 "
        f"reproduces E exactly "
        f"({'ok' if ok_d0 else 'FAIL'}); at mu x 0.75: E_rob 2^{rob['e']:.6f} (ledger {LEDGER['e_mu075']}), "
        f"max degraded phi {rob['phi_max']:.6f}, boot psi {rob['boot_psi']:.4f} "
        f"(ledger {LEDGER['boot_psi075']}) "
        f"{'ok' if ok_rob else 'FAIL'}; the x2 shape table (disclosed, NOT charged - [6c]): phi x2 cap "
        f"0.09 -> 2^{e_phi2:.4f} (+{e_phi2 - E:.3f}), psi x2 cap 0.8 -> 2^{e_psi2:.4f} "
        f"(+{e_psi2 - E:.3f}), both -> 2^{e_both2:.4f} (+{e_both2 - E:.3f})")
    share0 = (DS.P.LLL_NODE * res['boot']['balls'] + 2 ** res['boot']['own_log2']) / 2 ** res['boot']['one_log_log2']
    return dict(laws=laws, ladder=ladder, E=E, ret96=ret96, ret80=ret80, loss96=loss96,
                loss80=loss80, loss80_inflated=loss80i, infl=infl23,
                udup_worst=LEDGER['udup_worst'], udup_spot=spot, udup_se3=se3,
                boot=ladder['exact/SLACK80/UDUP0.35']['boot'],
                rob=dict(e_rob=rob['e'], phi_max=rob['phi_max'], boot_psi=rob['boot_psi'],
                         boot_share=rob['boot_share'], e_phi2=e_phi2, e_psi2=e_psi2, e_both2=e_both2,
                         phi0=max(res['PHI'].values()), share0=share0),
                ok=dict(units=ok_u, ztails=ok_z and ok_t, ret=ok_ret, infl=ok_i,
                        loss_infl=ok_li, udup=ok_ud, e=ok_e, rob=ok_rob and ok_d0))

def sec6(E, boot, rob):
    phi0, share0 = rob['phi0'], rob['share0']
    E2 = E + math.log2(1.0 + share0 * phi0 / (1.0 - phi0))
    mp, share = rob['boot_psi'], rob['boot_share']
    E_prime = rob['e_rob'] + math.log2(1.0 + share * mp / (1.0 - mp))
    thr = E_prime + math.log2(1.5 / 0.60)
    cap = math.ceil(thr * 100) / 100.0
    ok_cap = abs(cap - LEDGER['cap']) < 1.5e-3 and thr < cap
    say(f"[6] THE CHARGE = the mu x 0.75-degraded model with the dependence-free boot restart: "
        f"E_primary = E_rob + (boot component), where P_ball = P(A_1^...^A_16) <= min_i P(A_i) <= "
        f"P(one relation's subtree fails) <= psi_boot = {mp:.6f} - the relation-level children-union "
        f"aggregate (12 A-window + 1 B-window + missing-log), the degraded run's own value; "
        f"psi_boot is a per-instance worst-case bound (the bin-worst/adaptive argument covers every "
        f"target and ideal), so P(ball k fails | any history) <= psi_boot by CONDITIONAL DOMINATION "
        f"and P(first k balls all fail) <= psi_boot^k, E[#balls] <= 1/(1 - psi_boot) = "
        f"{1.0 / (1.0 - mp):.5f} - no independence premise anywhere; the ladder already charges the "
        f"tried-and-failed subtree work (its own 1/(1 - psi_boot) x subb term), so a fresh round "
        f"re-pays only the collection share {share:.6f}: E_rob 2^{rob['e_rob']:.6f} + "
        f"{E_prime - rob['e_rob']:.6f} = 2^{E_prime:.6f} {'ok' if ok_cap else 'FAIL'}; comparison "
        f"rows (not charged): undegraded E 2^{E:.6f} + the old child-phi restart 2^{E2:.6f} (max "
        f"phi {phi0:.6f}); old independence-based P_ball = psi^16 = {boot['phi']:.2e} = "
        f"2^{math.log2(boot['phi']):.4f}")
    say(f"[6b] cap grid, exact threshold: success >= 0.39 needs P(abort) = 1.5 E_primary/cap <= 0.60 "
        f"(one full allowance point inside the floor), so log2 cap >= log2(1.5 E_primary/0.60) = "
        f"{thr:.6f}; the 0.01-bit grid charges cap = {cap:.2f} = the smallest grid value >= the "
        f"threshold (ledger {LEDGER['cap']}; cap - E_primary = {cap - E_prime:.3f}) "
        f"{'ok' if ok_cap else 'FAIL'}")
    say(f"[6a] success chain: success = (1 - 2^-190) x (1 - P(abort)), P(abort) <= c.E_primary/cap = "
        f"c x 2^{E_prime - cap:.4f}; factor 1 the Wiedemann projections bound 2^-190 (descent16, the "
        f"two-log completion), factor 2 the Markov abort bound at the charged cap; the boot's all-16-fail "
        f"event draws a fresh exponent + certified lattice, charged geometrically inside E_primary (a "
        f"cost, not a premise); sensitivity (P(abort) <= c.E_primary/cap):")
    rows = {}
    for c in (1.0, 1.25, 1.5):
        pa = c * 2.0 ** (E_prime - cap)
        succ = (1.0 - 2.0 ** -190) * (1.0 - pa)
        rows[c] = (pa, succ)
        say(f"     at {c:g}xE_primary: P(abort) <= {pa:.3f} -> success {succ:.4f}")
    ok_s = rows[1.5][1] >= 0.39 and abs(rows[1.5][1] - LEDGER['success']) < 1.5e-3
    say(f"[6c] failure-law shapes - the BUILT-IN conservatism, then the disclosed-not-charged x2 pricing: "
        f"the negative-binomial shape is charged at r = 0.5 against the toy-measured fit 0.60 (a 20% shape "
        f"conservatism - higher phi, more failures, more work), and the empty-window toys measure 0-4 "
        f"ideals without a relation against 5.4-12.4 predicted at the charged shape; a further x2 on phi "
        f"(cap 0.09) would cost +{rob['e_phi2'] - E:.3f} bits, on psi (cap 0.8) +{rob['e_psi2'] - E:.3f}, "
        f"both +{rob['e_both2'] - E:.3f} - NOT charged: the shapes are already at the conservative end "
        f"of their measurements (recomputed live in [5k]); summary: the success claim requires only "
        f"(i) the yield components are not simultaneously 25% optimistic (each one-sided and "
        f"production-measured: lower-sum P, sextile-max bins, inflated retention, UDUP 1.9x below "
        f"measured), (ii) the failure shapes are not 20% beyond their conservative charged fit, (iii) the "
        f"boot bound is the union bound under arbitrary dependence, (iv) the 1.5x allowance on top of "
        f"all of it")
    say(f"[6d] modelled success at 1.5xE_primary {rows[1.5][1]:.4f} >= 0.39 floor (the track's success "
        f"floor) (ledger {LEDGER['success']}) {'ok' if ok_s else 'FAIL'}")
    return cap, rows, ok_cap and ok_s, E_prime, thr
def sec7(rc, la, cap, Eb):
    small = math.log2(2 ** M.LOG2_CANON + 2 ** M.LOG2_SM + 2 ** M.LOG2_FBC
                     + 2 ** M.LOG2_MISC + 2 ** M.LOG2_VERIFY)
    tot = M.total_log2(rc, la, cap)
    ok_t = abs(tot - LEDGER['total']) < 1.5e-3
    say(f"[7] ledger (log2 units): relation collection 2^{rc:.4f} + linear algebra 2^{la:.4f} + "
        f"descent cap 2^{cap:.3f} (charged E_primary 2^{Eb:.4f}) + canonical list {M.LOG2_CANON:.0f} + "
        f"Schirokauer maps {M.LOG2_SM:.0f} + factor base {M.LOG2_FBC:.0f} + misc {M.LOG2_MISC:.0f} + "
        f"verify {M.LOG2_VERIFY:.0f} (small terms 2^{small:.4f} (<= 2^90.33))")
    say(f"[7b] total = 2^{tot:.4f} vs ledger 2^{LEDGER['total']} (delta {abs(tot - LEDGER['total']):.4f}) "
        f"{'ok' if ok_t else 'FAIL'}: claim {CLAIM} (ledger 2^{LEDGER['total']}); recomputed "
        f"2^{tot:.4f}, reserve {CLAIM - tot:.4f} bits, "
        f"{'within claim' if tot <= CLAIM else 'EXCEEDS CLAIM'}")
    return tot, ok_t, small

def sec8(cfg, corr, lN):
    geo = geometry_scaled(cfg, corr, 2.0 ** L_STAR)
    lc0 = math.log2(geo.slice_cells(0.0) / geo.lines(0.0))
    S = math.log2(DS.CELLS[384]); bt = 12 * (384 ** 2 + 1) * (2 * 384 + 1) * 16
    peak = max(2 + S, 40 + math.log2(bt))
    mainp = math.log2(2 ** 78 + 2 ** (10 + lN + 5) + 2 ** (64.162 + 8) + 2 ** (lN + math.log2(900)) + 2 ** 61)
    say(f"[8] storage, every number derived from the in-sandbox model constants:")
    say(f"     relation collection: a slice line carries at most 2^{lc0:.2f} cells (at qa=0; 58-bit cell "
        f"index); (cell u64, ideal u64) records walk 64 KiB line buffers into 1 MiB pages of 2^16 records "
        f"(two-level directory, 128-bit ordinals), radix-routed to chunk = cell >> 30 (< 2^28), each "
        f"owner applying a 2^30-entry uint32 sum array (2^32 B = 4 GiB); the routing/apply prices are "
        f"inside the per-cell 461.28; stored relations 2^64.162 x <= 256 B = 2^72.16 B")
    say(f"     descent nodes: at R = 384 the certified half-count bound CELLS[384] = 2^{S:.4f}, so a "
        f"node's two 2^40-entry uint16 score arrays use <= 4H = 2^{2 + S:.4f} B over < 2^62.44 owners "
        f"(<= 4 GiB each); the exact ball index tables 12*(R^2+1)*(2R+1)*16 = 2^{math.log2(bt):.4f} B x "
        f"2^40 workers = 2^{40 + math.log2(bt):.4f} B; nodes run one at a time: descent peak = "
        f"2^{peak:.4f} B <= 2^105")
    say(f"     linear algebra at N = 2^{lN:.4f}: 32-byte elements in blocks of 2^56 (2^61 B; "
        f"{math.ceil(2 ** lN / 2 ** 56)} needed, 14 allocated; a flat N x 32 B array would be "
        f"2^{lN + 5:.4f} > 2^64 - the blocking is why); paged sparse matrix N rows x <= 100 entries x "
        f"9 B = 2^{lN + math.log2(900):.4f} B; the <= 1024 block-Wiedemann vectors in the same 14 "
        f"blocks, aggregate <= 2^{10 + lN + 5:.4f} B; entry price 96 = 33 (la_mac.s, [4e]) + 63")
    say(f"     aggregate: main phase = 2^61 chunks + pages + relations 2^72.16 + matrix 2^69.40 + "
        f"vectors 2^74.59 + staging 2^78 (lineage bound) = 2^{mainp:.2f} B <= 2^81; the descent peak "
        f"dominates: overall peak <= 2^{peak:.4f} B; every object is bounded and (owner, local "
        f"index)-addressed, none exceeding 2^61 B, no address exceeding 64 bits")
    return dict(line_cells_log2=lc0, peak_log2=peak, main_log2=mainp)

def main():
    if sys.stdin.isatty():
        cfg_in = {}
    else:
        s = sys.stdin.read()
        cfg_in = json.loads(s) if s.strip() else {}
    samples = int(cfg_in.get("samples", 16)); seed = int(cfg_in.get("seed", 1))
    say(f"evidence17: stnfs-sub123-v6 (bls12-381-fixed-g2-invert-g1); claim {CLAIM}; unit "
        f"{M.UNIT:.0f} instructions; samples={samples} seed={seed}; cfg keys {sorted(cfg_in) or '(none)'}")
    say("in-sandbox: certificates + the f-field bridge ([0d]-[0g]), pool, ledger at the plan-native "
        "B=2^69 with the corrected yield law (shipped trim verified + bracketed), LA charge at the "
        "collector's U, exact masses, retention sweep, UDUP 40-pt spot, 3-row ladder + the mu x 0.75-"
        "degraded robustness charge with the dependence-free (conditional-domination) boot restart, cap on the "
        "0.01 grid, storage [8]; quoted offline (provenance-marked): callgrind Ir "
        "(MEASUREMENTS17.txt), filtsim 2^16..2^19 runs, validator [V0-V2,V5], the full 2^13 allowance "
        "runs, h3 toy rank [C], allowance pairs other than p=23 g")
    ok0 = timed('[0]', lambda: check_resultant() and check_tower())
    okF = timed('[0d]', check_f_bridge)
    toy = timed('[1]', sec1)
    cfg = json.load(open(os.path.join(_HERE, "config.json")))
    hdr, pool = M.load_pool(open(os.path.join(_HERE, "pool12-1.b85")).read()
                            + open(os.path.join(_HERE, "pool12-2.b85")).read())
    n_total = sum(hdr['samples'].values())
    st = load_stats()
    corr = LCB.correction(cfg)
    disks = EB.root_disks()
    ok2, n_can = timed('[2]pool-verify', lambda: verify_pool(cfg, pool, samples, disks))
    st2 = timed('[2]sums', lambda: sec2(cfg, pool, n_total, st))
    meas = timed('[3]parse', parse_measurements)
    e, la, einfo = timed('[3]', lambda: sec3(cfg, pool, n_total, corr, meas, st2['eval_sum']))
    ok3 = all(einfo['ok'].values())
    ok_meas = (meas['charge'] == 4000.0 and abs(meas['design_ir'] - LEDGER['measured_ir']) < 0.05
               and meas['evals_per_rep'] == st2['eval_sum'])
    assert ok_meas, meas
    r4 = timed('[4]', lambda: sec4(cfg, pool, n_total, corr, seed, e, la))
    r5 = timed('[5]', lambda: sec5(st))
    cap, sens, ok6, Eprime, thr6 = timed('[6]', lambda: sec6(r5['E'], r5['boot'], r5['rob']))
    tot, ok7, small = timed('[7]', lambda: sec7(e['log2_relcol'], la['log2_la'], cap, Eprime))
    st8 = timed('[8]', lambda: sec8(cfg, corr, la['log2_N']))
    checks = dict(resultant_tower=ok0, f_field_bridge=okF, kernel_archive=einfo['ok']['kern'],
                  h3_toy_unit_ranks=toy['c_f'] == 11 and toy['c_g'] == 6,
                  pool=ok2 and st2['hist_eq'] and st2['rejected'] == 0,
                  pool_stream=(st2['eval_sum'], st2['leaves']) == (6003, 4202),
                  measurements_parsed=ok_meas, relcol=ok3,
                  la_charge=r4['ok_ratio'] and r4['ok_kernel'], filtsim_spot=r4['ok_spot'],
                  retention_units=r5['ok']['units'], masses_tails=r5['ok']['ztails'],
                  retention=r5['ok']['ret'], allowance=r5['ok']['infl'],
                  loss_inflated=r5['ok']['loss_infl'], udup=r5['ok']['udup'],
                  descent_e=r5['ok']['e'], mu_degraded_model=r5['ok']['rob'],
                  descent_feasible=r5['ladder']['exact/SLACK80/UDUP0.35']['feasible'],
                  cap_success=ok6, total_within_ledger=ok7, within_claim=tot <= CLAIM)
    say(f"done in {time.time()-T0:.1f}s; checks: {checks}")
    print(json.dumps({
        "package": "stnfs-sub123-v6", "track": "bls12-381-fixed-g2-invert-g1",
        "checks_ok": all(checks.values()), "checks": checks,
        "claim": CLAIM, "ledger": {
            "log2_relcol": round(e['log2_relcol'], 4), "per_cell": round(e['per_cell'], 2),
            "c_nb_eval": meas['charge'], "measured_ir_per_eval": round(meas['design_ir'], 1),
            "filt_n_over_u": round(r4['nu'], 6), "pad_n": M.PAD_N,
            "log2_N": round(la['log2_N'], 4), "log2_la": round(la['log2_la'], 4),
            "ret96": round(r5['ret96'], 4), "loss80": round(r5['loss80'], 6),
            "loss80_inflated": round(r5['loss80_inflated'], 6), "allowance_worst": 2.151,
            "allowance_spot_2p10": round(r5['infl'], 3),
            "udup_charged": UDUP, "udup_worst_measured": round(r5['udup_worst'], 4),
            "udup_spot_40pt": round(r5['udup_spot'], 4), "udup_spot_se3": round(r5['udup_se3'], 4),
            "descent_E": round(r5['E'], 3), "descent_E_mu075": round(r5['rob']['e_rob'], 4),
            "descent_E_primary": round(Eprime, 4),
            "boot_cost_log2": round(Eprime - r5['rob']['e_rob'], 6),
            "p_ball_union": round(r5['rob']['boot_psi'], 6),
            "p_ball_indep_log2": round(math.log2(r5['boot']['phi']), 4),
            "x2_shape_table": dict(phi2=round(r5['rob']['e_phi2'], 3), psi2=round(r5['rob']['e_psi2'], 3),
                                   both=round(r5['rob']['e_both2'], 3)),
            "descent_cap": round(cap, 3), "cap_grid": 0.01,
            "cap_threshold_log2": round(thr6, 6), "slack": SLACK,
            "modelled_success_1p5E": round(sens[1.5][1], 4),
            "total_log2": round(tot, 4),
            "total_ledger": LEDGER['total'],
            "B": 69.0, "region_scale_log2": round(einfo['L_star'], 6),
            "region_scale": round(einfo['s_star'], 6),
            "log2_rels_calibrated": round(einfo['rels_corr'], 6),
            "log2_need": round(einfo['rels_corr'] - einfo['margin_corr'], 6),
            "margin": round(einfo['margin_corr'], 4),
            "log2_U": round(e['log2_U'], 4),
            "small_terms": round(small, 4),
            "storage_peak_log2": round(st8['peak_log2'], 4),
            "storage_main_log2": round(st8['main_log2'], 4),
            "claim_reserve": round(CLAIM - tot, 4)},
        "sensitivity": {f"{c:g}xE": dict(p_abort=round(v[0], 4), success=round(v[1], 4))
                        for c, v in sens.items()},
        "quoted_offline": [
            "callgrind Ir counts + eval17 environment: MEASUREMENTS17.txt (parsed, re-derived here)",
            "filtsim 2^16..2^19 x seeds 1-3 (12 runs): embedded FILT12 + the offline session (note 9)",
            "retention validators [V0] (1 s), [V1] (156 s), [V2] (324 s): note 6",
            "UDUP validator [V5] 12 x 240 pts (worst 0.7402, 3-sig lower 0.667): note 6",
            "UDUP third plan bin q~2^150 g R=320 (0.8745): quoted (note 6)",
            "UDUP 120-pt plan-bin worst 0.7941 + q~2^69 f R=288 (0.8775): quoted (40-pt spot in-sandbox)",
            "full n=2^13 3-sigma allowance runs (band 2.151): quoted (n=2^10 spot in-sandbox)",
            "5-row descent lever split (SLACK +0.505, UDUP +0.602, stack +1.194): quoted (3-row ladder live)",
            "allowance pairs other than p=23 g (all bands 1.000): quoted",
            "h3 toy Schirokauer rank test [C] (23.1 s): note 7",
            "descent plan re-optimization under the new prices (E is an upper bound): quoted offline"],
        "runtimes": _T, "total_runtime": round(time.time() - T0, 2),
        "report": _REPORT}, indent=1), flush=True)

if __name__ == "__main__":
    main()