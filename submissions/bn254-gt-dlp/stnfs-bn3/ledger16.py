import base64, json, lzma, math, re
from fractions import Fraction as F
from math import log2, sqrt, floor, ceil, log
from pathlib import Path
import common12 as C, model16 as M, alpha12 as A, countmc as CM, lattice_count_bounds as LC, ideal_canonical as IC
U = C.U; P = 36 * U ** 4 + 36 * U ** 3 + 24 * U ** 2 + 6 * U + 1; R = 36 * U ** 4 + 36 * U ** 3 + 18 * U ** 2 + 6 * U + 1
PROFILE_P = "30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47"; PROFILE_R = "30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001"
HERE = Path(__file__).resolve().parent
DECL = json.loads((HERE / "declared.json").read_text())
GATE_TARGET = F(39, 100); P_CAP = 1e-9; NONLOG_CAP = 112.; WRAPPER_LOG2 = 68.; OUTPUT_LOG2 = DECL.get("output_instr_log2")
REC7_BASES, REC7_PREP = 4, 5.0e6      # reduced 7-dim bases per ideal (stretch classes, KERNELS.md) and instructions per basis
SIZING_COUNT = .9905          # count factor used by the sizing rule (fixed from a pilot count before the pools were sized)
SIZING_K = 2 * sqrt(2)
def disc_h():
    h = C.H; d = [i * c for i, c in enumerate(h)][1:]; n, m = 6, 5; S = [[0] * (n + m) for _ in range(n + m)]
    for i in range(m):
        for j, c in enumerate(reversed(h)): S[i][i + j] = c
    for i in range(n):
        for j, c in enumerate(reversed(d)): S[m + i][i + j] = c
    return abs(IC.determinant(S))      # |Res(h, h')| = |disc h| (h monic)
def maximal_order():
    """Dedekind's criterion for Z[y] at every prime whose square divides disc h."""
    good = True
    for p in (2, 5, 7):
        fac = []; rest = [c % p for c in C.H]; deg = 1; import itertools
        while len(rest) - 1 >= 2 * deg:
            for lo in itertools.product(range(p), repeat=deg):
                g = list(lo) + [1]
                while len(rest) > len(g) and not A.pmod(rest, g, p): rest = A.pdiv(rest, g, p); fac.append(g)
            deg += 1
        if len(rest) > 1: fac.append(rest)
        rad = [1]
        for g in {tuple(g) for g in fac}: rad = A._trim([sum(rad[i] * g[k - i] for i in range(len(rad)) if 0 <= k - i < len(g)) % p for k in range(len(rad) + len(g) - 1)])
        co = A.pdiv([c % p for c in C.H], rad, p); prod = [0] * (len(rad) + len(co) - 1)
        for i, x in enumerate(rad):
            for j, y in enumerate(co): prod[i + j] += x * y        # integer lift of rad x co
        t = [(a - b) // p for a, b in zip(prod + [0] * 7, C.H)]; assert all((a - b) % p == 0 for a, b in zip(prod + [0] * 7, C.H))
        good &= len(A.pgcd(A.pgcd(t, rad, p), co, p)) == 1 if len(co) > 1 else True
    return good
def archive():
    t = lzma.decompress(base64.b85decode("".join((HERE / "kernels_asm.b85").read_text().split()))).decode(); parts = re.split(r"(?m)^@@ (.*)\n", t); return {parts[i]: parts[i + 1] for i in range(1, len(parts), 2)}
def counts(pd, cfg, seed, replay=1 << 13):
    """Lattice points per volume: a (pool half-ellipsoid) and b (every slice radius of the region), each point at 6 SE; (lower, upper) factors."""
    A2 = cfg['production_A_squared']; lam = F(cfg['production_region_lam']); qa, da, ia = C.integer_metric(cfg['matrix_a']); qb0, db, _ = C.integer_metric(cfg['matrix_b']); qb = C.sieve_metric(qb0); ib = C.inverse(qb)
    ca, cb = pd['count']['a'], pd['count']['b']; Ra = int(A2 * F(cfg['production_theta']) * da); ra, sa = CM.ratio(qa, da, Ra, ca['cap'], ca['hits'], ca['n']); la, ua = ra * (1 - 6 * sa), ra * (1 + 6 * sa)
    tmax = int(F(cfg['production_theta']) * 4000); grid = [int((lam - F(t, 4000)) * A2 * db) for t in range(tmax, -1, -1)]; hits = [cb['h0']]
    for d in cb['dh']: hits.append(hits[-1] + d)
    assert len(hits) == len(grid); need = tmax - ceil(F(cfg['production_theta_region']) * 4000) - 1; rb = [CM.ratio(qb, db, g, cb['cap'], h, cb['n']) for g, h in zip(grid, hits)]
    lb = min(rb[j][0] * (1 - 6 * rb[j][1]) * (grid[j] / grid[j + 1]) ** 3 for j in range(need, len(grid) - 1)); ub = max(rb[j + 1][0] * (1 + 6 * rb[j + 1][1]) * (grid[j + 1] / grid[j]) ** 3 for j in range(need, len(grid) - 1))
    k = replay; rp = CM.replay(qa, ia, Ra, ca['seed'], 0, k, [Ra])[0] == ca['prefix'][1] and ca['prefix'][0] == k
    rp &= CM.replay(qb, ib, grid[-1], cb['seed'], 0, k, [grid[0], grid[-1]]) == cb['prefix'][1:] and cb['prefix'][0] == k
    return dict(a=(la, ua), b=(lb, ub), lower=la * lb, upper=max(1., ua) * ub, replay=rp, events=2 * (len(grid) - need) + 2, ra=ra)
def certified(cfg):
    """Fallback brackets by exact rational bounds (lattice_count_bounds.py), in bases reduced for the cell radius."""
    A2 = cfg['production_A_squared']; qa = [[F(x) for x in r] for r in cfg['matrix_a']]; qb = C.sieve_metric([[F(x) for x in r] for r in cfg['matrix_b']])
    B = [[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 2, -2, 0, 1, 0], [-1, 0, 1, 0, 0, 0], [2, -2, 0, 1, 0, 0], [-1, 0, 3, -2, 0, 1]]    # unimodular, a5 = last coordinate
    assert abs(IC.determinant(B)) == 1 and all(v[5] == 0 for v in B[:5]) and B[5][5] == 1
    qa2 = [[IC.inner(u, v, qa) for v in B] for u in B]; th = F(cfg['production_theta']); lam = F(cfg['production_region_lam'])
    a = LC.half_count_ratios(qa2, math.isqrt(int(th * A2))); b = LC.half_count_ratios(qb, math.isqrt(int((lam - th) * A2))); return float(a[0] * b[0]), float(a[1] * b[1]), float(a[0]), float(a[1])
def norm_sup(cfg):
    """Upper bounds on the region: power coordinates of a and b, log2 of N_f and N_g (per place: |x_e|^2 <= rho (u Q^-1 u + v Q^-1 v), x_e = sum x_j y_e^j = u.x + i v.x)."""
    import runs16
    A2 = cfg['production_A_squared']; out = {}; emb = {}
    for nm, key, rho in ('a', 'matrix_a', float(F(cfg['production_theta'])) * A2), ('b', 'matrix_b', float(F(cfg['production_region_lam'])) * A2):
        qi = [[float(x) for x in r] for r in C.inverse([[F(x) for x in r] for r in cfg[key]])]; out[nm] = max(sqrt(rho * qi[j][j]) for j in range(6)); emb[nm] = []
        for z in runs16._roots():
            u = [(z ** j).real for j in range(6)]; v = [(z ** j).imag for j in range(6)]; emb[nm].append(sqrt(rho * sum(u[i] * qi[i][j] * u[j] + v[i] * qi[i][j] * v[j] for i in range(6) for j in range(6))) * (1 + 1e-9))
    nf = ng = 0.
    for z, a, b in zip(runs16._roots(), emb['a'], emb['b']):
        w = abs(sum(c * z ** k for k, c in enumerate(C.W))) * (1 + 1e-9); nf += log2(sum(C.FC[k] * w ** k * a ** (8 - 2 * k) * b ** (2 * k) for k in range(5))); ng += log2(a * a + C.U * w * b * b)
    out.update(nf=nf, ng=ng); return out
def small_rows(st, geo, tJ):
    """log2 units of the fixed rows (allowances per object are stated in proof.md section 8)."""
    cfg = st['config']; lb = M.OP['lb']; li = lambda x: 2 ** M.log2_li(x); u = log2(M.UNIT)
    na = log2(geo.va * (float(F(cfg['production_theta_region']) / F(cfg['production_theta']))) ** 3)            # integer a with a5 >= 1, Q_a <= theta A^2 (count factor <= 1)
    pi_b = log2(1.25506 * 2 ** lb / (lb * M.LN2)); ret = M.retentions()[0]; rmin = min(ret.at(n / 4, s)['retention'] for n in range(4 * 300, 4 * 530) for s in (M.OP['slack_f'], M.OP['slack_f'] - 12))
    rels = tJ['log2_rels'] + 2 * tJ['yield_se_log2'] + log2(geo.upper / geo.lower) - log2(rmin)                # kept relations, upper
    bases = log2(li(lb) - li(M.P7)) + log2(REC7_BASES * REC7_PREP)
    L = ceil(lb - 6 + 0.53); pg = sum(2 ** (lb - 6 + 0.53 - j - 1) * (M.fft_F(2 ** (j + 1)) + 3 * 2 ** (j + 1)) for j in range(L))    # product tree of the primes <= B (2^(lb-6) x 1.443 limbs)
    return dict(canon=na + 42.8 - u, fb=pi_b + 23.4 - u, misc=rels + 6 + 24.4 - u, bases=bases - u, split=rels + log2(M.split_cost()) - u, gtree=log2(2 * pg) - u,
                cert=50., verify=11., wrapper=WRAPPER_LOG2 - u, rels_upper=rels, rmin=rmin, na=na, pi_b=pi_b)
def g2_output_bound():
    """Instructions of [x]Q0 on the twist over F_p2 (left-to-right double-and-add, Jacobian coordinates, mixed additions, one inversion):
    F_p multiplication <= 4096 and addition <= 256 instructions (loose four-limb Montgomery bounds), F_p2 product or square <= 3 mult + 5 add."""
    m2 = 3 * 4096 + 5 * 256; a2 = 2 * 256; dbl = 7 * m2 + 12 * a2; add = 11 * m2 + 12 * a2; inv = 2 ** 22 + 8 * m2
    return 254 * (dbl + add) + inv
def compute(st):
    """Every ledger number from the shipped statistics and module constants."""
    import la_bn as LA, descent_checks as DC
    cfg, pd = st['config'], st['pool_d']; op = M.OP; geo = M.Geometry(cfg, F(cfg['production_theta'])); ND = pd['draws_per_seed']; nS, nJ = ND * len(pd['pools']['S']), ND * len(pd['pools']['J'])
    tS, tJ = M.pool_totals(pd['sums']['S'], nS, geo), M.pool_totals(pd['sums']['J'], nJ, geo); lU = M.log2_fb(op['lb']); need = lU + log2(op['excess'])
    o = dict(tS=tS, tJ=tJ, log2_U=lU, need=need, margin=tJ['log2_rels'] - need, rels_2se=tJ['log2_rels'] + log2(1 - 2 * tJ['yield_se_log2'] * M.LN2), coll_est=tS['log2_est'], work_allow=1 + M.WORK_SE_K * tS['work_se'])
    rows = dict(collection=tS['log2_est'] + log2(M.CANON_WORK_ALLOWANCE)); la = LA.linear_algebra(LA.NU, LA.BW_M, LA.BW_N); rows.update(la=la['log2_la'], gen=la['log2_gen']); o['la'] = la
    t7, _ = DC.tree(); two, boot_fail = t7['two_log'], t7['boot_phi']; o.update(two_log=two, boot_fail=boot_fail, t7=t7); sr = small_rows(st, geo, tJ); o['small'] = sr
    rows.update(two_logs=two, canon=sr['canon'], sm=LA.schirokauer_cost(sr['rels_upper'])['log2_units'], misc=sr['misc'], fb=sr['fb'], bases=sr['bases'], split=sr['split'], gtree=sr['gtree'], cert=sr['cert'], verify=sr['verify'], wrapper=sr['wrapper'])
    if OUTPUT_LOG2: rows['output'] = log2(g2_output_bound()) - log2(M.UNIT)
    rows = {k: v if k == 'two_logs' else v + 1 for k, v in rows.items()}
    o['rows'] = rows; o['row_names'] = list(rows)
    return o
def totals(rows, cap, dn):
    """rows: declared (rounded-up) values. Expected run, every-run bound, success, gate, claim."""
    ex = log2(sum(2 ** float(v) for v in rows.values()))
    nonlog = sum(2 ** float(v) for k, v in rows.items() if k != 'two_logs')
    pnonlog = nonlog / 2 ** NONLOG_CAP
    rb = log2(1 + 2 ** NONLOG_CAP + 2 ** cap)
    pcap = 2 ** (float(rows['two_logs']) - cap)
    main_loss = 2 ** (dn['log2_N'] + 2 - 128) + P_CAP + 2 ** -48. + 2 / R
    succ = 1 - main_loss - pcap - pnonlog
    return dict(exp_run=ex, run_bound=rb, pcap=pcap, nonlog_abort=pnonlog, success=succ)
def gate(succ_decl, exp_run_decl):
    s = F(succ_decl); K = -(-GATE_TARGET * 2 ** 64 // s); assert K * s >= GATE_TARGET * 2 ** 64 > (K - 1) * s; return int(K), log2(1 + K / 2 ** 64 * 2 ** float(exp_run_decl))
def check_all(st, part, seed, say, ok):
    cfg, pd = st['config'], st['pool_d']; D = DECL; op = M.OP; out = {}; led = part in ('ledger', 'all'); geo = M.Geometry(cfg, F(cfg['production_theta']))
    c = counts(pd, cfg, seed); cert = certified(cfg); inside = cert[2] <= c['a'][0] and c['a'][1] <= cert[3] and cert[0] <= c['lower'] and c['upper'] <= cert[1]
    say(f"[3] points per volume: a [{c['a'][0]:.6f}, {c['a'][1]:.6f}] ({pd['count']['a']['n']:.3g} trials), b [{c['b'][0]:.6f}, {c['b'][1]:.6f}] on every slice radius ({pd['count']['b']['n']:.3g} trials, 6 SE each, {c['events']} one-sided events); "
        f"COUNT {M.COUNT} within ({c['lower']:.6f}, {c['upper']:.6f}); certified fallback ({cert[0]:.4f}, {cert[1]:.4f}) contains both; replays {ok('counts', M.COUNT[0] <= c['lower'] and c['upper'] <= M.COUNT[1] and inside and c['replay'])}")
    geo = M.Geometry(cfg, F(cfg['production_theta'])); ND = pd['draws_per_seed']; nS, nJ = ND * len(pd['pools']['S']), ND * len(pd['pools']['J'])
    al = pd['alpha']; sing = {}
    for p, (sf, sf2, sg, sg2) in al['raw'].items(): lo = lambda s, s2: s / al['n'] - 2 * sqrt(max(s2 / al['n'] - (s / al['n']) ** 2, 0) / al['n']); sing[int(p)] = (lo(sf, sf2), lo(sg, sg2))
    a2k = A.alpha(2000, sing); psum = al['partial']; good = all(abs(x - y) < 1e-9 for x, y in zip(a2k, psum['2000'])) and op['alpha_f'] >= psum['1000000'][0] + .05 and op['alpha_g'] >= psum['1000000'][1] + .05
    say(f"[14] alpha over odd primes (3, 5, 7 at pool means - 2 SE over {al['n']} primitive rows; Kummer-Dedekind root counts elsewhere): p < 2000 recomputed ({a2k[0]:.4f}, {a2k[1]:.4f}); stored partial sums "
        + "; ".join(f"{k}: ({v[0]:.3f}, {v[1]:.3f})" for k, v in psum.items()) + f"; charged ({op['alpha_f']}, {op['alpha_g']}) >= value at 10^6 + 0.05 {ok('alpha', good)}")
    vs = pd['vs']; dev = []
    for p, (sf, sf2, sg, sg2) in vs['sums'].items():
        if int(p) > 7:
            ef, eg = A.local(int(p)); n = vs['n']
            for m, s2, e in ((sf / n, sf2 / n, ef), (sg / n, sg2 / n, eg)): se = sqrt(max(s2 - m * m, 1e-12) / n); dev.append(abs(m - e) / max(se, 1e-4))
    h2 = al['v2hist']; k2 = max(M.kappa2(u) for u in (6.5, 8., 9.5)); good = set(h2) == {'0,0', '12,0', '0,2'} and max(dev) < 4 and max(M.KAP2) >= k2
    say(f"[14v] root-count law against exact valuations of {vs['n']} primitive pairs at the primes 11..47: largest deviation {max(dev):.2f} SE; exact 2-parts of all {al['n']} primitive rows: (v2 N_f, v2 N_g) in {{(0,0), (12,0), (0,2)}} {h2}; kappa2 <= {k2:.4f} <= {max(M.KAP2)} {ok('alpha', good)}")
    tS, tJ = M.pool_totals(pd['sums']['S'], nS, geo), M.pool_totals(pd['sums']['J'], nJ, geo); gs = M.Geometry(cfg, F(cfg['production_theta']), lower=SIZING_COUNT); t0, t1 = M.pool_totals(pd['sums']['S'], nS, gs), M.pool_totals(pd['sizing_prev'], nS, gs)
    lU = M.log2_fb(op['lb']); need = lU + log2(op['excess']); rule = lambda t: .06 + SIZING_K * t['yield_se_log2']; good = t0['log2_rels'] - need >= rule(t0) and t1['log2_rels'] - need < rule(t1)
    say(f"[4s] sizing on S: margin {t0['log2_rels'] - need:.6f} >= 0.06 + 2 sqrt(2) SE = {rule(t0):.6f} at theta = {cfg['production_theta_region']}; {t1['log2_rels'] - need:.6f} < {rule(t1):.6f} one step below: smallest {ok('sizing', good)}")
    mJ = tJ['log2_rels'] - need; se = tJ['yield_se_log2']; r2 = tJ['log2_rels'] + log2(1 - 2 * se * M.LN2); good = mJ >= max(M.MARGIN, M.MARGIN_SE * se) and r2 >= D['log2_rels_input'] and tJ['nf_rel'] <= D['filter_means'][0] and tJ['ng_rel'] <= D['filter_means'][1]
    say(f"[4] judged pool J (seeds {pd['pools']['J']}, {nJ} draws, {pd['sums']['J']['rows']} region rows): relations 2^{tJ['log2_rels']:.5f} against {op['excess']} U' = 2^{need:.5f} (U' = 2^{lU:.4f} tau-orbits): margin {floor(mJ * 1e4) / 1e4:.4f} >= max(0.06, 2 SE), SE {se:.5f}; "
        f"2-SE count 2^{floor(r2 * 1e4) / 1e4:.4f} >= descent input 2^{D['log2_rels_input']}; relation-weighted mean norms {tJ['nf_rel']:.2f} / {tJ['ng_rel']:.2f} <= filtering campaign {D['filter_means']} {ok('yield', good)}")
    alw = 1 + M.WORK_SE_K * tS['work_se']; coll = tS['log2_est'] + log2(M.CANON_WORK_ALLOWANCE); pc = pd['parts']
    say(f"[4a] work rows {pd['sums']['S']['wrows']} of S's {pd['sums']['S']['rows']} region rows: row SE {100 * tS['cond_se']:.3f}%, total SE {100 * tS['work_se']:.3f}%; allowance {M.CANON_WORK_ALLOWANCE} >= 1 + 2 SE = {alw:.5f} {ok('work', M.CANON_WORK_ALLOWANCE >= alw)}")
    say(f"[5] 2^{tS['log2_population']:.3f} cells at {tS['per_cell']:.2f} instructions (" + ", ".join(f"{k} {v:.2f}" for k, v in sorted(pc.items(), key=lambda x: -x[1])) + f"): collection 2^{tS['log2_est']:.6f} x {M.CANON_WORK_ALLOWANCE} = 2^{coll:.6f} units {ok('collection', coll <= float(D['rows']['collection']))}")
    rs = pd['run_stats']; cs = pd['contam']; bound = 1 - .05 ** (1 / cs['n'])
    say(f"[5r] run walk on runs of {M.RUN_T[0]}..{M.RUN_T[1]} a-values: {100 * rs['short']:.2f}% of the cells lie in shorter runs (slice walk), longest run {rs['T_max']}; contamination sample: {cs['n']} runs, {cs['points']} a-values, {cs['bad']} not canonical, "
        f"{cs['ends_bad']} runs end next to a canonical a (runs16.run is a sufficient test: that a starts its own run, as in the work rows): contaminated share <= {bound:.5f} <= RUN_ALLOW {M.RUN_ALLOW} (95%) {ok('runs', cs['bad'] == 0 and bound <= M.RUN_ALLOW and rs['T_max'] <= M.RUN_T[1] and rs['rejected'] == 0)}")
    out.update(collection=coll, margin=mJ, rels_2se=r2, tJ=tJ, tS=tS, need=need, log2_U=lU)
    lb = op['lb']; td, npr = M.tdiv_cost(); gb = M.gbatch_cost(M.pbits(0, lb)); s1 = M.gbatch_cost(M.pbits(0, M.BAND_TOP), c_leaf=2 ** 19); s2 = M.gbatch_cost(M.pbits(0, lb), c_leaf=2 ** 19); ar = archive()
    nk = max(int(x) for x in re.findall(r"max (\d+)", ar['norms.txt'])); sup = norm_sup(cfg); bad, fix = M.srt_toy(6, seed)
    fp = [x.split(":") for l in ar['fft.txt'].splitlines() if l.startswith("FFTPTS ") for x in l.split()[1:]]
    good = td <= M.C_TDIV and gb <= M.C_GBATCH and s1 <= M.C_BATCH and s2 <= M.C_BATCH2 and nk <= M.C_NORMS and ar['norms.txt'].count("mismatches 0") >= 1 and sup['a'] < 2 ** 12 and sup['b'] < 2 ** 13 and max(sup['nf'], sup['ng']) < M.NORM_BITS and pd['alpha']['maxbits'] <= M.NORM_BITS
    good &= bad == 0 and fix <= 2 and len(fp) >= 25 and all(float(c) <= M.fft_F(int(n)) for n, c in fp)
    say(f"[5t] g batch 2^{log2(gb):.3f} <= C_GBATCH 2^{log2(M.C_GBATCH):g}, strips 2^{log2(s1):.3f} / 2^{log2(s2):.3f} <= 2^{log2(M.C_BATCH):g} / 2^{log2(M.C_BATCH2):g}, trial division ({npr} primes) 2^{log2(td):.3f} <= 2^{log2(M.C_TDIV):g}; "
        f"norm kernel max {nk} <= C_NORMS {M.C_NORMS:.0f}; coordinates |a_j| <= {sup['a']:.0f} < 2^12, |b_j| <= {sup['b']:.0f} < 2^13; norms < 2^{sup['nf']:.1f}, 2^{sup['ng']:.1f} < 2^{M.NORM_BITS}; exact SRT on 2^6 planted leaves ({fix} corrections), {len(fp)} FFT points <= F {ok('batch', good)}")
    out.update(sup=sup)
    import kernstats as KS; k7 = KS.rec7_k(ar['rec7bn.txt']); nb = KS.nb_prices(ar['nb2bn.txt']); kmax = max(u for _, u, _ in k7)
    good = kmax <= min(k for _, k in M.REC7_K) and nb['kappa']['up'] <= M.C_NB_EVAL and nb['eta']['up'] <= M.C_NB_ROW and KS.prep_max(ar['prep7bn.txt']) <= REC7_PREP
    say(f"[10k] run walk at BN254: {len(k7)} strata, mean + 2 SE of Ir / carried price <= {kmax:.4f} <= REC7_K {M.REC7_K[0][1]}; norm-bound kernel at BN254: {nb['kappa']['up']:.1f} per evaluation <= C_NB_EVAL {M.C_NB_EVAL:.1f}, "
        f"{nb['eta']['up']:.1f} per row <= C_NB_ROW {M.C_NB_ROW:.1f}; basis preparation <= {KS.prep_max(ar['prep7bn.txt'])} <= {REC7_PREP:.0f} {ok('kernels', good)}")
    import la_bn as LA, descent_checks as DC
    if led:
        tab, fd = LA.filt_table(); fs = LA.filt_stats(tab); cr, ks = LA.charged_ratio(fs); k5 = fd['K5_RUNS']['16']; m5 = sum(k5) / len(k5); s5 = sqrt(sum((x - m5) ** 2 for x in k5) / (len(k5) - 1) / len(k5)); us = min(fd['UNIT_SHARE']['16'])
        good = LA.NU >= cr and max(v['maxc'] for v in fs.values()) <= 255 and m5 + 2 * s5 <= LA.C_ENTRY_K and us >= .86
        say("[6] filtering campaign (filtsim.py at the BN254 profile): " + "; ".join(f"2^{k}: n {v['n']} mean {v['mean']:.5f} SE {v['se']:.5f}" for k, v in fs.items()) + f"; NU = {LA.NU} >= max over 2^{ks[0]}, 2^{ks[1]} of mean + 2 SE = {cr:.6f}; "
            f"entries per row k5 {m5:.3f} + 2 SE <= {LA.C_ENTRY_K}, unit share >= {us} {ok('filter', good)}")
        la = LA.linear_algebra(LA.NU); kb = LA.kernel_bounds(); sg = LA.signatures(); un = LA.unramified()
        good = all(kb[k] for k in ('chi_sq', 'K_r_ge', 'D_lt', 'barrett_ok', 'row_mass_ok', 'mont256_ok', 'mont320_ok', 'kron_ok')) and kb['r_limbs'] == 4 and sg['sm'] == LA.N_SM and un['unramified'] and un['lead_ok']
        good &= la['log2_la'] <= float(D['rows']['la']) and la['log2_gen'] <= float(D['rows']['gen']) and [la['m'], la['n']] == D['bw']
        say(f"[6b] N = 2^{la['log2_N']:.4f} rows; block Wiedemann (m, n) = ({la['m']}, {la['n']}): {la['c_prod']:.2f} instructions per row-vector product, linear algebra 2^{la['log2_la']:.6f}, generator 2^{la['log2_gen']:.4f} units; "
            f"kernel bounds at r (4 limbs, Barrett slack {kb['barrett_slack']:.4f} < 2, chi {kb['chi_bits']} bits, chi^2 = -1); Schirokauer maps {sg['sm']} ({sg['sm_f']} on f, {sg['sm_g']} on g; signatures f {sg['f']}, g {sg['g']}); "
            f"r unramified in K_f, K_g (degree {un['deg_f']} / {un['deg_g']} absolute polynomials squarefree mod r) {ok('la', good)}")
        t7, res = DC.tree()
        good = t7['ok'] and t7['two_log'] <= float(D['rows']['two_logs'])
        good &= len(res['DET']) == 82
        good &= all(d['pool'] == {'f':12,'g':36}[side] and d['balls'] * d['det'] >= d['pool'] * (1-1e-12) and d['det'] >= .5 for (b,side),d in res['DET'].items())
        say(f"[7] randomized finite-family sieve: {t7['nodes']} node types, exact auxiliary-family sizes f12/g36, finite-schedule detection charged;  minimum modeled family-mean relations {t7['mu_min']:.6f}; two-log expectation 2^{t7['two_log']:.6f} {ok('descent', good)}")
        say('Conditional finite-family actual-mean positivity and boot smoothness are heuristic premises, not established by this arithmetic check. No negative-binomial, subtree-independence or DISC check is used.')
        o = compute(st); rows = o['rows']; dr = D['rows']; good = set(rows) == set(dr) and all(rows[k] <= float(dr[k]) for k in rows)
        tot = totals(dr, D['cap'], dict(log2_N=la['log2_N']))
        K, claim = gate(D['success'], D['exp_run'])
        good &= tot['success'] >= float(D['success']) and tot['exp_run'] <= float(D['exp_run']) and tot['run_bound'] <= float(D['run_bound']) and K == D['gate_K'] and claim <= float(D['claim'])
        good &= float(D['rows']['two_logs']) + 12 <= D['cap']
        say("[8] rows recomputed <= declared: " + ", ".join(f"{k} {rows[k]:.6f} <= {dr[k]}" for k in rows) + f"; expected run 2^{tot['exp_run']:.10f} <= 2^{D['exp_run']}, every run <= 2^{tot['run_bound']:.10f} (logs metered at 2^{D['cap']:g}); "
            f"success {tot['success']:.6f} >= {D['success']}; K = {K}; E[W] <= 1 + (K/2^64) 2^{D['exp_run']} = 2^{claim:.12f} <= 2^{D['claim']} {ok('ledger', good)}")
        out.update(rows=rows, success=tot['success'], claim=claim)
    if part in ('filter', 'all'):
        n = LA.replay(1)
        say(f"[6r] filtering simulator replay at 2^12 orbits (seed 1): N = {n} == {LA.REPLAY12} {ok('filter_replay', n == LA.REPLAY12)}")
    if part in ('descent', 'all'):
        rs = DC.residues(seed); po = DC.pools(seed); ss = DC.structure_spot()
        say(f"[11r] residue classes at 2, 3, 5, 7, 13 within the local bounds {rs[1]} {ok('residues', rs[0])}; "
            f"[11] fresh special-q and boot lattices inside the bins {po[1]} {ok('pools_d', po[0])}; local prime structure {ok('structure', ss)}")
    if part in ('laws', 'all'):
        lw = DC.laws(); say(f"[11b] omitted-mass laws recomputed from the local structure table (mean bits, P(X > 96)): {lw[1]} {ok('laws', lw[0])}")
    return out
