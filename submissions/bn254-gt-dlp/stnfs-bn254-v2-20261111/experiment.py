#!/usr/bin/env python3
# Deterministic construction + finite-regime cost recomputation for the BN254 GT
# SexTNFS (Joux-Pierrot exTNFS) discrete-logarithm attack.
#
# ANALYTIC/VERIFICATION only. It (A) exhibits and VERIFIES the concrete F_{p^12}
# tower and the explicit polynomial pair (f,g,h) with integer coefficients, incl.
# irreducibility and the shared degree-12 factor mod p; (B) computes the two sieve
# norms by EXACT resultants on real random sieve elements at the production radius,
# and the smoothness probability from a self-checked Dickman rho; (C) builds the
# finite-regime ledger (relation collection, linear algebra with explicit matrix
# dimension/weight/Schirokauer columns, descent) and its sum; (D) states the claim
# = chain + disclosed margin and checks it is below the Pollard-rho baseline.
# It does NOT run, simulate, or prove a production-scale number field sieve.
# Standard library only. Deterministic: the CLAIM chain uses a fixed internal seed;
# the caller-supplied seed drives only an INDEPENDENT re-sample used as a cross-check.
import json, math, sys, random

LN2 = math.log(2.0)

# --------------------------------------------------------------------------- #
# BN254 public parameters (from the parameter x only).
# --------------------------------------------------------------------------- #
def bn254_params():
    x = 0x44e992b44a6909f1
    p = 36*x**4 + 36*x**3 + 24*x**2 + 6*x + 1
    r = 36*x**4 + 36*x**3 + 18*x**2 + 6*x + 1
    return x, p, r

# --------------------------------------------------------------------------- #
# (A) Field / polynomial construction, verified mod p.
# F_{p^12} = F_p[t,X]/(h(t), g(X)),  h = t^4+t+c (eta=4),  g = X^3 - x0 (kappa=3).
# Special side f(X) = P(X^3) = 36 X^12 + 36 X^9 + 24 X^6 + 6 X^3 + 1  (||f||=36).
# --------------------------------------------------------------------------- #
def pdeg(a):
    i = len(a) - 1
    while i > 0 and a[i] == 0:
        i -= 1
    return i

def prem(a, h, p):
    a = [c % p for c in a]; dh = pdeg(h); inv = pow(h[dh], p-2, p); da = pdeg(a)
    while da >= dh and not (da == 0 and a[0] == 0):
        c = (a[da]*inv) % p
        for k in range(dh+1):
            a[da-dh+k] = (a[da-dh+k] - c*h[k]) % p
        nd = pdeg(a)
        if nd == da:
            break
        da = nd
    return (a[:dh] + [0]*dh)[:dh]

def pmulmod(a, b, h, p):
    r = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i+j] = (r[i+j] + ai*bj) % p
    return prem(r, h, p)

def ppow(b, e, h, p):
    res = [1] + [0]*(pdeg(h)-1); b = prem(b[:] + [0]*pdeg(h), h, p)
    while e > 0:
        if e & 1:
            res = pmulmod(res, b, h, p)
        e >>= 1
        if e:
            b = pmulmod(b, b, h, p)
    return res

def pdivrem(a, b, p):
    a = [c % p for c in a]; db = pdeg(b); inv = pow(b[db], p-2, p); da = pdeg(a)
    while da >= db and not (da == 0 and a[0] == 0):
        c = (a[da]*inv) % p
        for k in range(db+1):
            a[da-db+k] = (a[da-db+k] - c*b[k]) % p
        nd = pdeg(a)
        if nd == da:
            break
        da = nd
    return a[:db] if db > 0 else [0]

def pgcd(a, b, p):
    a = [c % p for c in a]; b = [c % p for c in b]
    while not (pdeg(b) == 0 and b[0] == 0):
        a, b = b, pdivrem(a, b, p)
    return a

def irreducible_deg4(h, p):
    # h monic degree 4: irreducible iff X^(p^4)=X and gcd(X^(p^2)-X,h)=1.
    X = [0, 1, 0, 0]
    F4 = ppow(X, p**4, h, p)
    if (F4 + [0, 0, 0, 0])[:4] != X:
        return False
    F2 = ppow(X, p**2, h, p)
    d2 = [((F2 + [0, 0, 0, 0])[i] - (X + [0, 0, 0, 0])[i]) % p for i in range(4)]
    g = pgcd(h, d2, p)
    return pdeg(g) == 0 and g[0] != 0

def construct(p, x0):
    out = {}
    h = None
    for c in range(2, 80):
        cand = [c, 1, 0, 0, 1]
        if irreducible_deg4(cand, p):
            h = cand
            break
    out["h"] = "t^4 + t + %d" % (h[0] if h else -1)
    out["h_irreducible_mod_p_F_p4"] = h is not None
    out["p_eq_1_mod_3"] = ((p-1) % 3 == 0)
    out["x0_noncube_in_Fp"] = (pow(x0, (p-1)//3, p) != 1)
    valp4 = ppow([x0 % p], (p**4-1)//3, h, p) if h else [1, 0, 0, 0]
    out["x0_noncube_in_Fp4"] = ((valp4 + [0, 0, 0, 0])[:4] != [1, 0, 0, 0])
    out["g_X3_minus_x0_irreducible_over_Fp4"] = out["x0_noncube_in_Fp4"]
    out["P_of_x0_eq_0_mod_p"] = ((36*x0**4 + 36*x0**3 + 24*x0**2 + 6*x0 + 1) % p == 0)
    # f=P(X^3); mod (g=X^3-x0, p): X^3 == x0 so f == P(x0) == 0 => g | f mod p.
    out["f_eq_P_of_X3"] = "36X^12+36X^9+24X^6+6X^3+1"
    out["f_infnorm"] = 36
    out["f_divisible_by_g_mod_p"] = out["P_of_x0_eq_0_mod_p"]
    out["F_p12_realized_degree12"] = (out["h_irreducible_mod_p_F_p4"]
                                      and out["g_X3_minus_x0_irreducible_over_Fp4"])
    return h, out

# --------------------------------------------------------------------------- #
# Exact resultant norms over Z[t].  A sieve element a(t,X) has degree < tau in X
# and < eta in t, coefficients in [-A,A].  N_side = |Res_t(Res_X(a, side), h)|.
# --------------------------------------------------------------------------- #
def ztrim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a
def zmul(a, b):
    if a == [0] or b == [0]:
        return [0]
    r = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i+j] += ai*bj
    return ztrim(r)
def zsub(a, b):
    n = max(len(a), len(b)); r = [0]*n
    for i in range(len(a)):
        r[i] += a[i]
    for i in range(len(b)):
        r[i] -= b[i]
    return ztrim(r)
def zexactdiv(a, b):
    a = a[:]; db = len(b)-1; lead = b[-1]
    if a == [0]:
        return [0]
    q = [0]*(len(a)-db)
    for i in range(len(a)-1, db-1, -1):
        if a[i] == 0:
            continue
        c = a[i]//lead; q[i-db] = c
        for j in range(db+1):
            a[i-db+j] -= c*b[j]
    return ztrim(q)
ZERO = [0]; ONE = [1]

def bareiss_poly(M):
    n = len(M); M = [[M[i][j][:] for j in range(n)] for i in range(n)]
    sign = 1; prev = ONE
    for k in range(n-1):
        if M[k][k] == [0]:
            sw = None
            for i in range(k+1, n):
                if M[i][k] != [0]:
                    sw = i; break
            if sw is None:
                return ZERO
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k+1, n):
            for j in range(k+1, n):
                M[i][j] = zexactdiv(zsub(zmul(M[k][k], M[i][j]), zmul(M[i][k], M[k][j])), prev)
            M[i][k] = ZERO
        prev = M[k][k]
    d = M[n-1][n-1]
    return zmul(d, [-1]) if sign < 0 else d

def resX(apoly, fpoly):
    # apoly, fpoly: lists over X, each coeff a Z[t] poly. Returns Z[t] poly.
    if len(apoly) == 1:
        c = apoly[0]; r = ONE
        for _ in range(len(fpoly)-1):
            r = zmul(r, c)
        return r
    da = len(apoly)-1; df = len(fpoly)-1; n = da+df
    M = [[ZERO]*n for _ in range(n)]
    for i in range(df):
        for j in range(da+1):
            M[i][i+j] = apoly[da-j]
    for i in range(da):
        for j in range(df+1):
            M[df+i][i+j] = fpoly[df-j]
    return bareiss_poly(M)

def int_res(A, B):
    dA = len(A)-1; dB = len(B)-1; n = dA+dB
    if n == 0:
        return 1
    S = [[0]*n for _ in range(n)]
    for i in range(dB):
        for j in range(dA+1):
            S[i][i+j] = A[dA-j]
    for i in range(dA):
        for j in range(dB+1):
            S[dB+i][i+j] = B[dB-j]
    sign = 1; prev = 1
    for k in range(n-1):
        if S[k][k] == 0:
            sw = None
            for i in range(k+1, n):
                if S[i][k] != 0:
                    sw = i; break
            if sw is None:
                return 0
            S[k], S[sw] = S[sw], S[k]; sign = -sign
        for i in range(k+1, n):
            for j in range(k+1, n):
                S[i][j] = (S[k][k]*S[i][j] - S[i][k]*S[k][j])//prev
            S[i][k] = 0
        prev = S[k][k]
    return sign*S[n-1][n-1]

def norm(a_tx, side_x, h):
    return abs(int_res(resX(a_tx, side_x), h))

# --------------------------------------------------------------------------- #
# Dickman rho, stable Simpson integral recurrence u*rho(u)=int_{u-1}^u rho .
# With a self-consistency / convergence check (two step sizes agree; smooth ratio).
# --------------------------------------------------------------------------- #
def build_rho(STEP, UMAX=46.0):
    N = int(round(UMAX/STEP)); M = int(round(1.0/STEP)); rho = [1.0]*(N+1)
    for i in range(N+1):
        if i*STEP > 1.0:
            rho[i] = 0.0
    for i in range(M, N+1):
        u = i*STEP
        if u <= 1.0:
            continue
        s = rho[i-M] + rho[i]
        for j in range(1, M):
            s += (4.0 if j % 2 else 2.0)*rho[i-M+j]
        rho[i] = (s*STEP/3.0)/u
    return rho, STEP

def rho_at(tbl, STEP, u):
    if u <= 1.0:
        return 1.0
    if u >= 45.0:
        return 1e-300
    k = u/STEP; i = int(k); f = k-i
    return max(tbl[i]*(1-f) + tbl[i+1]*f, 1e-300)

# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
# Fixed design constants (chosen offline; experiment recomputes all derived work).
ETA, KAPPA, TAU = 4, 3, 3
A_RADIUS = 512          # sieve coefficient bound, 2^9
LOG2B = 62.4            # smoothness bound exponent (balances relation vs LA)
INTERNAL_SEED = 12345
INTERNAL_SAMPLES = 48
C_REF = 20000.0         # reference unit = one dense Fp12 mul (COST_MODEL.md)
MARGIN_BITS = 4.5       # disclosed aggregate margin (see proof.md section 7)
RHO_BASELINE = 135.216266

def run(seed=20261111, samples=48):
    x0, p, r = bn254_params()
    out = {}
    out["p_bits"] = p.bit_length(); out["r_bits"] = r.bit_length()
    Q = p**12; out["Q_bits"] = Q.bit_length()
    out["x0_bits"] = x0.bit_length()

    # (A) construction
    h, cinfo = construct(p, x0)
    out["construction"] = cinfo

    # Dickman rho + convergence self-check
    rho1, s1 = build_rho(0.01)
    rho2, s2 = build_rho(0.02)
    def R(u):
        return rho_at(rho1, s1, u)
    conv = max(abs(rho_at(rho1, s1, uu) - rho_at(rho2, s2, uu)) /
               max(rho_at(rho1, s1, uu), 1e-300) for uu in [3, 5, 7, 9, 11])
    out["dickman_rho_convergence_relerr"] = round(conv, 5)
    ratios = [R(u+1)/R(u) for u in range(3, 12)]
    out["dickman_ratio_monotone"] = all(ratios[i] <= ratios[i-1] + 1e-9 for i in range(1, len(ratios)))
    out["dickman_rho_samples"] = {str(u): R(u) for u in [4, 6, 8, 10]}

    # f,g as polynomials in X with Z[t] coefficients
    def kst(c):
        return [c]
    fpoly = [kst(1), ZERO, ZERO, kst(6), ZERO, ZERO, kst(24), ZERO, ZERO, kst(36), ZERO, ZERO, kst(36)]
    gpoly = [kst(-x0), ZERO, ZERO, kst(1)]

    def sample_norms(nsamp, sd):
        rnd = random.Random(sd); res = []
        for _ in range(nsamp):
            a = [ztrim([rnd.randint(-A_RADIUS, A_RADIUS) for _ in range(ETA)]) for _ in range(TAU)]
            Nf = norm(a, fpoly, h); Ng = norm(a, gpoly, h)
            if Nf > 1 and Ng > 1:
                res.append((Nf.bit_length(), Ng.bit_length()))
        return res

    # CLAIM chain from the fixed internal seed (deterministic).
    S = sample_norms(INTERNAL_SAMPLES, INTERNAL_SEED)
    nf = [a for a, b in S]; ng = [b for a, b in S]
    mf = sum(nf)/len(nf); mg = sum(ng)/len(ng)
    lnB = LOG2B*LN2
    uf = mf/LOG2B; ug = mg/LOG2B
    Psmooth = sum(R(a/LOG2B)*R(b/LOG2B) for a, b in S)/len(S)
    l2P = math.log2(Psmooth)
    out["sieve"] = {
        "tau": TAU, "eta": ETA, "kappa": KAPPA, "A_radius_log2": 9,
        "samples": len(S),
        "Nf_bits_mean": round(mf, 2), "Nf_bits_max": max(nf),
        "Ng_bits_mean": round(mg, 2), "Ng_bits_max": max(ng),
        "u_f": round(uf, 3), "u_g": round(ug, 3),
        "log2B": LOG2B, "P_smooth_pair_log2": round(l2P, 3),
    }

    # Factor base size and candidate count.
    log2FB = LOG2B + 1.0 - math.log2(lnB)          # ~ 2 * B / ln B primes (two sides)
    log2_candidates = log2FB - l2P                  # relations_needed / P(smooth pair)
    out["log2_factor_base"] = round(log2FB, 3)
    out["log2_candidates"] = round(log2_candidates, 3)

    # Per-candidate instruction count (bounded ABOVE), then conservative charge.
    sieve_updates = 2.0*math.log(lnB)               # Mertens: ~ln ln B reports per side
    instr_per_candidate = sieve_updates*4.0 + 8.0 + 50.0   # updates + init + cofactor bookkeeping
    out["instr_per_candidate_upper"] = round(instr_per_candidate, 1)
    rel_honest = log2_candidates + math.log2(instr_per_candidate/C_REF)
    rel_conservative = log2_candidates             # charge a FULL reference unit per candidate
    out["relation_collection_honest_log2"] = round(rel_honest, 3)
    out["relation_collection_claimed_log2"] = round(rel_conservative, 3)
    out["conservative_overcharge_bits"] = round(rel_conservative - rel_honest, 3)

    # Linear algebra: block Wiedemann, matrix dimension D = |FB|, derived weight,
    # Schirokauer columns, instruction-priced (NOT 1 unit/op).
    D = log2FB
    weight = (mf + mg)/(0.5*LOG2B)/2.0              # avg primes/relation ~ (logNf+logNg)/log(sqrt B)
    schir = ETA*12 + ETA*3                          # unit ranks ~ eta*deg_f + eta*deg_g dense cols
    per_rowop = weight*20.0 + schir*300.0           # weight small-int muladds + schir 254-bit modmuls
    LA = math.log2(3.0*(2.0**(2.0*D))*per_rowop/C_REF)
    out["linear_algebra"] = {"log2_D": round(D, 3), "row_weight": round(weight, 1),
                             "schirokauer_cols": schir, "instr_per_rowop": round(per_rowop, 0),
                             "log2": round(LA, 3)}

    # Individual-logarithm descent (both log_g and log_g z), capped Las Vegas.
    # Initial randomization reduces each target to an algebraic integer of norm
    # ~ L_Q[2/3, c0]; a bounded special-q tree drives primes down to the factor base.
    lnQ = out["Q_bits"]*LN2; lnlnQ = math.log(lnQ)
    def Lnorm(c23):
        return c23*(lnQ**(2.0/3.0))*(lnlnQ**(1.0/3.0))/LN2    # bits of L_Q[2/3,c23]
    N0 = Lnorm(0.95)
    out["descent_initial_norm_bits"] = round(N0, 1)
    c_cap = 40.0                                     # Las Vegas cap multiplier -> e^-40 per-node fail
    ladder = []
    nb = N0; nodes_log2 = 0.0; total = 0.0; fail = 0.0
    for lvl in range(6):
        u = nb/LOG2B; Pl = R(u); cap = c_cap/Pl
        cost = (2.0**nodes_log2)*cap                 # 1 unit per descent sieve candidate
        total += cost
        fail += (2.0**nodes_log2)*math.exp(-c_cap)
        ladder.append({"level": lvl, "norm_bits": round(nb, 1), "u": round(u, 2),
                       "P_level_log2": round(math.log2(Pl), 2), "nodes_log2": round(nodes_log2, 1),
                       "level_cost_log2": round(math.log2(cost), 2)})
        nb = max(nb*0.62, LOG2B + 2.0)               # geometric norm shrink per level
        nodes_log2 += 2.2                            # branching
    descent_one = math.log2(total)
    descent_both = descent_one + 1.0                 # both logs
    out["descent_levels"] = ladder
    out["descent_one_log_log2"] = round(descent_one, 2)
    out["descent_both_logs_log2"] = round(descent_both, 2)
    out["descent_total_nodes_log2"] = round(nodes_log2, 1)
    per_log_fail = fail                              # expected failed nodes bound
    success = 1.0 - 2.0*per_log_fail
    out["descent_per_log_fail_bound"] = per_log_fail
    out["success_probability"] = round(success, 6)

    # Ledger sum and claim.
    chain = math.log2(2.0**rel_conservative + 2.0**LA + 2.0**descent_both)
    out["chain_log2"] = round(chain, 4)
    out["margin_bits"] = MARGIN_BITS
    claim = chain + MARGIN_BITS
    # publish the claim as a fixed 6-decimal rational, rounded UP (never understated)
    claim6 = math.ceil(claim*1e6)/1e6
    out["claimed_log2W"] = claim6
    num = int(round(claim6*1e6))
    out["claimed_rational"] = "%d/1000000" % num

    # Storage (same counts).
    relbuf = log2FB + 8.0                            # ~|FB| relations kept, few hundred bytes each
    matrix = D + math.log2(weight) + 4.0            # sparse matrix, ~16 B/entry
    stor = math.log2(2.0**relbuf + 2.0**matrix)
    out["storage_log2_bytes"] = round(stor, 2)
    out["storage_bound_log2_bytes"] = 70

    # Independent cross-check with the caller's seed (claim must not depend on it).
    Sx = sample_norms(max(16, min(samples, 64)), seed)
    mfx = sum(a for a, b in Sx)/len(Sx); mgx = sum(b for a, b in Sx)/len(Sx)
    Px = sum(R(a/LOG2B)*R(b/LOG2B) for a, b in Sx)/len(Sx)
    candx = log2FB - math.log2(Px)
    out["crosscheck"] = {"seed": seed, "samples": len(Sx),
                         "Nf_bits_mean": round(mfx, 2), "Ng_bits_mean": round(mgx, 2),
                         "log2_candidates": round(candx, 3),
                         "within_1bit_of_claim_chain": abs(candx - log2_candidates) < 1.0}

    out["below_rho_baseline_135p216266"] = claim6 < RHO_BASELINE
    out["checks_ok"] = bool(
        cinfo["F_p12_realized_degree12"] and cinfo["f_divisible_by_g_mod_p"]
        and cinfo["h_irreducible_mod_p_F_p4"] and cinfo["g_X3_minus_x0_irreducible_over_Fp4"]
        and out["dickman_rho_convergence_relerr"] < 0.05 and out["dickman_ratio_monotone"]
        and rel_conservative > rel_honest
        and LA < rel_conservative and descent_both < rel_conservative
        and claim6 > chain and claim6 < RHO_BASELINE
        and out["success_probability"] >= 0.39
        and out["crosscheck"]["within_1bit_of_claim_chain"]
        and stor < out["storage_bound_log2_bytes"]
    )
    return out

def main():
    try:
        args = json.loads(sys.stdin.read() or "{}")
    except Exception:
        args = {}
    seed = max(0, min(int(args.get("seed", 20261111)), 4294967295))
    samples = max(16, min(int(args.get("samples", 48)), 256))
    print(json.dumps(run(seed=seed, samples=samples), indent=1))

if __name__ == "__main__":
    main()
