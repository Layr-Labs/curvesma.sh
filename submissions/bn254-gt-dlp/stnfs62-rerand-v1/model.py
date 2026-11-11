"""Cost and yield model (ledger) for the capped-descent (6,2) special TNFS on BN254 GT.

Python standard library only.  Every number of proof.md and of the manifest is produced by ledger() from the
sample file data.json and the constants below.
"""
import json
import math
from pathlib import Path

import lat

LN2 = math.log(2.0)
C_UNIT = 20000.0
EDGES = [56, 60, 66, 74, 86, 104, 130, 165]          # top edges (bits) of the special-q bins

# ---------------------------------------------------------------- design parameters (the design, not fitted data)
CFG = dict(
    log2_B=54,            # factor-base bound
    A=330,                # cell box: 12 coefficients in [-A, A]
    slack_bits=47,        # exact-threshold slack (< log2 B)
    surplus=2.0,          # distinct relations required by the yield model / matrix dimension
    retention=0.95,       # charged fraction of relations surviving the slack (prime powers); checked by the experiment
    lam_min=1024.0,       # required model-expected usable relations in the ball of the worst sampled ideal
    sieve_max_child=66,   # nodes with child bound <= this use the exact sieve, the others ECM on every cell
    ecm_k=3.0,            # curves = ecm_k / rho(u_ecm)
    ecm_u=6.0,            # ln(bound) / ln(B1)
    boot_edge=130,        # boot smoothness bound 2^boot_edge
    ecm_delta=0.8,        # charged probability that the ECM schedule splits one smooth norm (H-ECM)
    boot_safety=0.5,      # charged boot trial success = model x this   (H-boot)
    inv_mult_flat=0.75,   # charged mean 1/m where the sample has no unit multiplicities (measured 0.787, 0.790)
    eps_empty=1.0e-3,        # H-desc: probability that a visited node has no usable relation, given the trial so far
    leaf_edge=74,         # ideals up to 2^leaf_edge are leaves: their relation must lie entirely in KNOWN
    spread_shift_bits=36.0,   # lambda is also evaluated with both norms enlarged by half of this (ideal-to-ideal spread)
    lam_shifted_min=256.0,    # required lambda at the shifted norms
    mean_allowance=8.0,   # H-indep: a cell is a usable relation with probability >= model / this
    cluster_allowance=2.0,    # expected distinct (non-associate) relations >= lambda / this  (measured mean 1/m >= 0.75)
    tree_safety=2.0,      # charged mean descent work = tree_safety x model
    radii=(2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256),
    meter_factor=16.0,    # global meter on both logarithms = meter_factor x charged expectation
    main_fail=1.0e-4,     # allowance for all coin failures outside the two logarithms (Wiedemann, kernel vector, wrapper)
    gate_num=27, gate_den=64,
)

# ---------------------------------------------------------------- instruction prices (derived in KERNELS.md)
I_HIT = 32            # one sieve hit (match, two-word cell index, queue append, owner add)
I_RES = 32            # one tabulated or looked-up residue
I_CELL = 8            # counter initialisation and threshold scan per cell
I_NORM = 2 ** 19      # one exact norm (N_f or N_g) of a cell
I_IDEALTEST = 128     # one (relation, ideal) divisibility test
I_SM = 2.0 ** 38      # Schirokauer maps of one relation (both sides)
I_ROOTS = 2.0 ** 24   # root data of one prime
I_LLL = 2.0 ** 46     # one 12-dimensional LLL (boot or special-q lattice)
I_EXP = 2.0 ** 25     # one exponentiation in F_(p^12)
I_WRAPPER = 2.0 ** 70
NORM_SLACK_BITS = 8      # a cell (or boot vector) whose norm exceeds the largest sampled norm by more is skipped
REL_CAP = 2.0 ** 16      # usable relations kept per node (more are ignored)
LA_PER_ROW = 96 * 24 + 60 * 128 + 512 + 64 + 128     # one row of one matrix-vector product, incl. the dot product
LA_BM_PER_D2 = 544                                  # Berlekamp-Massey, per D^2
SM_COLS = 60


def norm_instr(bits):
    """One exact norm whose value has at most `bits` bits: 2^19 up to 512 bits (8-limb intermediates), scaled by
    the square of the limb count above that."""
    return I_NORM * max(1.0, (bits / 512.0) ** 2)


def modmul_instr(bits):
    """Montgomery product of two k-limb residues (k = ceil(bits/64)), schoolbook CIOS: charged 24 k^2 + 64 k."""
    k = math.ceil(bits / 64)
    return 24 * k * k + 64 * k


def gcd_instr(bits):
    """Binary gcd of two numbers of at most `bits` bits: at most 2*bits iterations of 8k + 8 instructions."""
    k = math.ceil(bits / 64)
    return 2 * bits * (8 * k + 8)


def trialdiv_instr(bits):
    """Division of a k-limb number by each of the 564 primes below 4096, with up to 2 exact divisions each."""
    k = math.ceil(bits / 64)
    return 564 * 3 * (4 * k + 8)


# ---------------------------------------------------------------- Dickman rho
class Dickman:
    """rho by u rho(u) = int_(u-1)^u rho(t) dt (trapezoid, step 1/512).  The window sum is recomputed at every
    step (no running difference), so small values keep full relative precision: rho(20) = 2.46e-29."""

    def __init__(self, step=512, umax=44):
        self.st = step
        n = umax * step + 1
        rho = [1.0] * n
        for i in range(step + 1, n):
            u = i / step
            rho[i] = (rho[i - step] / 2 + sum(rho[i - step + 1:i])) / step / (u - 0.5 / step)
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


PRIMES_1000 = [q for q in range(2, 1000) if all(q % d for d in range(2, int(q ** 0.5) + 1))]
C1000 = sum(math.log(q) / (q - 1) for q in PRIMES_1000)   # mean ln of the 1000-smooth part of a random integer


def mean_se(v):
    n = len(v)
    m = sum(v) / n
    var = sum((x - m) ** 2 for x in v) / (n - 1)
    return m, math.sqrt(var / n)


def l2(x):
    return math.log2(x) if x > 0 else float("-inf")


def ecm_schedule(bound_bits, cfg):
    """Stage-1-only schedule for prime parts up to 2^bound_bits: (log2 B1, curves, mulmods per curve, gcds per curve)."""
    D = _D()
    lb1 = bound_bits / cfg["ecm_u"]
    b1 = 2.0 ** lb1
    curves = math.ceil(cfg["ecm_k"] / D(cfg["ecm_u"]))
    mulmods = 20.0 * b1 + 11           # 11 products per ladder step, at most 1.811 B1 ladder steps (KERNELS.md)
    gcds = 1.25506 * b1 / (lb1 * LN2) + 1
    return lb1, curves, mulmods, gcds


def ecm_instr(bound_bits, modulus_bits, cfg):
    _, curves, mulmods, gcds = ecm_schedule(bound_bits, cfg)
    return curves * (mulmods * modmul_instr(modulus_bits) + gcds * gcd_instr(modulus_bits) + 2.0 ** 21)


_DICK = None


def _D():
    global _DICK
    if _DICK is None:
        _DICK = Dickman()
    return _DICK


def shift_bits(R):
    """Shift of (log2 N_f, log2 N_g) of ball cells when the radius changes from the sampled 16 to R (measured law)."""
    t = math.log2(R / 16.0)
    if t >= 0:
        return 48.25 * t, 12.25 * t
    return 44.0 * t, 10.8 * t


def load(base):
    return json.loads((Path(base) / "data.json").read_text())


SMALL_SUMS = (0.7533493177969591, 1.186961270401491)


def ledger(data, cfg=CFG, small_sums=SMALL_SUMS):
    """small_sums: sum of 1/q over the degree-one ideals above the primes 5 <= q < 4096, q != 197, per side
    (recomputed by the experiment, check E3)."""
    D = _D()
    b = cfg["log2_B"]
    lnB = b * LN2
    A = cfg["A"]
    out = {"config": dict(cfg)}
    # ------------------------------------------------------------------ main phase
    smp = data["main"]
    assert smp["A"] == A
    n_cells = smp["n"]
    inv_m = cfg["inv_mult_flat"]
    # admissible cells: (log2 N_f, log2 N_g, 1000-smooth bits of N_f, of N_g)
    adm = [(c[1] / 10, c[2] / 10, 1.0 / inv_m, c[3] / 10, c[4] / 10) for c in lat.unpack_cells(smp["adm"], "33333")]
    sf_m, sf_se = mean_se([r[3] for r in adm])
    sg_m, sg_se = mean_se([r[4] for r in adm])
    alpha_f = max(C1000 - (sf_m - 3 * sf_se) * LN2, 0.0)
    alpha_g = max(C1000 - (sg_m - 3 * sg_se) * LN2, 0.0)
    nf_min = sorted(r[0] for r in adm)[len(adm) // 200]      # 0.5% quantiles: smaller norms are not counted
    ng_min = sorted(r[1] for r in adm)[len(adm) // 200]
    ys = []
    for nf, ng, m, _, _ in adm:
        ok = nf >= nf_min and ng >= ng_min
        ys.append(D((nf * LN2 + alpha_f) / lnB) * D((ng * LN2 + alpha_g) / lnB) / m if ok else 0.0)
    ys += [0.0] * (n_cells - len(adm))
    y_m, y_se = mean_se(ys)
    S = 2 * A + 1
    V = S ** 12 / 2.0
    n_ideals = 2.2 * 2.0 ** b / lnB
    dim = n_ideals + SM_COLS + 64
    rel_low = V * (y_m - 3 * y_se) * cfg["retention"] * 0.98
    hits = sum(small_sums) + 2 * 1.05 * math.log(lnB / math.log(4096.0))
    nf_max = max(r[0] for r in adm)
    ng_max = max(r[1] for r in adm)
    stage2_cap = V * 2.0 ** -20                       # at most this many cells reach stage 2 (else failure)
    tf, tg = nf_min - cfg["slack_bits"] - 1, ng_min - cfg["slack_bits"] - 1
    pass_f = (nf_max + 8 - tf) * LN2 * D(tf / b)      # model bound: N has a B-smooth divisor >= 2^tf
    pass_g = (ng_max + 8 - tg) * LN2 * D(tg / b)
    stage2_model = V * pass_f * pass_g
    per_cell = I_HIT * hits + 3 * I_RES * n_ideals / S ** 6 + I_CELL
    collection = V * per_cell + stage2_cap * 2 * I_NORM
    rel_cap = 8 * dim                                 # relations factored over the factor base
    factor_rel = rel_cap * (I_IDEALTEST * n_ideals + 2 * I_NORM)
    schirokauer = dim * I_SM
    factor_base = 2.0 ** b / lnB * I_ROOTS
    la = (3 * LA_PER_ROW + LA_BM_PER_D2) * dim * dim
    out["main"] = dict(
        admissible_fraction=len(adm) / n_cells, alpha_f=alpha_f, alpha_g=alpha_g,
        mean_nf=sum(r[0] for r in adm) / len(adm), mean_ng=sum(r[1] for r in adm) / len(adm),
        max_nf=nf_max, max_ng=ng_max, charged_inverse_multiplicity=inv_m,
        log2_cells=l2(V), log2_yield=l2(y_m), yield_rel_se=y_se / y_m, log2_relations_lower=l2(rel_low),
        log2_ideals=l2(n_ideals), log2_dim=l2(dim), relations_over_dim=rel_low / dim,
        hits_per_cell=hits, instr_per_cell=per_cell, log2_stage2_cap=l2(stage2_cap), log2_stage2_model=l2(stage2_model),
        nf_min=nf_min, ng_min=ng_min,
        la_instr_per_D2=3 * LA_PER_ROW + LA_BM_PER_D2)
    u_f = (out["main"]["mean_nf"] * LN2 + alpha_f) / lnB
    u_g = (out["main"]["mean_ng"] * LN2 + alpha_g) / lnB

    def kappa(u, t):
        return D(u - t) / D(u) if u - t > 0 else 1.0 / D(u)

    def missing(norm_bits, bound_bits, side_u):
        """Expected number of prime ideals of norm in (2^(b-12), 2^b] of a bound-smooth norm with no logarithm."""
        uu = norm_bits / bound_bits
        tot, steps = 0.0, 240
        for i in range(steps):
            lo = b - 12 + 12 * i / steps
            hi = lo + 12 / steps
            if hi > bound_bits:
                break
            dens = kappa(uu, hi / bound_bits) * math.log(hi / lo)
            lam = dim * kappa(side_u, lo / b) / 2.0 ** hi
            tot += dens * math.exp(-lam)
        return tot

    def children(norm_bits, bound_bits, lo_bits, hi_bits):
        hi_bits = min(hi_bits, bound_bits)
        if hi_bits <= lo_bits or norm_bits <= 0:
            return 0.0
        uu = norm_bits / bound_bits
        tot, steps = 0.0, 100
        for i in range(steps):
            lo = lo_bits + (hi_bits - lo_bits) * i / steps
            hi = lo_bits + (hi_bits - lo_bits) * (i + 1) / steps
            tot += kappa(uu, hi / bound_bits) * math.log(hi / lo)
        return tot

    # ------------------------------------------------------------------ descent nodes
    def pois_cdf(k, mu):
        t, acc = math.exp(-mu), 0.0
        for i in range(k + 1):
            acc += t
            t *= mu / (i + 1)
        return acc

    bins = data["desc"]["bins"]
    per_ideal = data["desc"]["cells_per_ideal"]
    W = {}            # model mean work of a node of the bin with its whole subtree (worse side)
    NODES = {}        # model mean number of nodes of the subtree, the node included (worse side)
    rows = []
    for idx, n in enumerate(EDGES):
        best_side = {}
        for side in "fg":
            samp = [lat.unpack_cells(it, "233") for it in bins[f"{side}{n}"]]
            cand = None
            leaf = n <= cfg["leaf_edge"]
            # leaf bins: the relation must lie entirely in KNOWN (child bound B, no children);
            # big bins: children are leaves (child bound at most 2^leaf_edge)
            for c in ([b] if leaf else [e for e in EDGES[:idx] if e <= cfg["leaf_edge"]]):
                for R in cfg["radii"]:
                    df, dg = shift_bits(R)
                    ps, ps_sh, mf, mg = [], [], 0.0, 0.0
                    sh = cfg["spread_shift_bits"] / 2.0    # each norm enlarged by this many bits (spread study)
                    for cl in samp:
                        acc, acc_sh = 0.0, 0.0
                        for _, nf10, ng10 in cl:
                            nf, ng = nf10 / 10 + df, ng10 / 10 + dg
                            if side == "f":
                                nf -= n - 1
                            else:
                                ng -= n - 1
                            nf, ng = max(nf, 1.0), max(ng, 1.0)
                            mf, mg = max(mf, nf), max(mg, ng)
                            acc += D((nf * LN2 + alpha_f) / (c * LN2)) * D((ng * LN2 + alpha_g) / (c * LN2))
                            acc_sh += D(((nf + sh) * LN2 + alpha_f) / (c * LN2)) * D(((ng + sh) * LN2 + alpha_g) / (c * LN2))
                        ps.append(acc / per_ideal)
                        ps_sh.append(acc_sh / per_ideal)
                    M = (2 * R + 1) ** 12 / 2.0
                    miss = missing(mf, c, u_f) + missing(mg, c, u_g)
                    if miss >= 0.5:
                        continue
                    ecm_mode = c > cfg["sieve_max_child"]
                    lam = M * min(ps) * (1 - miss) * (cfg["ecm_delta"] ** 2 if ecm_mode else cfg["retention"])
                    lam_sh = lam * min(ps_sh) / min(ps)
                    if lam < cfg["lam_min"] or lam_sh < cfg["lam_shifted_min"]:
                        continue
                    if ecm_mode:
                        cost = M * (norm_instr(mf + 8 + n) + norm_instr(mg + 8 + n) + trialdiv_instr(mf + 8) + trialdiv_instr(mg + 8)
                                    + ecm_instr(c, mf + 8, cfg) + ecm_instr(c, mg + 8, cfg)) + REL_CAP * (2.0 ** 34 + I_SM)
                    else:
                        n_c = 2.2 * 2.0 ** c / (c * LN2)
                        hits_c = sum(small_sums) + 2 * 1.05 * math.log(c * LN2 / math.log(4096.0))
                        cost = (3 * I_RES * (2 * R + 1) ** 6 * n_c / 2 + M * (I_HIT * hits_c + I_CELL)
                                + M * 2.0 ** -20 * (norm_instr(mf + 8 + n) + norm_instr(mg + 8 + n))
                                + REL_CAP * (I_IDEALTEST * n_c + norm_instr(mf + 8 + n) + norm_instr(mg + 8 + n) + I_SM))
                    cost += I_LLL
                    kids, cw, sub, lo = 0.0, 0.0, 0.0, b
                    for e in EDGES:
                        if lo >= c:
                            break
                        k = children(mf, c, lo, e) + children(mg, c, lo, e)
                        kids += k
                        cw += k * W[e]
                        sub += k * NODES[e]
                        lo = e
                    w = cost + cw
                    score = w
                    if cand is None or score < cand[0]:
                        cand = (score, c, R, "ecm" if ecm_mode else "sieve", l2(min(ps)), l2(sum(ps) / len(ps)), l2(lam),
                                kids, miss, mf, mg, cost, 1 + sub, w, l2(lam_sh),
                                0 if leaf else int((mf + mg + 2 * NORM_SLACK_BITS) // b))
            if cand is None:
                raise ValueError(f"no admissible plan for bin {side}{n}")
            best_side[side] = cand
            rows.append(dict(side=side, top_bits=n, child_bound_bits=cand[1], radius=cand[2], mode=cand[3],
                             log2_p_worst=cand[4], log2_p_mean=cand[5], log2_lambda_worst=cand[6], log2_lambda_shifted=cand[14], mean_children=cand[7],
                             missing=cand[8], max_nf=cand[9], max_ng=cand[10], log2_node_instr=l2(cand[11]),
                             log2_subtree_instr=l2(cand[13]), subtree_nodes=cand[12], max_children=cand[15]))
        W[n] = max(best_side["f"][13], best_side["g"][13])
        NODES[n] = max(best_side["f"][12], best_side["g"][12])
    out["nodes"] = rows
    # ------------------------------------------------------------------ boot
    n1 = cfg["boot_edge"]
    bc = [(c[0] / 10, c[1], c[2] / 10) for c in lat.unpack_cells(data["boot"]["cells"], "313")]
    badm = [r for r in bc if r[1]]
    bs_m, bs_se = mean_se([r[2] for r in badm])
    alpha_b = max(C1000 - (bs_m - 3 * bs_se) * LN2, 0.0)
    pb = [D((r[0] * LN2 + alpha_b) / (n1 * LN2)) if r[1] else 0.0 for r in bc]
    pb_m, pb_se = mean_se(pb)
    nb_max = max(r[0] for r in bc)
    miss_b = missing(nb_max, n1, u_g)
    kids_b, cw_b, sub_b, lo = 0.0, 0.0, 0.0, b
    for e in EDGES:
        if lo >= n1:
            break
        k = children(nb_max, n1, lo, e)
        kids_b += k
        cw_b += k * W[e]
        sub_b += k * NODES[e]
        lo = e
    q_model = (pb_m - 2 * pb_se) * (1 - miss_b) * cfg["ecm_delta"]
    q_charged = q_model * cfg["boot_safety"]
    boot_trial = I_LLL + I_EXP + norm_instr(nb_max + 8) + norm_instr(4096) + trialdiv_instr(nb_max + 8) + ecm_instr(n1, nb_max + 8, cfg)
    attempt_model = cw_b + 2.0 ** 34 + I_SM            # model mean of the work of one descent
    attempt_charged = cfg["tree_safety"] * attempt_model
    # deterministic bound on the number of nodes of a descent: every boot part exceeds 2^b, a big part exceeds
    # 2^leaf_edge and has at most kmax children (each child exceeds 2^b and the norms are below the plan limits)
    kmax = max(r["max_children"] for r in rows if r["top_bits"] <= n1)
    limit_bits = nb_max + NORM_SLACK_BITS
    node_cap = 0
    for x in range(0, int(limit_bits // cfg["leaf_edge"]) + 1):
        y = int((limit_bits - cfg["leaf_edge"] * x) // b)
        node_cap = max(node_cap, x * (1 + kmax) + y)
    # H-desc (a): at most a fraction eps of the k-th visits is empty, for every k  =>  pi_d >= (1 - eps)^node_cap
    pi_model = (1.0 - cfg["eps_empty"]) ** node_cap
    pi_d = 2.0 ** math.floor(math.log2(pi_model))          # charged: rounded down to a power of two
    eps_law = math.exp(-cfg["lam_shifted_min"] / (cfg["mean_allowance"] * cfg["cluster_allowance"]))
    one_log = boot_trial / (q_charged * pi_d) + attempt_charged / pi_d
    two_logs = 2 * one_log
    lb1, curves, mulmods, gcds = ecm_schedule(n1, cfg)
    out["boot"] = dict(samples=len(bc), admissible_fraction=len(badm) / len(bc), alpha_boot=alpha_b,
                       norm_mean=sum(r[0] for r in bc) / len(bc), norm_max=nb_max,
                       log2_q_model=l2(q_model), log2_q_charged=l2(q_charged), missing=miss_b, mean_parts=kids_b,
                       tree_nodes_model=sub_b, max_children_of_a_big_node=kmax, node_cap=node_cap, pi_d_from_eps=pi_model, pi_d_charged=pi_d,
                       eps_from_compound_law=eps_law,
                       log2_trial_instr=l2(boot_trial),
                       ecm_log2_B1=lb1, ecm_curves=curves, ecm_log2_mulmods_per_curve=l2(mulmods),
                       log2_attempt_model_instr=l2(attempt_model), log2_attempt_charged_instr=l2(attempt_charged),
                       log2_one_log_instr=l2(one_log))
    # ------------------------------------------------------------------ totals
    fixed = collection + factor_rel + schirokauer + factor_base + la + I_WRAPPER
    run = fixed + two_logs
    gate = cfg["gate_num"] / cfg["gate_den"]
    post = 1.0 - cfg["main_fail"] - 1.0 / cfg["meter_factor"] - 2.0 ** -100
    total_units = 1 + gate * run / C_UNIT
    out["ledger_log2_instructions"] = dict(
        relation_collection=l2(collection), relation_factoring=l2(factor_rel), schirokauer_maps=l2(schirokauer),
        factor_base=l2(factor_base), linear_algebra=l2(la), two_individual_logarithms=l2(two_logs),
        wrapper=l2(I_WRAPPER), run_total=l2(run))
    out["log2_run_units"] = l2(run / C_UNIT)
    out["log2_logs_meter_instr"] = l2(cfg["meter_factor"] * two_logs)
    out["gate"] = f'{cfg["gate_num"]}/{cfg["gate_den"]}'
    out["post_gate_success"] = post
    out["success"] = gate * post
    out["log2_expected_total_units"] = l2(total_units)
    return out


if __name__ == "__main__":
    import sys
    d = load(Path(__file__).parent)
    o = ledger(d)
    json.dump(o, sys.stdout, indent=1)
