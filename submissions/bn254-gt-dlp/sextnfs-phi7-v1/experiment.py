"""Declared experiment for the SexTNFS (6,2) claim over Q(zeta_7) on BN254 GT.

Runs inside the pinned verifier's Docker sandbox (python3 -I, stdlib only, no network).
It does NOT run a production sieve. It (1) verifies the tower, twist, common root and the
shipped f-irreducibility (Capelli) certificate at full size; (2) re-verifies a reviewer-chosen
handful of the shipped measurements (production-cell norms with unit multiplicity, one
production special-q lattice sample, the toy-box calibration region with a replay, a few ECM
trials); (3) recomputes every summary statistic from data.json and the complete ledger of
machine.md from the protected constants, and reports the rounded-up score.

Every check is reported as a boolean in "checks"; the program never aborts on a numeric
mismatch, so a transient issue cannot masquerade as a mathematical refutation.
"""
import importlib.util
import json
import math
import random
import sys
import time
from decimal import Decimal, ROUND_CEILING, getcontext
from pathlib import Path

T0 = time.time()
base = Path("/snapshot/files")
if not base.is_dir():
    base = Path(__file__).parent
spec = importlib.util.spec_from_file_location("tnfs_common", base / "tnfs_common.py")
tc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tc)
req = json.load(sys.stdin)
seed = int(req["seed"])
recheck_cells = int(req["recheck_cells"])
spq_index = int(req["spq_index"])
toy_index = int(req["toy_index"])
ecm_index = int(req["ecm_index"])
coins = random.Random(seed)
profile = json.loads((base / "target-profile.json").read_text())
data = json.loads((base / "data.json").read_text())
for extra in ("data_cells.json", "data_spq.json", "data_toy.json"):   # data is split to stay under the per-file byte limit
    data.update(json.loads((base / extra).read_text()))
L2 = math.log2
LN2 = math.log(2)
checks = {}

def ok(name, cond):
    checks[name] = bool(cond)
    return bool(cond)

# ---------------------------------------------------------------- protected constants and tower
p = int(profile["field_modulus"], 16)
r = int(profile["subgroup_order"], 16)
u = int(profile["parameter_x"], 16)
ok("bn_P", tc.bn_P(u) == p and tc.bn_R(u) == r)
ok("gt_cyclotomic", (p**4 - p**2 + 1) % r == 0 and all((p**k - 1) % r != 0 for k in (1, 2, 3, 4, 6)))
EPS = data["parameters"]["EPS"]
ok("eps_value", EPS == [2, 1, 0, 0, 0, 0])
CTX = tc.make_context(p, u, EPS)
ok("tower_and_twist", tc.check_context(CTX))          # Phi_7 inert, u*eps non-square, u square, common root
ok("eps_norm_43", tc.zy_norm(EPS) == 43)
ok("rho_table", tc.rho_check())
cert = data["f_irreducibility_certificate"]
cq, cs, cth = cert
es = tc.zy_eval(EPS, cs, cq)
ok("capelli_prime", cq % 7 == 1 and tc.is_prime(cq) and pow(cs, 7, cq) == 1 and cs != 1)
ok("capelli_theta_is_P_root", (36*cth**4 + 36*cth**3 + 24*cth**2 + 6*cth + 1) % cq == 0)
ok("capelli_nonresidue", tc.legendre(es*cth, cq) == -1)   # eps*theta a non-square in K_h(theta) => f irreducible deg 8
A = data["parameters"]["A"]
ok("A_255", A == 255)
T_tower = time.time() - T0

# ---------------------------------------------------------------- re-verify shipped data (handfuls)
cells = data["production_cells"]
_U = None
def unit_list():
    global _U
    if _U is None:
        import itertools
        u1 = [1, 1, 0, 0, 0, 0]; u2 = [1, 1, 1, 0, 0, 0]
        def inv(w):
            for v in itertools.product(range(-2, 3), repeat=6):
                if tc.zy_mul(w, list(v)) == tc.ONE:
                    return list(v)
            raise AssertionError
        u1i, u2i = inv(u1), inv(u2)
        out = []
        for k in range(7):
            zk = tc.ONE
            for _ in range(k):
                zk = tc.zy_mul(zk, tc.Y)
            for i in range(-2, 3):
                for j in range(-2, 3):
                    w = zk
                    for _ in range(abs(i)):
                        w = tc.zy_mul(w, u1 if i > 0 else u1i)
                    for _ in range(abs(j)):
                        w = tc.zy_mul(w, u2 if j > 0 else u2i)
                    out.append(w)
        _U = out
    return _U
cell_ok = True
for i in coins.sample(range(len(cells)), min(recheck_cells, len(cells))):
    c = cells[i]
    nf, ng = tc.norm_pair(CTX, c["a"], c["b"])
    m = sum(1 for w in unit_list() if max(abs(v) for v in tc.zy_mul(w, c["a"])) <= A and max(abs(v) for v in tc.zy_mul(w, c["b"])) <= A)
    cell_ok = cell_ok and abs(L2(nf) - c["lf"]) < 0.01 and abs(L2(ng) - c["lg"]) < 0.01 and m == c["mult"] and all(abs(v) <= A for v in c["a"] + c["b"])
ok("production_cell_norms", cell_ok)
T_cells = time.time() - T0 - T_tower

spq = data["production_special_q"]
rec = spq[spq_index % len(spq)]
rnd = random.Random(rec["seed"])
q, s, t = tc.special_q(CTX, rec["bits"], rec["side"], rnd)
spq_ok = (str(q), str(s), str(t)) == (rec["q"], rec["s"], rec["t"]) and q % 7 == 1 and tc.is_prime(q)
if rec["side"] == "g":
    spq_ok = spq_ok and (t*t - u*tc.zy_eval(EPS, s, q)) % q == 0
else:
    spq_ok = spq_ok and sum(tc.zy_eval(cf, s, q)*pow(t, k, q) for k, cf in CTX["FCO"].items()) % q == 0
R = tc.lll_int(tc.spq_lattice(q, s, t))
spq_ok = spq_ok and max(abs(x) for v in R for x in v).bit_length() == rec["basis_max_bits"]
spq_ok = spq_ok and all(tc.in_spq_lattice(q, s, t, v[:6], v[6:]) for v in R)
rows = []
for E in (1, 2, 4, 8, 16, 32, 64):
    for _ in range(30):
        cvec = [rnd.randint(-E, E) for _ in range(12)]
        if not any(cvec):
            cvec[0] = 1
        a, b = tc.combine(R, cvec)
        nf, ng = tc.norm_pair(CTX, a, b)
        if rec["side"] == "g":
            spq_ok = spq_ok and ng % q == 0; ng //= q
        else:
            spq_ok = spq_ok and nf % q == 0; nf //= q
        rows.append([E, round(L2(max(nf, 1)), 2), round(L2(max(ng, 1)), 2)])
spq_ok = spq_ok and rows == rec["rows"]
ok("production_special_q_sample", spq_ok)
T_spq = time.time() - T0 - T_tower - T_cells

# toy-box calibration: re-verify one region's shipped sample and replay to match its counts
TP = data["parameters"]
TOY = tc.make_context(tc.bn_P(TP["toy_u"]), TP["toy_u"], TP["toy_eps"])
ok("toy_tower", tc.bn_P(TP["toy_u"]) == TP["toy_p"] and tc.is_prime(TP["toy_p"]) and tc.check_context(TOY))
P20 = tc.primes_upto(1 << 20)
treg = data["toy_box"][toy_index % len(data["toy_box"])]
rnd = random.Random(treg["seed"])
toy_ok = True
npre = treg["prefix_n"]
counts = {int(b): 0 for b in TP["toy_B_bits"]}
shipped_sample = {i: s_ for i, s_ in enumerate(treg["sample"])}
for i in range(npre):                            # replay only the shipped prefix (keeps the sandbox fast)
    a = [rnd.randint(-1, 1) for _ in range(6)]; b = [rnd.randint(-1, 1) for _ in range(6)]
    if not any(b):
        b[0] = 1
    nf, ng = tc.norm_pair(TOY, a, b)
    mf, okf = tc.max_prime_factor_bits(nf, P20); mg, okg = tc.max_prime_factor_bits(ng, P20)
    if i in shipped_sample:
        sm = shipped_sample[i]
        toy_ok = toy_ok and sm["a"] == a and sm["b"] == b and sm["mf"] == (mf if okf else -1) and sm["mg"] == (mg if okg else -1)
    for bb in TP["toy_B_bits"]:
        if okf and okg and mf <= bb and mg <= bb:
            counts[int(bb)] += 1
toy_ok = toy_ok and all(counts[int(bb)] == treg["counts_prefix"][str(bb)] for bb in TP["toy_B_bits"])
ok("toy_box_prefix", toy_ok)
T_toy = time.time() - T0 - T_tower - T_cells - T_spq

eset = data["ecm"][ecm_index % len(data["ecm"])]
rnd = random.Random(eset["trials"][0]["p1"] ^ seed)   # choose a few trials to re-run deterministically
replay_seed = 13000 + (ecm_index % 12) if eset["fbits"] == 24 else 14000 + (ecm_index % 6)
rnd = random.Random(replay_seed)
pl = tc.primes_upto(eset["B1"] + 1)
ecm_ok = True
for j, trial in enumerate(eset["trials"]):
    p1 = next(x for x in range((1 << (eset["fbits"]-1)) + (rnd.getrandbits(eset["fbits"]-1) | 1), 1 << (eset["fbits"]+1), 2) if tc.is_prime(x))
    p2 = next(x for x in range((1 << 59) + (rnd.getrandbits(59) | 1), 1 << 61, 2) if tc.is_prime(x))
    match = (p1, p2) == (trial["p1"], trial["p2"])
    if j < 4:                       # re-run only a few curves to stay within the time budget
        g, used = tc.ecm_schedule(p1*p2, eset["B1"], eset["curves"], pl, rnd)
        match = match and (g is not None) == trial["found"] and used == trial["curves"] and (g in (p1, p2) if g else True)
    else:
        for _ in range(trial["curves"]):
            rnd.randrange(6, p1*p2)   # advance the stream to mirror the shipped run
    ecm_ok = ecm_ok and match
ok("ecm_trials", ecm_ok)
T_ecm = time.time() - T0 - T_tower - T_cells - T_spq - T_toy

# ---------------------------------------------------------------- statistics from the shipped data
def mean(v):
    return sum(v)/len(v)
def se(v):
    mu = mean(v); return math.sqrt(sum((x-mu)**2 for x in v)/(len(v)-1)/len(v)) if len(v) > 1 else 0.0
lf = [c["lf"] for c in cells]; lg = [c["lg"] for c in cells]
alpha_f = [c["af"] for c in cells]; alpha_g = [c["ag"] for c in cells]
alpha_f_ch = max(0.0, mean(alpha_f) + 2*se(alpha_f)); alpha_g_ch = max(0.0, mean(alpha_g) + 2*se(alpha_g))
p_prim = sum(1 for c in cells if c["prim"])/len(cells)
inv_mult = mean([1.0/c["mult"] for c in cells])
mult_edge = max(c["edge"] for c in cells)
boot_bits = [b["best_bits"] for b in data["production_boot"]]
spq_stat = {}
for recd in spq:
    for E, a_, b_ in recd["rows"]:
        spq_stat.setdefault((recd["bits"], recd["side"], E), []).append((a_, b_))

# toy calibration: realised/predicted ratio (with a 2-sigma lower bound) and the dispersion index
tb = data["toy_box"]
Bcal = max(TP["toy_B_bits"])
reg_counts = [t_["counts"][str(Bcal)] for t_ in tb]
reg_pred = [t_["predicted"][str(Bcal)] for t_ in tb]
tot_c = sum(reg_counts); tot_p = sum(reg_pred)
ratio = tot_c/tot_p
ratio_se = math.sqrt(tot_c)/tot_p          # Poisson count error on the realised total
cal_yield = min(1.0, ratio - 2*ratio_se)   # never credit a yield bonus; penalise only if rho over-predicts
mu_c = mean(reg_counts); var_c = sum((x-mu_c)**2 for x in reg_counts)/(len(reg_counts)-1)
dispersion = var_c/mu_c
DISP = max(2.0, 2*dispersion)              # charged upper bound on Var/mean for the node zero-probability
zero_frac = sum(1 for x in reg_counts if x == 0)/len(reg_counts)
ecm_meas = {}
for e_ in data["ecm"]:
    ecm_meas.setdefault((e_["fbits"], e_["B1"]), []).extend(e_["trials"])
ecm_summary = {}
for (fb, b1), trials in ecm_meas.items():
    found = sum(1 for t_ in trials if t_["found"])
    rate = found/max(1, sum(t_["curves"] for t_ in trials))
    ecm_summary["%db_B1_%d" % (fb, b1)] = dict(trials=len(trials), found=found, per_curve_measured=round(rate, 5),
                                               per_curve_predicted=round(tc.rho(fb/L2(b1)), 5))
ecm_cal = min(1.0, min(v["per_curve_measured"]/v["per_curve_predicted"] for v in ecm_summary.values() if v["per_curve_predicted"] > 0))

# ---------------------------------------------------------------- ledger (scalar AArch64 instructions)
PAR = dict(logB=54, excess=1.25, SM=28, sample_margin_bits=0.5, survivor_mult=16,
           instr_hit=12, instr_hit_large=40, instr_setup=256, instr_scan=4, instr_cofactor_mulmod=64,
           instr_entry=20, instr_sm=600, instr_rowfix=160, matvec_factor=3, w_small=130,
           boot_logB1_bits=34, ecm_fail_bits=10, mu_node_bits=10, retries=4, levels=[(200, 150), (150, 110), (110, 80), (80, 54)])
logB = PAR["logB"]; B = 2.0**logB; lnB = logB*LN2
FB_side = B/math.log(B)*(1 + 1/math.log(B) + 2/math.log(B)**2)
D = 2*FB_side + PAR["SM"]
yjoint = mean([tc.rho((lf[i] + alpha_f_ch/LN2)/logB)*tc.rho((lg[i] + alpha_g_ch/LN2)/logB) for i in range(len(cells))])
yield_cell = yjoint*p_prim*inv_mult*cal_yield*2.0**(-PAR["sample_margin_bits"])
cells_total = (2*A+1)**12/2
relations = cells_total*yield_cell
need = PAR["excess"]*D
ok("relations_sufficient", relations >= need)
hits_side = math.log(lnB) + 0.2615 + 0.3
hits_large = math.log(logB/36.0)
instr_cell = 2*(hits_side*PAR["instr_hit"] + hits_large*PAR["instr_hit_large"]) + PAR["instr_scan"]

def mulmod_instr(limbs):
    return 4.0*limbs*limbs + 32.0*limbs + 64

def ecm_cost(fac_bits, b1_bits, limbs):
    """stage-1 ECM: per-curve success rho(fac_bits/b1_bits)*measured-calibration, curves for 2^-ecm_fail_bits failure."""
    pc = max(1e-300, tc.rho(fac_bits/b1_bits)*ecm_cal)
    curves = math.ceil(PAR["ecm_fail_bits"]*LN2/pc)
    per_curve = 1.44*2.0**b1_bits*11*mulmod_instr(limbs) + 2.0**16
    return curves, curves*per_curve

def ecm_cost_opt(fac_bits, limbs):
    """ECM with the near-optimal stage-1 bound for a factor of fac_bits bits."""
    f = fac_bits*LN2
    b1_bits = max(16.0, math.sqrt(f*math.log(f)/2)/LN2)
    return ecm_cost(fac_bits, b1_bits, limbs)

blocks = (2*A+1)**6/2
setup = blocks*D*PAR["instr_setup"]
sieve = cells_total*instr_cell + setup
lll_fb = D*2.0**20
# survivors = cells whose sieved residue is below the smoothness threshold; charged survivor_mult x relations
survivors = PAR["survivor_mult"]*relations
curves54, sched54 = ecm_cost(54, 20, 8)
rounds_main = math.ceil((max(lf) + max(lg))/20)
ecm_main = survivors*(rounds_main*sched54 + 2.0**20)
sm_cost = relations*(2.0**32.4 + 2.0**27.7)
collect = sieve + lll_fb + ecm_main + sm_cost
row = PAR["w_small"]*PAR["instr_entry"] + PAR["SM"]*PAR["instr_sm"] + PAR["instr_rowfix"]
la = PAR["matvec_factor"]*D*D*row
# boot: smooth the g-side lift of the target to bound 2^B_boot (= top descent level). Rigorous bound on
# the g-side boot norm: LLL on the dim-12 lattice of det p^6 gives coords <= 2^(2.75)*p^(1/2); each of the
# 6 conjugates of a^2 - u*eps*b^2 is <= (6*2^(2.75+127))^2 * (1 + 3u).
BOOT_RIG = 6*(2*(L2(6) + 2.75 + 127.0) + L2(1 + 3*u))
ok("boot_lll_within_rigorous", max(boot_bits) < BOOT_RIG)
B_boot_bits = 200                                   # top descent level; boot large primes are <= 2^200
limbs_boot = math.ceil(BOOT_RIG/64)
p_boot = tc.rho(BOOT_RIG/B_boot_bits)*cal_yield
K_boot = math.ceil(8*LN2/max(p_boot, 1e-300))       # candidates for 2^-8 boot failure
curves_boot, sched_boot = ecm_cost_opt(B_boot_bits, limbs_boot)
large_primes_boot = math.ceil(BOOT_RIG/logB)        # every prime factor above the factor-base bound 2^54 is a descent child
cand_boot = 2.0**22 + large_primes_boot*sched_boot
boot_total = 2*K_boot*cand_boot                     # two boots (g and z)
boot_fail = 2*2.0**-8                               # each boot fails with probability <= (1-p)^K <= exp(-K p) <= 2^-8

# descent: per level, pick the worst side; node = sieve the ball + ECM on every cell above the factor-base level
ladder = []
nodes = 2*large_primes_boot                       # the two boots seed this many top-level special-q's
desc_total = 0.0; desc_fail = 0.0; total_nodes = 0.0
for frm, to in PAR["levels"]:
    worst = None
    for side in ("g", "f"):
        chosen = None
        for E in (1, 2, 4, 8, 16, 32, 64):
            smp = spq_stat.get((frm, side, E))
            if not smp:
                continue
            fb = max(x[0] for x in smp) + 2.0; gb = max(x[1] for x in smp) + 2.0
            pcell = mean([tc.rho(x[0]/to)*tc.rho(x[1]/to) for x in smp])*cal_yield
            cells_E = (2*E+1)**12
            if cells_E*pcell >= 2.0**PAR["mu_node_bits"]:
                chosen = (E, fb, gb, pcell, cells_E); break
        if chosen is None:
            E = 64; smp = spq_stat[(frm, side, E)]
            fb = max(x[0] for x in smp) + 2.0; gb = max(x[1] for x in smp) + 2.0
            pcell = mean([tc.rho(x[0]/to)*tc.rho(x[1]/to) for x in smp])*cal_yield
            cells_E = (2*E+1)**12
        else:
            E, fb, gb, pcell, cells_E = chosen
        mu = cells_E*pcell
        rounds = math.ceil((fb + gb)/logB)    # the ball is sieved to 2^54 first; each ECM round strips one factor above 2^54
        if to > 54:
            b1 = max(16, int(round(to*0.17))); limbs = math.ceil(max(fb, gb)/64)
            cv, sched = ecm_cost(to, b1, limbs)
            per_surv = rounds*sched + 2.0**20
        else:
            cv, sched, b1 = 0, 0.0, 0
            per_surv = instr_cell + 2.0**10
        survs = cells_E if to > logB else PAR["survivor_mult"]*max(mu, 1.0)   # levels above 2^54: ECM on EVERY cell (no survivor filter)
        attempt = cells_E*instr_cell + survs*per_surv + 2.0**24
        p_fail1 = min(0.5, DISP/max(mu/4.0, 1.0))   # one attempt = a disjoint translate of the box, charged at a quarter of the centred mean count
        node = PAR["retries"]*attempt         # fixed retry cap -> hard per-node work bound on every run
        p_fail = p_fail1**PAR["retries"]      # all retries fail
        # children: every prime factor above the factor-base bound 2^54 is descended; each exceeds 2^54, so at most (fb+gb)/54
        children = max(1, math.ceil((fb + gb)/logB))
        info = dict(side=side, E=E, cells=round(L2(cells_E), 2), max_log2_Nf=round(fb, 1), max_log2_Ng_cof=round(gb, 1),
                    p_cell=round(L2(max(pcell, 1e-300)), 2), mu=round(L2(max(mu, 1e-300)), 2), B1_bits=b1, curves_per_survivor=cv,
                    per_survivor=round(L2(per_surv), 2), node=round(L2(node), 2), children=children, p_fail=round(p_fail, 4))
        if worst is None or node > worst["_node"]:
            worst = dict(info, _node=node, _children=children, _pfail=p_fail)
    lvl = {k: v for k, v in worst.items() if not k.startswith("_")}
    lvl["level"] = "%d->%d" % (frm, to); lvl["nodes"] = round(L2(nodes), 2); lvl["cost"] = round(L2(nodes*worst["_node"]), 2)
    desc_total += nodes*worst["_node"]
    desc_fail += nodes*worst["_pfail"]
    total_nodes += nodes
    ladder.append(lvl)
    nodes = nodes*worst["_children"]
fixed = 2.0**50
R_total = collect + la + boot_total + desc_total + fixed
W = R_total/20000.0
score = L2(W)
getcontext().prec = 50
score_up = Decimal(score).quantize(Decimal("0.000001"), rounding=ROUND_CEILING)
# 0.05: allowance for linear-algebra rank deficiency (H4) and missed ECM factors; 0.02: Schirokauer-map deficiency (H5)
success = 1 - (1 - cal_yield if cal_yield < 1 else 0) - 0.05 - boot_fail - desc_fail - 0.02
CLAIM = Decimal("108.000000")                     # the submitted bound (submission.json); margin disclosed in proof.md
ok("derived_score_within_claim", score_up <= CLAIM)
ok("success_at_least_0.39", success >= 0.39)
out = {
    "inputs": req,
    "checks": checks,
    "all_checks_pass": all(checks.values()),
    "tower": {"ord_7_p": 6, "eps": "2+y", "norm_eps": 43, "capelli_certificate": cert},
    "reverified": {"production_cells": min(recheck_cells, len(cells)),
                   "special_q": {"index": spq_index % len(spq), "bits": rec["bits"], "side": rec["side"], "q_bits": q.bit_length()},
                   "toy_box": {"index": toy_index % len(tb), "cells": treg["n"]},
                   "ecm": {"index": ecm_index % len(data["ecm"]), "trials": len(eset["trials"]), "curves_rerun": min(4, len(eset["trials"]))},
                   "times_s": {"tower": round(T_tower, 2), "cells": round(T_cells, 2), "spq": round(T_spq, 2), "toy": round(T_toy, 2), "ecm": round(T_ecm, 2)}},
    "production_cells": {"count": len(cells), "A": A, "mean_log2_Nf": round(mean(lf), 2), "mean_log2_Ng": round(mean(lg), 2),
                         "max_log2_Nf": round(max(lf), 2), "max_log2_Ng": round(max(lg), 2),
                         "alpha_f_charged_nats": round(alpha_f_ch, 3), "alpha_g_charged_nats": round(alpha_g_ch, 3),
                         "primitive_fraction": round(p_prim, 4), "mean_inverse_unit_multiplicity": round(inv_mult, 4),
                         "max_unit_exponent_in_box": mult_edge},
    "production_boot": {"lll_best_norm_bits": boot_bits, "rigorous_bits_charged": round(BOOT_RIG, 1)},
    "production_special_q": {"%d_%s_E%d" % k: {"n": len(v), "mean_f": round(mean([x[0] for x in v]), 1), "max_f": round(max(x[0] for x in v), 1),
                                               "mean_gcof": round(mean([x[1] for x in v]), 1), "max_gcof": round(max(x[1] for x in v), 1)}
                             for k, v in sorted(spq_stat.items())},
    "toy_calibration": {"p_bits": TP["toy_p"].bit_length(), "u": TP["toy_u"], "B_bits": TP["toy_B_bits"], "regions": len(tb), "cells_per_region": TP["toy_cells"],
                        "calibration_B_bits": Bcal, "realised_total": tot_c, "predicted_total": round(tot_p, 2),
                        "ratio": round(ratio, 3), "ratio_2sigma_low": round(ratio - 2*ratio_se, 3), "charged_yield_factor": round(cal_yield, 4),
                        "mean_count": round(mu_c, 2), "dispersion_index": round(dispersion, 3), "charged_dispersion": DISP, "zero_fraction": round(zero_frac, 4)},
    "ecm_calibration": {"measured": ecm_summary, "charged_factor": round(ecm_cal, 4)},
    "ledger_log2_instructions": {
        "factor_base_per_side": round(L2(FB_side), 3), "columns_D": round(L2(D), 3),
        "yield_joint_rho": round(L2(yjoint), 3), "yield_cell_charged": round(L2(yield_cell), 3),
        "cells_in_region": round(L2(cells_total), 3), "relations": round(L2(relations), 3), "relations_needed": round(L2(need), 3),
        "instr_per_cell": round(instr_cell, 1), "sieve": round(L2(sieve), 3), "setup_included": round(L2(setup), 3),
        "ecm_survivors": round(L2(ecm_main), 3), "curves_per_54b_survivor": curves54, "schirokauer": round(L2(sm_cost), 3), "collection": round(L2(collect), 3),
        "la_row_weight_instr": row, "linear_algebra": round(L2(la), 3),
        "boot": {"norm_bits_rigorous": round(BOOT_RIG, 1), "norm_bits_measured_max": round(max(boot_bits), 1), "B_boot_bits": B_boot_bits,
                 "candidates_log2": round(L2(K_boot), 2), "p_candidate_log2": round(L2(max(p_boot, 1e-300)), 2), "curves_per_cofactor": curves_boot,
                 "large_primes_each": large_primes_boot, "per_candidate_log2": round(L2(cand_boot), 2), "total": round(L2(boot_total), 3)},
        "descent": {"levels": ladder, "nodes_total_log2": round(L2(total_nodes), 2), "total": round(L2(desc_total), 3)},
        "fixed_setup_verify_output": 50.0,
        "total_instructions_log2": round(L2(R_total), 4), "W_total_units_log2": round(score, 6), "score_rounded_up": str(score_up),
        "claimed_bound": str(CLAIM), "disclosed_margin_bits": round(float(CLAIM) - score, 3),
    },
    "success_lower_bound": round(success, 4),
    "elapsed_s": round(time.time() - T0, 2),
    "scope": "Full-size tower/twist/irreducibility checks, exact norms of sampled production cells and special-q balls, toy-scale realised joint-smoothness calibration (ratio and dispersion), small-scale ECM detection, and ledger recomputation. No production sieve, DLP or verifier acceptance is claimed.",
}
print(json.dumps(out, indent=1))
