# Sign-law counts for proof.md sections 2.2 and 4: the secp256k1 sigma walk at A = 8, 64, 256, 1024 adders
# (not production size), with the forced-return split. Re-derives every number in sign-law.json from
# ../measurement/walk_ref.py and writes it to OUTDIR; compare with the committed file by cmp.
#
# Run from any directory:  python3 -B -I sign_law.py [OUTDIR]     (OUTDIR defaults to this folder)
import hashlib, importlib.util, json, math, multiprocessing, os, random, sys

sys.dont_write_bytecode = True
os.environ.update(ADDER_BITS="6", B="8", DP_BITS="30")
HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("walk_ref", os.path.join(HERE, "..", "measurement", "walk_ref.py"))
wr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(wr)
p, r, BETA, BETA2, G = wr.p, wr.r, wr.BETA, wr.BETA2, wr.G
LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
GAM = [1, r - 1, LAM, (r - LAM) % r, LAM * LAM % r, (r - LAM * LAM % r) % r]
GI = {g: n for n, g in enumerate(GAM)}
MINUS = 1  # position of -1 in GAM
STEPS = 1500
# (A, NW walks, seed, split): the runs quoted in proof.md; a smaller NW with the same seed is a prefix, so it is not listed
JOBS = [(8, 2048, 1, False), (64, 2048, 1, True), (64, 1024, 2, False), (64, 1024, 3, False),
        (256, 2048, 1, False), (256, 1024, 2, False), (256, 1024, 3, False),
        (1024, 4096, 1, False), (1024, 4096, 2, False)]


def eps(P):
    k, s, kk = wr.canon_info(P)
    return ((r - 1) if s else 1) * pow(LAM, k, r) % r, kk


def sha(obj):
    return hashlib.sha256(repr(obj).encode()).hexdigest()[:16]


def job(args):
    A, NW, seed, split = args
    rng = random.Random(seed)
    x = rng.randrange(1, r)
    gp = wr.comb(G)
    gq = wr.comb(wr.mul_comb(gp, x))
    ADD = [wr.add(wr.mul_comb(gp, rng.randrange(r)), wr.mul_comb(gq, rng.randrange(r))) for _ in range(A)]
    START = [wr.add(wr.mul_comb(gp, rng.randrange(r)), wr.mul_comb(gq, rng.randrange(r))) for _ in range(NW)]
    SR = [wr.canon_info(R)[:2] for R in ADD]  # (rotation k, sign s) of each adder
    c = dict(fresh=0, pairs=0, num=0, eq=0, eq_num=0, all=0, eq_all=0)
    sp = {"sign_same": [0, 0], "sign_diff": [0, 0], "rot_same": [0, 0], "rot_diff": [0, 0]}
    walk_eq, hist = [], {}
    for W0 in START:
        W, seen, seq = W0, set(), []
        for _ in range(STEPS):
            k, s, kk = wr.canon_info(W)
            fresh = kk not in seen
            seen.add(kk)
            i = (kk >> 1) & (A - 1)
            X, Y = ADD[i]
            img = ((X, BETA2 * X % p, BETA * X % p)[k], (p - Y) % p if s else Y)
            tau, _ = eps(wr.add(wr.canon(W), ADD[i]))
            seq.append((fresh, i, GI[tau]))
            W = wr.add(W, img)
        c["fresh"] += sum(1 for f, _, _ in seq if f)
        n_eq = 0
        for t in range(len(seq) - 1):
            f0, i0, u0 = seq[t]
            f1, i1, u1 = seq[t + 1]
            if not (f0 and f1):
                continue
            c["all"] += 1
            c["eq_all"] += i1 == i0
            if u0 != MINUS:
                continue
            m = u1 == MINUS
            c["pairs"] += 1
            c["num"] += m
            if i1 == i0:
                c["eq"] += 1
                c["eq_num"] += m
                n_eq += 1
            elif split:
                for name, a, b in (("sign", SR[i0][1], SR[i1][1]), ("rot", SR[i0][0], SR[i1][0])):
                    cell = name + ("_same" if a == b else "_diff")
                    sp[cell][0] += 1
                    sp[cell][1] += m
        walk_eq.append(n_eq)
        hist[n_eq] = hist.get(n_eq, 0) + 1
    out = dict(A=A, NW=NW, steps=STEPS, seed=seed, fresh_steps=c["fresh"], pairs=c["pairs"], m1m1=c["num"],
               eq_pairs=c["eq"], eq_m1m1=c["eq_num"], eq_violations=c["eq"] - c["eq_num"],
               nonEq_pairs=c["pairs"] - c["eq"], nonEq_m1m1=c["num"] - c["eq_num"],
               all_fresh_pairs=c["all"], all_eq_pairs=c["eq_all"],
               walk_eq_hist={str(k): v for k, v in sorted(hist.items())},
               walk_eq_sha256=sha(walk_eq), adders_sha256=sha(ADD), starts_sha256=sha(START))
    if split:
        out["split"] = sp
    return out


def stats(n, k, q):
    """Share k/n against q with the normal approximation, as used in proof.md."""
    sd = math.sqrt(q * (1 - q) / n)
    return dict(n=n, k=k, share=round(k / n, 6), excess_pp=round(100 * (k / n - q), 4),
                z=round((k / n - q) / sd, 2), ci95_pp=[round(100 * (k / n - q - 1.96 * sd), 4),
                                                       round(100 * (k / n - q + 1.96 * sd), 4)])


def dump(path, obj):
    with open(path, "w", encoding="utf-8") as f:
        f.write(json.dumps(obj, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    outdir = sys.argv[1] if len(sys.argv) > 1 else HERE
    with multiprocessing.get_context("fork").Pool(min(len(JOBS), os.cpu_count() or 1)) as pool:
        runs = pool.map(job, JOBS, chunksize=1)
    p6 = 1 / 6
    out = dict(purpose="Sign-law counts for proof.md sections 2.2 and 4, from sign_law.py. Small-A replica of the "
                        "sigma walk, not production size. A: per adder count A. split_A64: the A = 64 seed 1 run, "
                        "with its non-equal pairs split by the sign and the rotation of the two adders.",
               A={}, split_A64=None)
    for A in sorted({j[0] for j in JOBS}):
        rr = [o for o in runs if o["A"] == A]
        tot = {k: sum(o[k] for o in rr) for k in ("pairs", "m1m1", "eq_pairs", "eq_m1m1", "eq_violations",
                                                 "nonEq_pairs", "nonEq_m1m1", "fresh_steps")}
        n, num, eq = tot["pairs"], tot["m1m1"], tot["eq_pairs"]
        e = stats(n, num, p6)
        q = eq / n
        derived = dict(
            excess_share_based=dict(e, E_times_A=round(e["excess_pp"] / 100 * A, 4),
                                    E_times_A_ci95=[round(x / 100 * A, 4) for x in e["ci95_pp"]]),
            equal_index=dict(q=round(q, 6), A_times_q=round(A * q, 4),
                             A_times_q_ci95=[round(A * (q - 1.96 * math.sqrt(q * (1 - q) / n)), 4),
                                             round(A * (q + 1.96 * math.sqrt(q * (1 - q) / n)), 4)],
                             z_vs_1_over_A=round((q - 1 / A) / math.sqrt(q * (1 - q) / n), 2)),
            non_equal_index=stats(tot["nonEq_pairs"], tot["nonEq_m1m1"], p6),
            prediction_5_over_6A=round(5 / (6 * A), 8),
            excess_minus_prediction_z=round((e["excess_pp"] / 100 - 5 / (6 * A))
                                            / math.sqrt(p6 * (1 - p6) / n), 2))
        out["A"][str(A)] = dict(runs=[{k: v for k, v in o.items() if k != "split"} for o in rr],
                                pooled=tot, derived=derived)
        if A == 64:
            s = [o for o in rr if o["seed"] == 1][0]["split"]
            cells = {}
            for name in ("sign_same", "sign_diff", "rot_same", "rot_diff"):
                cells[name] = dict(stats(s[name][0], s[name][1], p6), k_m1m1=s[name][1])
            out["split_A64"] = dict(A=64, seed=1, pairs=[o for o in rr if o["seed"] == 1][0]["pairs"],
                                    equal_index_pairs=[o for o in rr if o["seed"] == 1][0]["eq_pairs"],
                                    equal_index_share=1.0, non_equal=cells)
    path = os.path.join(outdir, "sign-law.json")
    dump(path, out)
    print("written to", path)
