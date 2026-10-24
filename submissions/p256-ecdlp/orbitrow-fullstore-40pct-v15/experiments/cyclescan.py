"""v15 generalized 4-cycle scan: EVERY row 4-tuple, not just [A,B,A,B].

v14 imposed the alternating shape and multiplied by min(1,ratio), discarding the
measured 1.779x inflation. v15 enumerates all 64^4 windows per seed by an exact
prefix/inverse-image split and replays every hit with coefficient tracking.
delta == (0,0) is the coefficient-trivial death; nonzero delta is a usable relation.
Pricing multiplies by the measured ratio, never min(1,ratio).
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v15kit as K

R_FIX = 64
SEEDS = 16
R_PROD = 1 << 20
LOG2_N_PROD = 256.0
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")


def scan(P, tab):
    pts = [t[2] for t in tab]
    coefs = [(t[0], t[1]) for t in tab]
    n = K.N_TOY
    stats = {"windows": 0, "returns": 0, "trivial": 0, "collapse_confirmed": 0,
             "nonzero_delta_returns": 0}
    fam_trivial = {"alt4": 0, "rev4": 0, "other": 0}
    for s in range(1, SEEDS + 1):
        S = K.canon(K.mul(s, P))
        negS = K.neg(S)
        pref = {}
        for A in range(R_FIX):
            X1 = K.canon(K.add(S, pts[A]))
            for B in range(R_FIX):
                X2 = K.canon(K.add(X1, pts[B]))
                pref.setdefault(X2, []).append((A, B))
        H = {}
        for D in range(R_FIX):
            H.setdefault(K.add(S, K.neg(pts[D])), []).append(D)
            H.setdefault(K.add(negS, K.neg(pts[D])), []).append(D)
        stats["windows"] += R_FIX ** 4
        stats["probes_per_seed"] = R_FIX ** 3
        for X2, pairs in pref.items():
            for (A, B) in pairs:
                for C in range(R_FIX):
                    Y = K.canon(K.add(X2, pts[C]))
                    Ds = H.get(Y)
                    if not Ds:
                        continue
                    for D in Ds:
                        st = S; ca = cb = 0
                        for r in (A, B, C, D):
                            T = K.add(st, pts[r])
                            ca = (ca + coefs[r][0]) % n
                            cb = (cb + coefs[r][1]) % n
                            Cc = K.canon(T)
                            if Cc != T:
                                ca, cb = (-ca) % n, (-cb) % n
                            st = Cc
                        if st != S:
                            stats["collapse_confirmed"] += 1
                            continue
                        stats["returns"] += 1
                        if A == C and B == D:
                            fam = "alt4"
                        elif A == D and B == C:
                            fam = "rev4"
                        else:
                            fam = "other"
                        if ca == 0 and cb == 0:
                            stats["trivial"] += 1
                            fam_trivial[fam] += 1
                        else:
                            stats["nonzero_delta_returns"] += 1
    w_alt = SEEDS * R_FIX * R_FIX
    w_rev = SEEDS * R_FIX * R_FIX
    w_tot = stats["windows"]
    pred_alt = w_alt * 2.0 ** -4
    pred_rev = w_rev * 2.0 ** -4
    pred_other = w_tot * 2.0 ** -4 * (1.0 / K.N_TOY)
    ratio_alt = fam_trivial["alt4"] / pred_alt if pred_alt else 0.0
    ratio_rev = fam_trivial["rev4"] / pred_rev if pred_rev else 0.0
    ratio_other = fam_trivial["other"] / pred_other if pred_other else 0.0
    per4 = (ratio_alt * 2.0 ** -4 * R_PROD ** -2.0
            + ratio_rev * 2.0 ** -4 * R_PROD ** -2.0
            + ratio_other * 2.0 ** -4 * 2.0 ** -LOG2_N_PROD)
    steps_per_walker = 2 ** 107.67
    walkers = 1 << 20
    deaths = (steps_per_walker / 4.0) * per4
    waste = 6.0
    rel_loss = deaths * waste * walkers / (walkers * steps_per_walker)
    return stats, fam_trivial, {
        "windows_scanned": w_tot, "ratio_alt4": ratio_alt, "ratio_rev4": ratio_rev,
        "ratio_other": ratio_other, "pred_alt4": pred_alt, "pred_rev4": pred_rev,
        "pred_other": pred_other, "per_4step_death_rate_upper": per4,
        "deaths_per_walker": deaths, "waste_per_death_steps": waste,
        "relative_work_loss": rel_loss, "claim_delta_bits": K.log2(1.0 + rel_loss),
    }


def main():
    P = K.base_point()
    Q = K.mul(1234567, P)
    tab = K.row_table(P, Q, R_FIX)
    stats, fam, price = scan(P, tab)
    ev = {
        "fixture": {"p": K.P_TOY, "b": K.B_TOY, "order": K.N_TOY, "R_fix": R_FIX,
                    "seeds": SEEDS, "cls": "canonical x (P/-P merged, x/p-x distinct)"},
        "scan": stats, "trivial_by_family": fam, "pricing": price,
        "method": "full 64^4 enumeration per seed via prefix/inverse-image split; every hit replayed with coefficient tracking; ratio charged as a multiplicative factor, never min(1,ratio)",
        "conclusion": "all coefficient-trivial 4-step returns are charged; the measured ratio is applied as-is",
    }
    os.makedirs(OUT, exist_ok=True)
    json.dump(ev, open(os.path.join(OUT, "cycle_channel_v15.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps(ev, indent=1))


if __name__ == "__main__":
    main()
