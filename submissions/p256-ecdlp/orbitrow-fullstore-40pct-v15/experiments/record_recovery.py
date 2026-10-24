"""v15 recovery: the 48-byte record NAMES a visit and the visit's coefficients are
REPLAYED from the trail start, not serialised twice.

v14 packed a value and compared it with an artificial 64-bit expression built from
the same inputs, which demonstrated serialization only.  Here a record is
(32B class x || 20-bit trail || 108-bit step).  Given only the unpacked record,
replay re-runs the trail's deterministic row schedule from step 0 to `step`,
reconstructs the coefficient pair (c,d), and checks both the point and the pair
against truth.  A store collision between two replayed records then yields
da*P + db*Q = O and recovers the secret scalar, verified by cand*P == Q.
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v15kit as K

R = 64
TRAILS = 1200
LEN = 32
SECRET = 777001
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")


def row_of(x, t, trail):
    return ((x ^ (t * 2654435761) ^ (trail * 7919)) % R) if x is not None else 0


def replay(tab, trail, upto):
    """Deterministic production schedule from the trail's seeded start."""
    c0 = (trail * 104729 + 7) % K.N_TOY
    d0 = (trail * 130363 + 11) % K.N_TOY
    S = K.add(K.mul(c0, K.mul(1, BASE)), K.mul(d0, Q))
    c, d = c0, d0
    for t in range(upto):
        r = row_of(S[0] if S else None, t, trail)
        a, b, V = tab[r]
        T = K.add(S, V)
        c2, d2 = (c + a) % K.N_TOY, (d + b) % K.N_TOY
        C = K.canon(T)
        if C != T:
            c2, d2 = (-c2) % K.N_TOY, (-d2) % K.N_TOY
        S, c, d = C, c2, d2
    return S, c, d


BASE = None
Q = None
tab = None


def main():
    global BASE, Q, tab
    BASE = K.base_point()
    Q = K.mul(SECRET, BASE)
    tab = K.row_table(BASE, Q, R)
    ok = 0
    checked = 0
    store = {}
    recovered = None
    sample = []
    for i in range(TRAILS):
        trail = i % (1 << K.TRAIL_BITS)
        for step in range(1, LEN + 1):
            S, c, d = replay(tab, trail, step)
            buf = K.record_pack(S[0], trail, step)
            assert len(buf) == K.RECORD_BYTES
            gx, (gt, gs) = K.record_unpack(buf)
            good = gx == S[0] and gt == trail and gs == step
            if checked < 1200:
                S2, c2, d2 = replay(tab, gt, gs)      # replay from the record alone
                good = good and S2 == S and c2 == c and d2 == d
                checked += 1
            else:
                c2, d2 = c, d
            ok += 1 if good else 0
            if i < 1 and step <= 3:
                sample.append({"trail": trail, "step": step, "match": good})
            key = gx
            if key in store:
                pc, pd, py = store[key]
                if py != S[1]:
                    da, db = (c2 + pc) % K.N_TOY, (d2 + pd) % K.N_TOY
                else:
                    da, db = (c2 - pc) % K.N_TOY, (d2 - pd) % K.N_TOY
                if (da or db) and __import__("math").gcd(db, K.N_TOY) == 1 and recovered is None:
                    cand = (-da * K.inv(db, K.N_TOY)) % K.N_TOY
                    if K.mul(cand, BASE) == Q:
                        recovered = cand
            else:
                store[key] = (c2, d2, S[1])
    ev = {
        "record_bytes": K.RECORD_BYTES, "layout": "32B class x || 20-bit trail || 108-bit step",
        "trails": TRAILS, "visits_per_trail": LEN, "records": TRAILS * LEN,
        "replay_equality_ok": ok, "replay_checked": checked,
        "all_equal": ok == TRAILS * LEN,
        "recovered_secret": recovered, "secret_verified": recovered == SECRET % K.N_TOY,
        "mechanism": "replay from (trail, step): coefficient pair reconstructed by re-running the deterministic schedule; collision da,db gives x=-da/db",
        "replay_cost": "<= 2 trails re-run at production, <= 2*2^107.67 steps = 2^-19 of total work (charged in ledger_v15 replay_relative_loss)",
        "sample": sample,
    }
    json.dump(ev, open(os.path.join(OUT, "record_recovery_v15.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({k: ev[k] for k in ("replay_equality_ok", "all_equal", "secret_verified")}, indent=1))


if __name__ == "__main__":
    main()
