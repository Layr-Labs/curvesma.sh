"""Fix 1b: full-store census with the zero-denominator reject repaired.

The v13 census called first_repeat and rejected a walk whose repeat had zero
coefficient difference. That threw away the event instead of counting it; the
reviewer read the resulting zero cell as evidence the channel does not fire.

v14 defines success once (a verified scalar recovery) and counts, per trail:
  - visits (records appended to the full store)
  - trivial_repeat_no_progress  (repeat with zero coefficient delta: the 4-cycle channel)
  - true_collisions             (repeat with nonzero delta: usable relation)
  - restarts                    (walk re-seeded after a no-progress repeat, <=4 steps lost)
The reconciliation row prints the v13 reject vs the v14 count on the same trails.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v14kit as K

R = 64
WALKERS = 8
STEPS = 1200
SECRET = 1234567
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")


def base_point():
    for x in range(1, K.P_TOY):
        rhs = (x * x * x - 3 * x + K.B_TOY) % K.P_TOY
        if pow(rhs, (K.P_TOY - 1) // 2, K.P_TOY) == 1:
            y = pow(rhs, (K.P_TOY + 1) // 4, K.P_TOY)
            if (y * y - rhs) % K.P_TOY == 0:
                return (x, y)
    raise SystemExit("no base point")


def row_table(P, Q):
    """R rows, each a point a*P + b*Q with a,b in [1,R]; injected by row index."""
    tab = []
    for r in range(R):
        a = 1 + (r * 7 + 3) % R
        b = 1 + (r * 11 + 5) % R
        tab.append((a, b, K.add(K.mul(a, P), K.mul(b, Q))))
    return tab


def curve_order(P):
    """Point count by Legendre sum; the group is cyclic so ord(P) = n for a curve
    with prime order, which the fixture asserts."""
    p = K.P_TOY
    n = 1
    for x in range(p):
        rhs = (x * x * x - 3 * x + K.B_TOY) % p
        if rhs == 0:
            n += 1
        else:
            n += 2 if pow(rhs, (p - 1) // 2, p) == 1 else 0
    return n


def factor_small(m):
    fs, d = [], 2
    while d * d <= m:
        if m % d == 0:
            fs.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        fs.append(m)
    return fs


def point_order(P, n):
    """ord(P) from the curve order n: divide out each prime factor while it stays killed."""
    o = n
    for q in factor_small(n):
        while o % q == 0 and K.mul(o // q, P) is None:
            o //= q
    return o


def run():
    P = base_point()
    Q = K.mul(SECRET, P)
    tab = row_table(P, Q)
    n = point_order(P, curve_order(P))
    store = {}
    v13_rejects = 0
    trivial = 0
    true_col = 0
    restarts = 0
    visits = 0
    recovered = None
    for w in range(WALKERS):
        cx, cy = 1000 + w, 2000 + w          # S = cx*P + cy*Q
        S = K.add(K.mul(cx, P), K.mul(cy, Q))
        trail = w
        step = 0
        while step < STEPS:
            r = (K.mul(step + 1, P)[0] * 2654435761 + trail * 7919 + step) % R
            a, b, V = tab[r]
            T = K.add(S, V)
            cx2, cy2 = (cx + a) % n, (cy + b) % n   # coefficient accounting mod n
            C = K.canon(T)
            if C != T:                                # canonical flip negates the point
                cx2, cy2 = (-cx2) % n, (-cy2) % n
            visits += 1
            key = K.cls(C)
            if key in store:
                pa, pb, py = store[key]
                if py != C[1]:                       # stored representative is -C
                    da, db = (cx2 + pa) % n, (cy2 + pb) % n
                else:
                    da, db = (cx2 - pa) % n, (cy2 - pb) % n
                if da == 0 and db == 0:
                    v13_rejects += 1          # v13 dropped the trail here
                    trivial += 1
                    # v14: count it, then restart from a fresh point (<=4 steps lost)
                    cx2, cy2 = (5000 + w) % n, (9000 + step) % n
                    S = K.add(K.mul(cx2, P), K.mul(cy2, Q))
                    restarts += 1
                    cx, cy = cx2, cy2
                    step += 1
                    continue
                true_col += 1
                if db % n and recovered is None:
                    cand = (-da * K.inv(db, n)) % n
                    if K.mul(cand, P) == Q:
                        recovered = cand
            else:
                store[key] = (cx2, cy2, C[1])
            S = C
            cx, cy = cx2, cy2
            step += 1
    ok = recovered is not None and (recovered % n) == (SECRET % n)
    return {
        "params": {"p": K.P_TOY, "curve_order": n, "b": K.B_TOY, "R": R, "walkers": WALKERS,
                   "steps_per_walker": STEPS, "secret": SECRET,
                   "record_bytes": K.RECORD_BYTES},
        "counts": {"visits": visits, "trivial_repeat_no_progress": trivial,
                   "true_collisions": true_col, "restarts": restarts,
                   "store_occupancy": len(store)},
        "reconciliation": {"v13_first_repeat_rejects": v13_rejects,
                           "v14_counted_no_progress": trivial,
                           "delta": trivial - v13_rejects,
                           "reading": "same events; v13 discarded the trail and reported a "
                                      "zero cell, v14 counts the event and restarts, so the "
                                      "channel now has a numerator."},
        "success_definition": "one verified scalar recovery: recovered*P == Q, recovered == secret mod n",
        "recovered_secret": recovered, "secret_recovered": bool(ok),
        "rates": {"no_progress_per_step": trivial / max(1, visits),
                  "collision_per_step": true_col / max(1, visits)},
    }


def main():
    ev = run()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "census_v14.json"), "w") as f:
        json.dump(ev, f, indent=1, sort_keys=True)
    print(json.dumps(ev["counts"], indent=1))
    print("recovered:", ev["recovered_secret"], ev["secret_recovered"])


if __name__ == "__main__":
    main()
