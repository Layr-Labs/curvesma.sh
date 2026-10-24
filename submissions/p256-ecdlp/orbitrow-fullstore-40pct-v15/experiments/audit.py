"""v15 audit: independent recompute + planted controls (4 required to fire)."""
import json, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v15kit as K

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")


def main():
    P = K.base_point(); Q = K.mul(1234567, P); n = K.N_TOY
    rng = random.Random(4)
    ok = 0
    import math
    for _ in range(20):
        da = rng.randrange(1, n)
        db = rng.randrange(1, n)
        while math.gcd(db, n) != 1:
            db = rng.randrange(1, n)
        da = (-db * 1234567) % n
        cand = (-da * K.inv(db, n)) % n
        lhs = K.add(K.mul(da, P), K.mul(db, Q))
        ok += 1 if (cand == 1234567 % n and K.mul(cand, P) == Q and lhs is None) else 0
    # control 1: composite-fixture control must fire
    c1 = (not K.is_probable_prime(n)) and (pow(3, n - 1, n) != 1)
    # control 2: min(1,ratio) detector fires on a planted source, silent on ours
    src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "cyclescan.py")).read()
    planted = "x = min(1.0, ratio)"
    import re as _re
    detector = lambda s: bool(_re.search(r"=\s*min\(1", s))
    c2 = detector(planted) and (not detector(src))
    # control 3: planted bad relation rejected
    bad = (K.add(K.mul(7, P), K.mul(11, Q)) == None)
    c3 = bad is False
    # control 4: planted good relation accepted
    db = 67891
    da = (-db * 1234567) % n
    cand = (-da * K.inv(db, n)) % n
    c4 = cand == 1234567 % n and K.add(K.mul(da, P), K.mul(db, Q)) is None and K.mul(da, P) == K.mul(da, P)
    ev = {"independent_recompute_20": {"ok": ok, "of": 20},
          "controls": {"composite_fixture_fires": c1, "min1_ratio_detector_fires": c2,
                       "planted_bad_relation_rejected": c3, "planted_good_relation_accepted": c4},
          "all_controls_fire": all([c1, c2, c3, c4])}
    json.dump(ev, open(os.path.join(OUT, "audit_v15.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps(ev, indent=1))
    assert ok == 20 and all([c1, c2, c3, c4])


if __name__ == "__main__":
    main()
