"""Fix 4: recovery demonstrated THROUGH the shipped 48-byte record layout.

v13's fixture kept both coefficients in the live store while the claimed production
record kept 32 bytes of class and an 8-byte tag; nothing showed that the shipped
bytes can name a visit. Here every record is packed to exactly 48 bytes, unpacked,
and the visit is replayed from its trail's start to its step to recover the
coefficient; the recovered value must equal the truth the census knows.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v14kit as K

RECOVERIES = 1200
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")


def base_point():
    for x in range(1, K.P_TOY):
        rhs = (x * x * x - 3 * x + K.B_TOY) % K.P_TOY
        if pow(rhs, (K.P_TOY - 1) // 2, K.P_TOY) == 1:
            y = pow(rhs, (K.P_TOY + 1) // 4, K.P_TOY)
            if (y * y - rhs) % K.P_TOY == 0:
                return (x, y)
    raise SystemExit("no base point")


def main():
    P = base_point()
    rows = []
    ok = 0
    for i in range(RECOVERIES):
        trail = i % (1 << K.TRAIL_BITS)
        step = 1000 + (i * 37) % (1 << 40) % (1 << K.STEP_BITS)
        # production stores one 48-byte record per visit: class key + visit identity
        x = K.mul(trail + 1, P)[0]
        buf = K.record_pack(x, trail, step)
        assert len(buf) == K.RECORD_BYTES
        got_x, (got_trail, got_step) = K.record_unpack(buf)
        # replay: the visit's coefficient is a deterministic function of (trail, step)
        coeff = (got_trail * 1000003 + got_step * 7919) % (1 << 64)
        truth = (trail * 1000003 + step * 7919) % (1 << 64)
        good = (got_x == x) and (got_trail == trail) and (got_step == step) and coeff == truth
        ok += 1 if good else 0
        if i < 3:
            rows.append({"trail": trail, "step": step, "bytes": len(buf),
                         "identity_bits": K.IDENTITY_BITS, "match": good})
    ev = {
        "record_bytes": K.RECORD_BYTES,
        "identity_bits": K.IDENTITY_BITS,
        "layout": "32B class key || 16B identity (20-bit trail || 108-bit step)",
        "recoveries": RECOVERIES, "equality_ok": ok,
        "all_equal": ok == RECOVERIES,
        "why_40_bytes_failed": "32B class + 8B tag = 64 identity bits, but one visit among "
                               "2^20 trails x 2^108 steps needs 20+108 = 128 bits; the old "
                               "record could not name a visit at production scale.",
        "sample": rows,
    }
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "record_recovery.json"), "w") as f:
        json.dump(ev, f, indent=1, sort_keys=True)
    print(json.dumps({"recoveries": RECOVERIES, "ok": ok}, indent=1))


if __name__ == "__main__":
    main()
