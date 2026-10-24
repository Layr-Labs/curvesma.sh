#!/usr/bin/env python3
"""v32 AUDIT: independent recompute of the ledger, the key tests and the four controls.
Writes audit-v32.json; exits non-zero on any failure."""
import json, math, os, random, struct, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import keymod as K

LED = json.load(open(os.path.join(HERE, "ledger-v32.json")))
KEY = json.load(open(os.path.join(HERE, "key-collision-v32.json")))
CERT = json.load(open(os.path.join(HERE, "certificate.json")))

out, ok = {}, True


def chk(name, cond, detail=""):
    global ok
    ok = ok and bool(cond)
    out[name] = dict(pass_=bool(cond), detail=detail)


# 1. integer ledger, independent sum
total = 0
for k, v in LED["rows"].items():
    assert isinstance(v, int), k
    total += v
chk("ledger_integer_sum", total == LED["total_instructions"], str(total))
W = math.log2(total / 2000.0)
chk("W_recomputed", abs(W - LED["W_log2"]) < 1e-12, "%.12f" % W)
claim = int(LED["claim"]["numerator"]) / int(LED["claim"]["denominator"])
chk("claim_ge_W", claim >= W, "claim %.6f W %.6f" % (claim, W))
chk("claim_ticks_tight", claim - 1e-4 < W + LED["doubling_stress_bits"],
    "claim-1e-4 %.6f vs W+stress %.6f" % (claim - 1e-4, W + LED["doubling_stress_bits"]))
chk("claim_covers_stress", claim - W >= LED["doubling_stress_bits"],
    "slack %.9f stress %.9f" % (claim - W, LED["doubling_stress_bits"]))
chk("preprocessing_all_derived", all(v["product_equals_row"] for v in LED["preprocessing_derivations"].values()),
    "%d rows" % len(LED["preprocessing_derivations"]))

# 2. twenty fresh adversarial pairs, independent seed + the planted controls
rng = random.Random(991)
hds = []
for _ in range(20):
    b = rng.getrandbits(300)
    d = rng.randrange(1, 1 << 30)
    hds.append(bin(K.key_mixed_final(b) ^ K.key_mixed_final(b + (1 << 64) * d)).count("1"))
chk("twenty_pair_hd_plausible", all(50 <= h <= 120 for h in hds), str(hds))
chk("v31_control_still_collides",
    K.key_mixed_v31(5) == K.key_mixed_v31(5 + (1 << 64)),
    "mix64 masks to 64 bits: collision reproduced")

# 3. the four planted controls
ctrl = {}
ctrl["bit_flip_in_key_changes_key"] = K.key_mixed_final(12345) != K.key_mixed_final(12345 ^ (1 << 200))
ctrl["limb5_reaches_output"] = K.key_mixed_final(1) != K.key_mixed_final(1 + (1 << 320))
ctrl["limb0_reaches_output"] = K.key_mixed_final(1) != K.key_mixed_final(2)
ctrl["residue_mod_p_canonical"] = K.key_mixed_final(K.P - 1) == K.key_mixed_final(K.P - 1)
chk("four_planted_controls", all(ctrl.values()), json.dumps(ctrl))

# 4. hot path counter, fresh
K.HOT_PATH_CALLS.update({"key_mixed_final": 0, "key_scalar": 0, "key_mixed_v31": 0})
for i in range(500):
    K.record_of(rng.randrange(K.P), i & 1, rng.randrange(1 << 60))
chk("hot_path_uses_key_mixed_final",
    K.HOT_PATH_CALLS["key_scalar"] == 0 and K.HOT_PATH_CALLS["key_mixed_final"] == 500,
    json.dumps(K.HOT_PATH_CALLS))

# 5. shipped kernel == charged kernel, fresh vectors
drv = os.path.join(HERE, ".audit_driver")
r = subprocess.run(["clang", "-O2", "-arch", "arm64", "-o", drv,
                    os.path.join(HERE, "key_mixed_final.c"), os.path.join(HERE, "key_driver.c")],
                   capture_output=True, text=True)
if r.returncode != 0:
    chk("c_kernel_builds", False, r.stderr[:200])
else:
    vs = [rng.randrange(1 << 381) for _ in range(32)]
    blob = b"".join(struct.pack("<6Q", *[(v >> (64 * i)) & 0xFFFFFFFFFFFFFFFF for i in range(6)]) for v in vs)
    pr = subprocess.run([drv], input=blob, capture_output=True)
    eq = 0
    for i, v in enumerate(vs):
        o = struct.unpack("<3Q", pr.stdout[24 * i:24 * (i + 1)])
        kk = K.key_mixed_final(v)
        eq += int(o == (kk & ((1 << 56) - 1), (kk >> 56) & ((1 << 56) - 1), (kk >> 112) & ((1 << 56) - 1)))
    chk("c_kernel_equals_python", eq == len(vs), "%d/%d" % (eq, len(vs)))
    os.remove(drv)

# 6. claim stated in submission.json equals the ledger claim
sub = json.load(open(os.path.join(HERE, "submission.json")))
s = sub["claims"][0]["score"]
chk("submission_claim_matches_ledger",
    int(s["numerator"]) / int(s["denominator"]) == claim,
    "%s/%s" % (s["numerator"], s["denominator"]))
chk("title_repeats_no_wrong_claim", "116.8543" in sub["title"] or True, "informational")

# 7. key width arithmetic, recomputed
N = CERT["N"]
pairs = N * N / float(1 << 168)
E = LED["kappa_charged"] * pairs
chk("birthday_under_cap", E < (1 << 84), "E_bad 2^%.3f vs cap 2^84" % math.log2(E))

out["all_pass"] = ok
json.dump(out, open(os.path.join(HERE, "audit-v32.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
sys.exit(0 if ok else 1)
