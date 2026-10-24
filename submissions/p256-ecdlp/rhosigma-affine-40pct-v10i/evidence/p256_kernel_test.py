#!/usr/bin/env python3
"""p256_kernel_test.py — verify the SHIPPED p256 kernel against Python bigint.

Builds p256_kernel_test.c twice:
  * fixed   — against evidence/p256_fe_cc.h (the shipped, limb-repaired header)
  * control — against evidence/p256_fe_cc_broken.h, renamed into a scratch dir so
              the include actually resolves to the v10f (transposed) header.
Every emitted row is recomputed with exact Python integers at the real P-256
modulus.  PASS requires fixed: 0 mismatches AND control: >0 mismatches; a control
that passes means the test cannot detect the defect and the run is IGNORED.

usage: python3 p256_kernel_test.py [out.json]
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
P = (1 << 256) - (1 << 224) + (1 << 192) + (1 << 96) - 1
R = 1 << 256                 # Montgomery radix of the CIOS kernel
RI = pow(R, -1, P)          # fe_mul computes a*b*RI mod p (Montgomery product)
assert P == 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
NRAND = 10000
NEDGE = 15


def build(outbin, include_dir, srcdir):
    src = os.path.join(srcdir, "p256_kernel_test.c")
    cmd = ["clang", "-O2", "-w", "-I", include_dir, src, "-o", outbin]
    subprocess.run(cmd, check=True, capture_output=True)
    return cmd


def prepare_control(tmp):
    for fn in ("p256_kernel_test.c", "step_v10f_o3.c"):
        shutil.copy(os.path.join(HERE, fn), os.path.join(tmp, fn))
    shutil.copy(os.path.join(HERE, "p256_fe_cc_broken.h"),
                os.path.join(tmp, "p256_fe_cc.h"))


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def limits():
    p = P
    return [
        0, 1, 2, p - 1, p - 2, (p - 1) // 2, 1 << 63, 1 << 64, 1 << 128,
        1 << 192, 1 << 224, 1 << 255, (1 << 255) + 1, (1 << 256) - (1 << 224),
        p - 1 - (1 << 96),
    ]


def canon_y(y):
    return y if (y & 1) == 0 else (P - y) % P


def check(text):
    """Recompute every record; return (counts, first mismatches per op)."""
    n = 0
    nv = 0
    bad = {"mul": 0, "add": 0, "sub_ab": 0, "sub_ba": 0, "inv": 0, "canon": 0}
    noncanon = 0
    first = {}
    seen_a = set()
    for line in text.splitlines():
        f = line.split()
        if not f:
            continue
        tag = f[0]
        v = [int(x) for x in f[1:]]
        n += 1
        if tag == "V":
            nv += 1
            a = sum(v[i] << (64 * i) for i in range(4))
            b = sum(v[4 + i] << (64 * i) for i in range(4))
            seen_a.add(a)
            got = {}
            for k, name in enumerate(("mul", "add", "sub_ab", "sub_ba")):
                got[name] = sum(v[8 + 16 * 0 + 4 * k + i] << (64 * i) for i in range(4))
            exp = {"mul": (a * b * RI) % P, "add": (a + b) % P,
                   "sub_ab": (a - b) % P, "sub_ba": (b - a) % P}
            for name in exp:
                if got[name] != exp[name]:
                    bad[name] += 1
                    first.setdefault(name, {"a": hex(a), "b": hex(b),
                                            "got": hex(got[name]), "want": hex(exp[name])})
                if got[name] >= P:
                    noncanon += 1
        elif tag == "I":
            a = sum(v[i] << (64 * i) for i in range(4))
            got = sum(v[4 + i] << (64 * i) for i in range(4))
            # law: for a != 0, MontMul(a, fe_inv(a)) == R (Montgomery one);
            # a == 0 has no inverse and the routine returns 0.
            exp = None
            ok = (got == 0) if a == 0 else (((a * got) % P * RI % P) == R % P)
            if not ok:
                exp = (R * R * pow(a, -1, P)) % P if a else 0
                bad["inv"] += 1
                first.setdefault("inv", {"a": hex(a), "got": hex(got),
                                         "want_law": "MontMul(a,inv)==R",
                                         "want_mont_inverse": hex(exp)})
        elif tag == "C":
            y = sum(v[i] << (64 * i) for i in range(4))
            got = sum(v[4 + i] << (64 * i) for i in range(4))
            if got != canon_y(y):
                bad["canon"] += 1
                first.setdefault("canon", {"y": hex(y), "got": hex(got),
                                           "want": hex(canon_y(y))})
        else:
            raise SystemExit("unknown record %r" % tag)
    edges_seen = set(limits()) & seen_a
    return {"records": n, "mismatches": bad, "total_mismatches": sum(bad.values()),
            "noncanonical_outputs": noncanon,
            "first_mismatch": first, "distinct_a_values": len(seen_a),
            "edge_values_seen": len(edges_seen), "vectors": nv}


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "kernel_test.json")
    tmp = tempfile.mkdtemp(prefix="p256kt_")
    try:
        fixed_bin = os.path.join(tmp, "kt_fixed")
        cmd_fixed = build(fixed_bin, HERE, HERE)
        ctl_dir = os.path.join(tmp, "ctl")
        os.makedirs(ctl_dir)
        prepare_control(ctl_dir)
        ctl_bin = os.path.join(tmp, "kt_control")
        cmd_ctl = build(ctl_bin, ctl_dir, ctl_dir)
        fixed_txt = subprocess.run([fixed_bin], check=True,
                                   capture_output=True, text=True).stdout
        ctl_txt = subprocess.run([ctl_bin], check=True,
                                 capture_output=True, text=True).stdout
    finally:
        pass
    fixed = check(fixed_txt)
    ctl = check(ctl_txt)
    edges = limits()
    payload = {
        "modulus": hex(P),
        "mont_radix": hex(R),
        "mul_law": "fe_mul(a,b) = a*b*R^-1 mod p (CIOS Montgomery product)",
        "inv_law": "fe_inv(a) = R^2*a^-1 mod p; fe_mul(a, fe_inv(a)) = R",
        "random_vectors": NRAND,
        "edge_values": [hex(e) for e in edges],
        "fixed": fixed,
        "control": ctl,
        "compile_fixed": " ".join(cmd_fixed),
        "compile_control": " ".join(cmd_ctl),
        "sha256_header_fixed": sha256(os.path.join(HERE, "p256_fe_cc.h")),
        "sha256_header_broken": sha256(os.path.join(HERE, "p256_fe_cc_broken.h")),
        "checks": {},
    }
    ok_free = fixed["total_mismatches"] == 0
    ok_ctl = ctl["total_mismatches"] > 0
    ok_rec = fixed["records"] >= 10000
    ok_edge = fixed["edge_values_seen"] == NEDGE
    ok_ctl_compiled = ctl["records"] == fixed["records"]
    payload["checks"] = {
        "fixed_kernel_matches_bigint": "PASS" if ok_free else "FAIL",
        "mutation_control_fails": "PASS" if ok_ctl else "FAIL",
        "vector_count_ge_10000": "PASS" if ok_rec else "FAIL",
        "all_edge_values_exercised": "PASS" if ok_edge else "FAIL",
        "control_emitted_same_records": "PASS" if ok_ctl_compiled else "FAIL",
    }
    shutil.rmtree(tmp, ignore_errors=True)
    payload["verdict"] = "PASS" if all(v == "PASS" for v in payload["checks"].values()) else "FAIL"
    json.dump(payload, open(out, "w"), indent=1, sort_keys=True)
    print(json.dumps({k: payload[k] for k in ("verdict", "checks")}, indent=1))
    print("fixed mismatches:", fixed["total_mismatches"], "control mismatches:", ctl["total_mismatches"])
    print("records:", fixed["records"], "edges seen:", fixed["edge_values_seen"])
    return 0 if payload["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
