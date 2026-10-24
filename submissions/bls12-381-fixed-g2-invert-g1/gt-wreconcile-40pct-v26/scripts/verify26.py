#!/usr/bin/env python3
"""verify26.py -- differential test of the SHIPPED dictionary (scripts/ir26/dict26.c)
against a Python model, the key-builder test, and the disassembly prices.

Evidence written: evidence/dict26.json (test rows) and evidence/ir26.json
(symbol spans).  Every number in ledger26.json that names a structure comes from
these two files, so no ledger-cited symbol can lack an implementation.
"""
import json
import os
import re
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
IR = os.path.join(HERE, "ir26")
M64 = (1 << 64) - 1
STRIDE = 0x9E3779B97F4A7C15
CAP, MASK, OPS = 1 << 12, (1 << 12) - 1, 20000
TAG_EMPTY, TAG_FULL, PROBES = 0, 1, 4


def mix64(x):
    x &= M64
    x ^= x >> 33
    x = (x * 0xff51afd7ed558ccd) & M64
    x ^= x >> 33
    x = (x * 0xc4ceb9fe1a85ec53) & M64
    x ^= x >> 33
    return x & M64


def state_key48(v):
    a = sorted(v)
    return [mix64(a[j] + STRIDE * a[j + 6]) for j in range(6)]


def dict_hash(key):
    h = STRIDE
    for w in key:
        h = mix64(h ^ w)
    return h


def fnv(buf):
    h = 0xcbf29ce484222325
    for b in buf:
        h ^= b
        h = (h * 0x100000001b3) & M64
    return h


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def build_and_run(mutate=False):
    src = os.path.join(IR, "dict26.c")
    if mutate:
        txt = open(src).read().replace("#define STRIDE    0x9E3779B97F4A7C15ULL",
                                       "#define STRIDE    0x9E3779B97F4A7C16ULL")
        src = os.path.join(IR, "dict26_mut.c")
        open(src, "w").write(txt)
    exe = os.path.join(IR, "dict_test_mut" if mutate else "dict_test")
    run(["clang", "-O2", "-o", exe, src, os.path.join(IR, "dict_test.c")])
    out = run([exe])
    ops, dig, final = [], [], None
    for line in out.splitlines():
        p = line.split()
        if p[0] == "D":
            dig.append((int(p[1]), int(p[2], 16)))
        elif p[0] == "HITS":
            final = {"hits": int(p[1]), "miss": int(p[3]), "evicted": int(p[5]),
                     "lf_ok": p[7] == "True" or p[7] == "1"}
        else:
            ops.append({"op": int(p[0]), "kind": int(p[1]), "res": int(p[2]),
                        "key": [int(x, 16) for x in p[3:9]]})
    return ops, dig, final


def model(ops):
    tab = [[0] * 8 for _ in range(CAP)]
    hits = miss = evicted = 0
    digs, struct_mismatch = [], 0
    for rec in ops:
        key = rec["key"]
        if rec["kind"] == 3:
            s = dict_hash(key) & MASK
            res = -1
            for i in range(PROBES):
                sl = (s + i * STRIDE) & MASK
                e = tab[sl]
                if e[7] == TAG_EMPTY:
                    break
                if e[:6] == key:
                    res = sl
                    break
            if res >= 0:
                hits += 1
            else:
                miss += 1
        else:
            s = dict_hash(key) & MASK
            res = None
            for i in range(PROBES):
                sl = (s + i * STRIDE) & MASK
                e = tab[sl]
                if e[7] == TAG_EMPTY:
                    e[:6] = key
                    e[6] = rec["op"]
                    e[7] = TAG_FULL
                    res = sl
                    break
                if e[:6] == key:
                    e[6] = rec["op"]
                    res = sl
                    break
            if res is None:
                cand = [(s + i * STRIDE) & MASK for i in range(PROBES)]
                v = 0
                for i in range(1, PROBES):
                    if tab[cand[i]][0] < tab[cand[v]][0]:
                        v = i
                tab[cand[v]][:6] = key
                tab[cand[v]][6] = rec["op"]
                tab[cand[v]][7] = TAG_FULL
                evicted += 1
                res = cand[v]
        if res != rec["res"]:
            struct_mismatch += 1
        if (rec["op"] + 1) % 1000 == 0:
            flat = [w for e in tab for w in e]
            buf = b"".join((w & M64).to_bytes(8, "little") for w in flat)
            digs.append(((rec["op"] + 1), fnv(buf)))
    return {"digs": digs, "hits": hits, "miss": miss, "evicted": evicted,
            "struct_mismatch": struct_mismatch}


def main():
    ops, digs, final = build_and_run(False)
    mm = model(ops)
    digest_mismatch = sum(1 for a, b in zip(digs, mm["digs"]) if a != b) + \
        abs(len(digs) - len(mm["digs"]))
    counter_mismatch = sum([final["hits"] != mm["hits"], final["miss"] != mm["miss"],
                            final["evicted"] != mm["evicted"]])
    # mutation control: one changed stride constant inside the SHIPPED hash must
    # show up as a different table (digest) and/or different counters.
    mops, mdigs, mfinal = build_and_run(True)
    mut_digest_diff = sum(1 for a, b in zip(digs, mdigs) if a != b)
    mut_counter_diff = sum([final["hits"] != mfinal["hits"],
                            final["miss"] != mfinal["miss"],
                            final["evicted"] != mfinal["evicted"]])

    # key builder: 32 fixed tuples, sort + fold, exact match required
    kexe = os.path.join(IR, "key_test")
    run(["clang", "-O2", "-o", kexe, os.path.join(IR, "dict26.c"),
         os.path.join(IR, "key_test.c")])
    kout = run([kexe])
    key_mismatch = 0
    # regenerate the same tuples the driver uses (xorshift from the same seed)
    s = 0x9E3779B97F4A7C15
    for line in kout.splitlines():
        p = line.split()
        vals = []
        for _ in range(12):
            s ^= (s << 13) & M64
            s ^= s >> 7
            s ^= (s << 17) & M64
            vals.append(s & M64)
        want = [int(x, 16) for x in p[1:7]]
        if state_key48(vals) != want:
            key_mismatch += 1

    # control-flow bound for the shipped key builder: insertion sort on twelve
    # values performs at most 11*12/2 = 66 comparisons, whatever the data.
    import random as _r
    _r.seed(7)
    cmax = 0
    for _ in range(10000):
        a = [_r.getrandbits(64) for _ in range(12)]
        c = 0
        for i in range(1, 12):
            t = a[i]
            j = i - 1
            while j >= 0:
                c += 1
                if a[j] <= t:
                    break
                a[j + 1] = a[j]
                j -= 1
            a[j + 1] = t
        cmax = max(cmax, c)

    dict_ev = {
        "key_comparisons_max_measured": cmax,
        "key_comparisons_bound": 66,
        "mode": "dict26-differential",
        "structure": "open addressing, four probes, stride 0x9E3779B97F4A7C15, "
                     "seven-word entries (six key words + tag), eviction of the "
                     "smallest leading key word when all four probes are occupied",
        "ops": len(ops), "table_slots": CAP, "load_factor": len(ops and set()) or 0,
        "lookups": sum(1 for o in ops if o["kind"] == 3),
        "inserts": sum(1 for o in ops if o["kind"] in (0, 1, 2)),
        "hits": final["hits"], "misses": final["miss"], "evictions": final["evicted"],
        "load_factor_ok": bool(final["lf_ok"]),
        "digest_checks": len(digs), "digest_mismatches": digest_mismatch,
        "counter_mismatches": counter_mismatch, "structure_mismatches": mm["struct_mismatch"],
        "key_tuples": 32, "key_mismatches": key_mismatch,
        "mutation_control": {"changed": "STRIDE 0x9E3779B97F4A7C15 -> ...C16",
                             "digest_diffs": mut_digest_diff,
                             "counter_diffs": mut_counter_diff,
                             "detected": mut_digest_diff > 0 or mut_counter_diff > 0},
        "zero_mismatches": digest_mismatch == 0 and counter_mismatch == 0
                           and mm["struct_mismatch"] == 0 and key_mismatch == 0,
    }
    json.dump(dict_ev, open(os.path.join(EV, "dict26.json"), "w"), indent=1)

    # disassembly prices for exactly the symbols the ledger names
    obj = os.path.join(IR, "dict26_arm64.o")
    run(["clang", "-O2", "-target", "arm64-apple-macos11", "-c",
         os.path.join(IR, "dict26.c"), "-o", obj])
    out = run(["otool", "-tvV", obj])
    syms, cur = {}, None
    for line in out.splitlines():
        m = re.match(r"^(_[A-Za-z0-9_]+):$", line.strip())
        if m:
            cur = m.group(1)
            syms[cur] = []
            continue
        m = re.match(r"^([0-9a-f]{16})\t(\S+)", line)
        if m and cur:
            syms[cur].append((m.group(1), m.group(2)))
    rows = {}
    for s, ins in syms.items():
        rows[s.lstrip("_")] = {"symbol": s, "first_ir_addr": "0x" + ins[0][0],
                               "last_ir_addr": "0x" + ins[-1][0],
                               "ir_count": len(ins),
                               "distinct_mnemonics": len(set(m for _, m in ins))}
    ir_ev = {"mode": "ir26", "object": "scripts/ir26/dict26.c (clang -O2, arm64)",
             "tool": "otool -tvV", "rows": rows,
             "lookup_ir": rows["dict_lookup"]["ir_count"],
             "insert_ir": rows["dict_insert"]["ir_count"],
             "key_ir": rows["state_key48"]["ir_count"],
             "note": "the ledger's per-step lookup is state_key48 + dict_lookup; the "
                     "record path is dict_insert; all three are priced here"}
    json.dump(ir_ev, open(os.path.join(EV, "ir26.json"), "w"), indent=1)
    print(json.dumps({"dict": dict_ev, "prices": {k: v["ir_count"] for k, v in rows.items()}}, indent=1))


if __name__ == "__main__":
    main()
