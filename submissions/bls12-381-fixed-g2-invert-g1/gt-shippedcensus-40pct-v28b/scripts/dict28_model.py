#!/usr/bin/env python3
"""dict28_model.py -- differential test of the shipped wide-address dictionary
against an independent Python model, including the >2^64-slot SPARSE path.

The model re-implements the address arithmetic in exact integers (128-bit
masks): h = mix128(key), a = h + salt*(probe+1) mod 2^128, slot = mulhi128(a,cap)
= (a*cap) >> 128.  If the C and the model agree on every status, probe count,
returned 128-bit reference and eviction victim over 20000 census-traffic ops,
the address widening and the eviction policy are both exercised; the sparse test
then drives cap = 2^70+12345 so the 127-bit slot field is load-bearing.
"""
import ctypes
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(os.path.dirname(HERE), "evidence")
IR = os.path.join(HERE, "ir28")
M128 = (1 << 128) - 1
M64 = (1 << 64) - 1
STRIDE_HI = 0x9E3779B97F4A7C15
STRIDE_LO = 0xC2B2AE3D27D4EB4F
PROBES, ENTRY = 4, 9


def lib(name):
    L = ctypes.CDLL(os.path.join(IR, name))
    u64p = ctypes.POINTER(ctypes.c_uint64)
    L.dict28_lookup.argtypes = [u64p, ctypes.c_uint64, ctypes.c_uint64,
                                ctypes.c_uint64, u64p, u64p, u64p]
    L.dict28_insert.argtypes = [u64p, ctypes.c_uint64, ctypes.c_uint64, ctypes.c_uint64,
                                u64p, u64p, u64p, u64p]
    L.dict28_addr_export.argtypes = [u64p, ctypes.c_uint64, ctypes.c_uint64,
                                     ctypes.c_uint64, u64p, u64p]
    return L


class Model:
    def __init__(self, cap, stride_hi=STRIDE_HI):
        self.cap, self.stride_hi = cap, stride_hi
        self.tab = {}

    def mix128(self, key):
        h = (self.stride_hi << 64) | STRIDE_LO
        for w in key:
            h ^= (w + self.stride_hi + ((h << 6) & M128) + (h >> 2)) & M128
            h = (h * 0xFF51AFD7ED558CCD) & M128
            h ^= h >> 67
        h = (h * 0xC4CEB9FE1A85EC53) & M128
        h ^= h >> 71
        return h & M128

    def addr(self, key, probe):
        salt = (self.stride_hi << 64) | STRIDE_LO
        a = (self.mix128(key) + salt * (probe + 1)) & M128
        return (a * self.cap) >> 128

    def lookup(self, key):
        for i in range(PROBES):
            s = self.addr(key, i)
            e = self.tab.get(s)
            if e is None:
                return 0, i + 1, None
            if e[0] == tuple(key):
                return 1, i + 1, e[1]
        return 0, PROBES, None

    def insert(self, key, ref, evict=True):
        slots = [self.addr(key, i) for i in range(PROBES)]
        for i, s in enumerate(slots):
            e = self.tab.get(s)
            if e is None:
                self.tab[s] = (tuple(key), ref)
                return 0, None
            if e[0] == tuple(key):
                self.tab[s] = (tuple(key), ref)
                return 0, None
        if not evict:
            return 2, None
        best = 0
        for i in range(1, PROBES):
            if self.tab[slots[i]][0][0] < self.tab[slots[best]][0][0]:
                best = i
        victim = self.tab[slots[best]][1]
        self.tab[slots[best]] = (tuple(key), ref)
        return 1, victim


def main():
    out = {"mode": "dict28-differential"}
    # ---- differential, 20000 ops, evicting build, cap 4096 (load factor 0.625)
    LE = lib("dict28.dylib")
    LN = lib("dict28_nolru.dylib")
    rng = random.Random(0xD28)
    cap = 4096
    tab = (ctypes.c_uint64 * (cap * ENTRY))()
    keys = [[rng.getrandbits(64) for _ in range(6)] for _ in range(8192)]
    model = Model(cap)
    kb = (ctypes.c_uint64 * 6)()
    ref = (ctypes.c_uint64 * 2)()
    ev = (ctypes.c_uint64 * 2)()
    np_ = ctypes.c_uint64(0)
    mism = {"status": 0, "probes": 0, "ref": 0, "victim": 0, "key": 0}
    evictions = inserts = lookups = hits = misses = 0
    digest_c = 0
    for op in range(20000):
        key = keys[rng.randrange(8192)]
        for i in range(6):
            kb[i] = key[i]
        if op % 2:
            r = rng.getrandbits(64)
            ref[0], ref[1] = r, r >> 1
            sc = LE.dict28_insert(tab, cap, 0, cap, kb, ref, ev, ctypes.byref(np_))
            sm, vic = model.insert(key, (r, r >> 1))
            if sc != sm:
                mism["status"] += 1
            if sc == 1 and (ev[0], ev[1]) != vic:
                mism["victim"] += 1
            evictions += sc == 1
            inserts += 1
            digest_c = (digest_c * 1000003 + sc + (ev[0] & 0xFF)) & M64
        else:
            sc = LE.dict28_lookup(tab, cap, 0, cap, kb, ref, ctypes.byref(np_))
            sm, pr, rf = model.lookup(key)
            if sc != sm:
                mism["status"] += 1
            if np_.value != pr:
                mism["probes"] += 1
            if sc == 1 and (ref[0], ref[1]) != rf:
                mism["ref"] += 1
            hits += sc == 1
            misses += sc == 0
            lookups += 1
            digest_c = (digest_c * 1000003 + sc + (ref[0] & 0xFF)) & M64
    # ---- no-evict control differential (2 = drop)
    tab2 = (ctypes.c_uint64 * (cap * ENTRY))()
    m2 = Model(cap)
    rng2 = random.Random(0xD29)
    drops = 0
    for op in range(8000):
        key = keys[rng2.randrange(8192)]
        for i in range(6):
            kb[i] = key[i]
        r = rng2.getrandbits(64)
        ref[0], ref[1] = r, r >> 1
        sc = LN.dict28_insert(tab2, cap, 0, cap, kb, ref, ev, ctypes.byref(np_))
        sm, _ = m2.insert(key, (r, r >> 1), evict=False)
        if sc != sm:
            mism["status"] += 1
        drops += sc == 2

    # ---- WIDE SPARSE path: cap = 2^70 + 12345, slots above 2^64 must route
    capw = (1 << 70) + 12345
    above = 0
    trunc_ambiguous = 0
    n = 200000
    mw = Model(capw)
    rng3 = random.Random(0xD2A)
    slo = ctypes.c_uint64(0)
    shi = ctypes.c_uint64(0)
    agree = 0
    for _ in range(n):
        key = [rng3.getrandbits(64) for _ in range(6)]
        for i in range(6):
            kb[i] = key[i]
        for probe in range(PROBES):
            LE.dict28_addr_export(kb, probe, capw & M64, capw >> 64,
                                  ctypes.byref(slo), ctypes.byref(shi))
            cslot = (shi.value << 64) | slo.value
            mslot = mw.addr(key, probe)
            if cslot == mslot:
                agree += 1
            if cslot >= (1 << 64):
                above += 1
            if cslot >= capw:
                trunc_ambiguous += 1
    # a slot above 2^64 has no 64-bit representation: truncating it changes the slot
    trunc_rows = 0
    for key in keys[:2000]:
        for i in range(6):
            kb[i] = key[i]
        LE.dict28_addr_export(kb, 0, capw & M64, capw >> 64, ctypes.byref(slo), ctypes.byref(shi))
        cslot = (shi.value << 64) | slo.value
        if cslot >> 64 and (cslot & M64) != cslot:
            trunc_rows += 1

    out.update({
        "differential": {"ops": 20000, "lookups": lookups, "inserts": inserts,
                         "hits": hits, "misses": misses, "evictions": evictions,
                         "mismatches": mism, "digest": digest_c,
                         "control_ops": 8000, "control_drops": drops},
        "wide_path": {"cap": str(capw), "keys": n, "probe_pairs": n * PROBES,
                      "c_model_agreement": agree, "slots_above_2_64": above,
                      "slots_above_cap": trunc_ambiguous,
                      "truncation_changes_slot_rows": trunc_rows,
                      "note": "a 64-bit slot field cannot represent these slots, so the 2-word "
                              "reference and the 128-bit address are load-bearing"},
        "zero_mismatches": all(v == 0 for v in mism.values()),
        "mutation_control": None,
    })
    # ---- mutation control: change the stride in the model, mismatches must appear
    mm = {"probe": 0}
    model3 = Model(cap, stride_hi=STRIDE_HI + 1)
    tab3 = (ctypes.c_uint64 * (cap * ENTRY))()
    rng4 = random.Random(0xD2B)
    diffs = 0
    for op in range(4000):
        key = keys[rng4.randrange(8192)]
        for i in range(6):
            kb[i] = key[i]
        r = rng4.getrandbits(64)
        ref[0], ref[1] = r, r >> 1
        LE.dict28_insert(tab3, cap, 0, cap, kb, ref, ev, ctypes.byref(np_))
        model3.insert(key, (r, r >> 1))
    for op in range(4000):
        key = keys[rng4.randrange(8192)]
        for i in range(6):
            kb[i] = key[i]
        sc = LE.dict28_lookup(tab3, cap, 0, cap, kb, ref, ctypes.byref(np_))
        sm, pr, rf = model3.lookup(key)
        if sc != sm or np_.value != pr or (sc == 1 and (ref[0], ref[1]) != rf):
            diffs += 1
    out["mutation_control"] = {"changed": "model mix128 stride 0x9E37...C15 -> ...C16",
                               "differential_rows_changed": diffs,
                               "detected": diffs > 0}
    json.dump(out, open(os.path.join(EV, "dict28.json"), "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("differential", "wide_path", "zero_mismatches",
                                          "mutation_control")}, indent=1)[:1600])


if __name__ == "__main__":
    main()
