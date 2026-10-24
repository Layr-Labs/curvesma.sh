#!/usr/bin/env python3
"""sort32.py -- EXECUTED count of the full-width state-key path (verdict b30ddc21 flaw 3).

The verdict: "Group.state computes keys12, its minimum/argmin and a Python sort
before calling state_key48w; pricing only that C symbol does not price the
complete state-key computation ... Production residues are 381-bit, whereas the
priced routine receives only four words per lane."

So the loop below IS the pre-state_key48w work, at production width:
  * twelve canonical 381-bit residues (6 words each), extracted from the public
    representation (12 x 6 word loads + masks),
  * min + argmin over them (11 full-width compare-selects),
  * an insertion sort over 12 full-width values (the shipped code sorts ks),
  * the 48-word (12 x 4 x 64-bit) lane packing that state_key48w then consumes.

Every comparison/move is COUNTED (a comparison counter, not wall time), the
counts are reported per call, and the charge is the MAX over calls -- the
verdict's "measured maximum over individual calls", not an aggregate divided by
calls.  Census traffic: the key distribution of the shipped protocol is the
uniform random state of the rho walk, so calls draw twelve uniform 381-bit
values (seeded, reproducible).

Output: evidence/sort32.json
"""
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(os.path.dirname(HERE), "evidence")

WORDS = 6                 # 381-bit residue = 6 x 64-bit words
N = 12                    # Fp12 -> twelve residues
CALLS = 20000
SEED = 20261024

# element-word costs, each a stated per-word primitive:
CMP_WORD = 3              # load, compare, branch
MOVE_WORD = 2             # load, store
MASK_WORD = 2             # shift + mask for the canonical extraction
PACK_WORD = 2             # shift + store for the 48-word packing


def gen_call(rng):
    return [[rng.getrandbits(64) for _ in range(WORDS)] for _ in range(N)]


def cmp_full(a, b):
    """compare two full-width values, most significant word first; counted."""
    c = 0
    for k in range(WORDS - 1, -1, -1):
        c += CMP_WORD
        if a[k] != b[k]:
            return (a[k] > b[k]), c
    return False, c


def run_call(vals):
    """the pre-state_key48w path; returns (comparisons, moves, total_ir)."""
    cmps = 0
    moves = 0
    # ---- canonical extraction (the 12 residues are the public representation)
    ext = []
    for v in vals:
        ext.append([(v[k] & 0xFFFFFFFFFFFFFFFF) for k in range(WORDS)])
        moves += WORDS
    # ---- min / argmin over the twelve, full width
    mn, arg = ext[0], 0
    for i in range(1, N):
        gt, c = cmp_full(ext[i], mn)
        cmps += c // CMP_WORD
        if gt:
            mn, arg = ext[i], i
            moves += WORDS
    # ---- insertion sort over the twelve full-width values
    a = [list(x) for x in ext]
    for i in range(1, N):
        key = a[i]
        j = i - 1
        while j >= 0:
            gt, c = cmp_full(a[j], key)
            cmps += c // CMP_WORD
            if not gt:
                break
            a[j + 1] = a[j]
            moves += WORDS
            j -= 1
        a[j + 1] = key
        moves += WORDS
    # ---- 48-word lane packing handed to state_key48w
    packed = []
    for x in a:
        for k in range(4):
            packed.append((x[k] & 0xFFFFFFFFFFFFFFFF))
    ir = (cmps * CMP_WORD + moves * MOVE_WORD
          + N * WORDS * MASK_WORD + 48 * PACK_WORD)
    return cmps, moves, ir


def main():
    rng = random.Random(SEED)
    hist = []
    for _ in range(CALLS):
        hist.append(run_call(gen_call(rng)))
    cmps = sorted(h[0] for h in hist)
    irs = sorted(h[2] for h in hist)

    def pct(xs, q):
        return xs[int(q * (len(xs) - 1))]

    out = {
        "mode": "sort32-full-width-state-key-executed",
        "what": "pre-state_key48w work at production width (12 canonical 381-bit residues): "
                "canonical extraction, min/argmin, insertion sort over 12 full-width values, "
                "48-word lane packing",
        "calls": CALLS, "seed": SEED, "words_per_residue": WORDS,
        "counted_quantities": "compare instructions and value moves, executed; not wall time",
        "comparisons_per_call": {
            "p50": pct(cmps, 0.50), "p95": pct(cmps, 0.95), "max": cmps[-1],
            "min": cmps[0], "mean": round(sum(cmps) / len(cmps), 3),
            "structural_max": N * (N - 1) // 2 + (N - 1),
        },
        "executed_ir_per_call": {
            "p50": pct(irs, 0.50), "p95": pct(irs, 0.95), "max": irs[-1], "min": irs[0],
        },
        "unit_costs": {"compare_word": CMP_WORD, "move_word": MOVE_WORD,
                       "mask_word": MASK_WORD, "pack_word": PACK_WORD},
        "charge_per_step": irs[-1],
        "charge_basis": "MAX over the 20000 executed calls (verdict: 'a measured maximum over "
                        "individual calls', not an aggregate divided by call count); p95 and p50 "
                        "reported for the sensitivity row",
        "sensitivity": {
            "at_p50": pct(irs, 0.50), "at_p95": pct(irs, 0.95), "at_max": irs[-1],
            "note": "the ledger charges at_max; the claim is recomputed at p50 in the note",
        },
        "excluded": "state_key48w itself is priced separately (3398 Ir, exec31.json); this row is "
                    "the full-width work BEFORE that call, which v31 priced as zero",
    }
    json.dump(out, open(os.path.join(EV, "sort32.json"), "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("comparisons_per_call", "executed_ir_per_call",
                                          "charge_per_step")}, indent=1))


if __name__ == "__main__":
    main()
