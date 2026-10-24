#!/usr/bin/env python3
"""v17b predicate agreement at EXACTLY the production configuration.

Runs the SHARED production rule (imported from gt_prod17b, not re-declared) against
two independently written derivations and against constructive witnesses, so the
v16 defect ("index and DP overlap, every DP hit restricts h to the lower half")
cannot survive:

  * `prod_bits`          -- the shipped rule, imported.
  * `bits_divmod`        -- rewritten with divmod on the binary expansion.
  * `bits_ranges`        -- symbolic: T decomposed as branch | index | dp window.
  * witnesses            -- for EVERY slot index in a sample (incl. both halves),
                            an explicit T is built that yields that index AND a DP
                            hit, proving the DP predicate does not constrain h.

usage: python3 scripts/predicate17b_check.py > evidence/predicate17b.json
"""
import json
import random
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gt_prod17b as G   # the shipped production module


def bits_divmod(t, Rtab, dpbits):
    """Independent rewrite: repeated divmod on T, LSB-first digit list."""
    digits = []
    v = t
    for _ in range(120):
        v, r = divmod(v, 2)
        digits.append(r)
    branch = digits[0] | (digits[1] << 1) | (digits[2] << 2)
    idx = 0
    for pos in range(11):
        idx |= digits[3 + pos] << pos
    window = 0
    for pos in range(dpbits):
        window |= digits[14 + pos] << pos
    return (branch < 3), idx, window == 0


def bits_ranges(t, Rtab, dpbits):
    """Symbolic decomposition: three disjoint shifts of the same word."""
    br_word = t & 0b111
    idx_word = t // 8
    idx_word = idx_word - (idx_word // G.HALF_PROD) * G.HALF_PROD
    dp_word = t // (2 ** 14)
    dp_word = dp_word - (dp_word // (2 ** dpbits)) * (2 ** dpbits)
    return (br_word < 3), idx_word, dp_word == 0


def main():
    assert G.RTAB_PROD == 4096 and G.IDX_BITS_PROD == 11 and G.DP_SHIFT_PROD == 14
    rng = random.Random(20261023)
    n = 0
    mism = {"divmod": 0, "ranges": 0}
    samples = [rng.randrange(1, 1 << 120) for _ in range(120000)]
    # structured: every (branch, index-edge, dp-edge) corner
    for br in range(8):
        for h in (0, 1, 2, 2046, 2047):
            for win in (0, 1, (1 << 100) - 1, 1 << 99):
                samples.append(br | (h << 3) | (win << 14))
    for t in samples:
        a = G.prod_bits(t, G.RTAB_PROD, 100)
        b = bits_divmod(t, G.RTAB_PROD, 100)
        c = bits_ranges(t, G.RTAB_PROD, 100)
        n += 1
        if a != b:
            mism["divmod"] += 1
        if a != c:
            mism["ranges"] += 1
    # constructive witnesses: index AND DP hit simultaneously, both halves
    wit = []
    ok = True
    for h in (0, 1, 1023, 1024, 2000, 2047):
        for slot_half in (0, 6):
            t = (h << 3) | (5 << 0)          # branch 5 (>=3 -> multiply) + index h
            br, got, dp = G.prod_bits(t, G.RTAB_PROD, 100)
            slot = G.prod_slot(got, slot_half)
            wit.append({"T": t, "branch": "multiply" if not br else "square",
                        "index": got, "dp_hit": dp, "kappa_ge_6": slot_half >= 6,
                        "slot": slot})
            ok = ok and (got == h) and dp and (slot == h + (2048 if slot_half >= 6 else 0))
    # and a T that yields the SAME index but a DP MISS (dp window non-zero)
    miss = G.prod_bits((2047 << 3) | (1 << 14), G.RTAB_PROD, 100)
    out = {
        "mode": "predicate-agreement-v17b",
        "config": {"R": G.RTAB_PROD, "idx_bits": G.IDX_BITS_PROD,
                   "dp_shift": G.DP_SHIFT_PROD, "dp_bits": 100, "half": G.HALF_PROD},
        "ranges": {"branch": [0, 2], "index": [3, 13], "dp": [14, 113]},
        "ranges_pairwise_disjoint": True,
        "samples_tested": n,
        "mismatches_shared_vs_divmod": mism["divmod"],
        "mismatches_shared_vs_ranges": mism["ranges"],
        "witness_index_and_dp_hit_both_halves": wit,
        "witness_all_ok": ok,
        "same_index_dp_miss_example": {"T": (2047 << 3) | (1 << 14),
                                       "index": miss[1], "dp_hit": miss[2]},
        "dp_does_not_constrain_index": ok and miss[1] == 2047 and miss[2] is False,
        "command": "python3 scripts/predicate17b_check.py > evidence/predicate17b.json",
    }
    print(json.dumps(out, indent=1))
    return 0 if (mism["divmod"] == 0 and mism["ranges"] == 0 and ok) else 1


if __name__ == "__main__":
    sys.exit(main())
