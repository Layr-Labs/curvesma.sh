# Argument, v15q

Claim `51313/400` = 128.2825 at success >= 0.39, track secp256k1-ecdlp. This file
states the method and points at the machine-emitted numbers; it introduces no
number of its own.

## What is priced

`evidence/ledger-v15q.json` holds one row per phase, each an exact rational in
protected-normalised instructions (C = 1000). The step is 1853 Ir: the point
listing's head+tail per iteration (744), the alignment symbol plus its two field
products (430), the label symbol plus its two Montgomery products (601), and the
scalar detector (78). The report routine is a separate row (four-probe bound,
one per delivered record).

Each row is a CONSERVATIVE PATH UPPER BOUND, never an equality with an executed
count: an acyclic symbol can only skip instructions through forward branches, so
its static length bounds every path, and a symbol carrying a backward branch is
priced body x trip count plus the remainder with the trip count and its C source
line printed in `evidence/check-executed-v15q.json:symbols`. Branch targets are
recognised in both Darwin `LBB` and ELF `.L` form. The counted translation unit
is built with the vectoriser disabled, so `simd_instructions_in_counted_symbols`
is 0 and no SIMD-to-scalar translation is assumed anywhere.

## What moved the claim

1. The refill row: the priced step updates a prefix/suffix inversion schedule
   (`PROD *= d`, `SUF *= d`), which requires one 256-bit inversion per refill
   window. One inversion (267 Montgomery products) per 67-step window per walker
   is 905 Ir/step, +0.575 bits.
2. The scalar rebuild of the detector and the exact acyclic counts of align and
   label (430 and 601 against 435 and 604 in v15p).

## What the index fix does to the claim's structure

With the slot removed from the address, `idx` is a function of the stored
canonical point, so the transition is a map on the `(r-1)/6` classes and the
partner row can be tested on STORED points: all six representatives (including
the three negatives) must store the same point and the same index, which the
shipped harness reports with zero mismatches while the v15p rule reports 20/20
and 60/60. `slot` survives only as the label's `lambda^slot` factor, which the
canonicalisation forces.

## What is not established

The success side is not closed on the shipped object: the census through the
compiled pool returned more key matches than the class-collision law predicts,
so its recovery path is not trustworthy and `evidence/census-v15q.json:floor_met`
is `false`. The claim is therefore a priced HOLD: the accounting and the map
structure are the deliverable, not the success floor.

## Reproduction

    python3 scripts/build_evidence_v15q.py     # checker, ledger, certificate
    python3 scripts/harness_v15q.py            # invariant + partner + 6 controls
    V15Q_TRIALS=512 python3 scripts/census_v15q.py

No production execution is claimed: the harness runs the real curve for 500-step
walks and the census runs the shipped routines on the q25 fixture.
