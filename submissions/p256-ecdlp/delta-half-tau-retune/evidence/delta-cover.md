# Additive debit halved: delta = 2^-21 covering record

This revision halves the additive collision debit from 2^-20 to 2^-21
(delta = 4.76837158203125e-7) and does not import the predecessor's
§2.6 (its bundle omits that section). The debit's covered failures are
identified and tabulated here with their bound references and
evaluated magnitudes at this S (tau = 4943017/10000000).

## Covered failures

| failure | bound reference | magnitude at this S |
|---|---|---|
| finite-horizon truncation of the birthday sum below T steps | §2.2 birthday truncation | exponentially small in T (below 1e-30) |
| stopped-process padding after the first G violation | §2.3 exposure process, padding | exponentially small in S (below 1e-30) |
| merge-resolution residue below the DELTA_MERGE bound | §2.8 Lemma M tail | at most DELTA_MERGE = 1.2520e-14 |
| good-event violations (within-trail coincidences, long closures) | §2.2, DELTA_CUT | DELTA_CUT = 5.8542e-17 |

Covered sum: below 4.728e-14 (dominated by the merge residue; the
truncation/padding terms are exponentially small and contribute
nothing at this precision). The halved 2^-21 debit exceeds the
covered sum by a factor above 1e7.

## Why halving is the scored move

The old 2^-20 debit exceeded the covered sum by a factor above 2e7 —
pure unpriced slack sitting in the success ledger. Halving frees
4.768e-7 of success probability, which moves the operating point nine
1e-7-grid steps down (tau 4943026 to 4943017) at the re-declared
1.0e-8 floor: model term 0.3900031009624625, p_claim
0.39000005189416408, margin 5.1894e-8 above 0.39. Nine steps of tau
relief buy ~1.3e-6 of float score (ledger float 127.007262672),
clearing the next rounding step where the floor-only ladder (one
step, ~1.5e-7) could not. No other success term moves: eps, every
DELTA_* debit, KM and the 2^-100 term are unchanged.

## Zero-product fallback cross-reference

Predecessor §1.3 requires the unbatched complete formula when a batch
product is zero; the harness aborts. Production replacement: the
affected walk recomputes via `f_single` (F_MAX = 41190, forced maximum
in measurement/forced/summary.json); remaining lanes keep PLAIN/ANCHOR.
Covered by W_ENTRY inside the KM*Eres reserve (proof.md fallback
section); the ledger asserts the covering inequality. No success
debit: the walk continues with the correct successor.
