# H1 collision margins as charged (in-tree quantitative record)

This revision does not import the predecessor's §2.6 (its bundle omits
that section). The collision law H1 is charged here as follows.

## Model term

`1 - exp(-tau*(1-eps))`, the ordinary random-mapping finite-horizon
collision lower bound over T counted main-phase steps, with the eps
fraction of the count discounted before it is applied. At this S:
tau_exact = T^2/(2N) with tau_target = 4943025/10000000,
eps = 5.960475e-7, model term = 0.3900035889594956.

## Eps losses (discounted collision exposure)

| term | value at A=2^24 | identification |
|---|---|---|
| 4/A (2-cycle loss) | 2.384186e-7 | fresh-source 2-cycle entries at exact rate 1/(2A) under L-ori; each removes the entry step and its escape partner (2 steps x 2 orientations in §2.3's exposure process) |
| 2*(3/A) (orientation/closure margin) | 3.576279e-7 | at most 3/A of continuations lose one orientation bit each through the L-ori conditional sign and the canon representative choice; discounted, not re-randomised |
| 2*192*long_rate (longer-entry closure margin) | 2.046364e-12 | Lemma E1(c) rate 3/(4A^2)+2/A^3 = 2.664536e-15 per window; each 3..64 entry wastes 192 states (64-step window x 3-window look-ahead overlap); counted entry- and exit-side |
| total eps | 5.960475e-7 | ledger `eps`, asserted in accounting.py |

## Delta debit (additive 2^-20 = 9.536743e-7)

Covers the failures outside eps: (a) finite-horizon truncation of the
birthday sum below T steps; (b) stopped-process padding after the first
G violation; (c) merge-resolution residue below the DELTA_MERGE bound.
Each is bounded separately in §§2.2-2.4/2.8; their sum at this S is
below 9.536743e-7 (DELTA_CUT = 5.8542e-17, DELTA_MERGE = 1.2520e-14,
truncation/padding terms exponentially small in T), so the uniform
2^-20 debit exceeds them jointly. The escaped, restarted, capped trails
enter only through T; meter aborts and caps convert non-G outcomes into
the DELTA_* debits below, reducing claimed probability, not the law.

## Full success bound at this S (tau = 4943025/10000000)

model 0.3900035889594956 - delta 2^-20 - DELTA_ABANDON 7.2116e-8 -
DELTA_CUT 5.8542e-17 - DELTA_MERGE 1.2520e-14 - DELTA_CAP 2.4393e-18 -
DELTA_PROBE 2.4923e-44 - DELTA_HASH 1.1509e-10 - 1/400000 - 2^-100 =
p_claim 0.39000006305398043, margin 6.3054e-8 above 0.39.

## Zero-product fallback cross-reference

Predecessor §1.3 requires the unbatched complete formula when a batch
product is zero; the harness aborts. Production replacement: the
affected walk recomputes via `f_single` (F_MAX = 41190, forced maximum
in measurement/forced/summary.json); remaining lanes keep PLAIN/ANCHOR.
Covered by W_ENTRY inside the KM*Eres reserve (proof.md fallback
section); the ledger asserts the covering inequality. No success
debit: the walk continues with the correct successor.
