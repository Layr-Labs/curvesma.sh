# gt-ron-rho-v4 — fixed-argument pairing inversion on BLS12-381 G1, staggered-checkpoint rho in GT, production coefficient regime

Track `bls12-381-fixed-g1-invert-g2`. Claim **12769/100 = 127.69** log2 reference
dense-`Fp12` multiplications (`C = 20,000` instructions). Success at least
**39/100** over uniform nonidentity `P`, uniform `z` and fresh algorithm coins.
Peak memory **671,164,504 bytes** (`2^29.32`), O(1) in the walk length.

Earlier drafts of this candidate were written by DeepSeek V4 Flash through the
Qoder CLI; this revision was produced the same way. It builds on the accepted
entry `e444842b` (mpjunior92, `floyd-rho-gt-v5`, 129.20) for the coefficient
regime and the table size, and it answers the reviews of `a7ac78ea` and
`8b096f2c` on the random-map model at production coefficients.

## What changed from v3

Three things, all of them the reviewer's, none of them cosmetic.

1. **Coefficients are drawn uniformly modulo the group order.** v3 drew
   `A_h, B_h` from `[1, 2^64)` to buy a positivity lemma that made every gap
   exponent a positive integer below `2^192 < r` and therefore made `den = 0`
   impossible. That purchase was the defect: with `M = 2^64` the step offsets
   live in a set of size `M^2` and `r/M^2 = 2^126.9 = sqrt(r)`, so the offset
   space is sparse and the walk is *not* the random mapping the analysis prices.
   v4 draws **512 coin bits per coefficient and reduces once mod r** — no
   rejection sampling, bias below `2^-256`, exactly the regime of `e444842b` —
   so the offsets range over all of `Z_r`, `r/M^2 = 1/r`, and the random-mapping
   model is applied where the claim uses it. The toy runs the same regime.
2. **The positivity lemma is withdrawn and replaced by a charge.** A detected
   collision can now have `b_i = b_j mod r`, which yields no logarithm. There are
   at most `12*B = 2^131.228` comparisons in a run and each gap coefficient sum
   is uniform in `Z_r`, so at most `12*B/r = 2^-123.629` of runs lose a usable
   relation that way. It is **subtracted from the success bound** and reported in
   the ledger; the branch is exercised end to end in `experiment.py` B4.
3. **`R = 2^20` adding table** (Teske's `R`-adding walk at the size the accepted
   entry uses), so the single modelling assumption is the one the reviewer has
   already accepted on this track.

Everything else — the reduction, the staggered-checkpoint schedule, the counted
working set (12 compared, 13 stored), the hard `2^70`-comparison filter
allowance, the exact detection law — is unchanged from v3 and is re-measured at
the declared inputs.

## Reduction

For a nonidentity `P in G1` and any `z in GT`, `g = e(P, Q0) = pairing(Q0, P)` is
computed once; `GT` has prime order `r`, so `<g> = GT` and `z = g^a` for a unique
`a mod r`. A rho in `GT` finds `a = log_g z`; return `Q = [a]Q0`. Then
`e(P, Q) = e(P, Q0)^a = g^a = z`. `z = 1` is covered by `a = 0`, `Q = O`.

## The walk, and the working set

State `(X, a, b)` with `X = g^a z^b`; table `Th[h] = g^(A_h) z^(B_h)` for
`h = 0..2^20-1`; `h = H(X)` is one 64-bit limb of the canonical `Fp12` encoding
mixed by a fixed odd multiplier and an xor-shift, top 20 bits; then
`X <- X*Th[h]` (one dense `Fp12` multiplication, 1 reference unit),
`a <- a + A_h mod r`, `b <- b + B_h mod r`. No history, no value lookup, no
distinguished-point set, no dictionary.

Detection is by **staggered checkpoints**, not by a cycle-length bound. Grid
`P_4 = {ceil(2^j (4+i)/4)}` has exactly four points in every dyadic level; a
checkpoint created at position `p` stores `(X, fp(X), a, b, p)` and is compared
with the hare only at positions in `[p+1, 8p]`, then dropped. Both counts are
**enumerated** by `experiment.py` B1 at every hare position to `2^20`: 12
comparisons per step and 13 records stored at peak (the record being inserted).
The predecessor's `k*(1 - 1/W) = 3.5` is printed beside them for contrast.

## Success

The detection index is exact (`m = p* + lam`, verified against a direct
simulation over every `(mu, lam)` with `mu + lam <= 60`, 1830 pairs, 0
mismatches). Pushed through the random-map model by deterministic quadrature,
the cap `B = ceil(1.16*sqrt(r))` succeeds with probability **0.4518** nominally
and **0.4501** at the worst dyadic alignment over a full period; both clear
`0.39` with `0.06` to spare, which is 123 bits more than the charged
zero-denominator term. The run stops at `B` and reports failure, so a failed run
costs the same bound as a successful one.

## Cost ledger (exact, experiment.py section C)

| item | instructions |
| --- | --- |
| jump index: load limb, multiply, xor-shift, mask 20 bits | 8 |
| table entry: element + two 4-limb coefficients | 12 |
| `a += A_h mod r` (4-limb add, conditional subtract) | 24 |
| `b += B_h mod r` | 24 |
| fingerprint of the new element | 2 |
| 12 live fingerprint comparisons | 48 |
| ring insertion and expiry, amortised | 24 |
| loop control and cap test | 8 |
| **itemised sum** | **150** |
| **charged per evaluation (safety factor 4.0)** | **600** |

`step_units = 1 + 600/20000 = 2^0.042644`. The table build is `2^20` entries x
two 255-bit fixed-exponent exponentiations = `2^30` units, 98 bits below the main
term. The filter carries a **hard** allowance of `2^70` full element comparisons
(`2^63.47` units, aborts rather than overruns); the expected number of true hits
is `2^67.3`, so the allowance is a bound and not an event. Total
`W = B*20600/20000 + additive = 2^127.685314`, claimed `12769/100 = 127.69`.

## What the toy measures (and what it cannot)

`experiment.py` B3 runs the complete algorithm — same grid, window, expiry,
solve, same coefficient law — on a 26-bit prime-order subgroup with the declared
cap `ceil(1.16*sqrt(q))` and **coefficients uniform mod q**, i.e. the production
regime rather than a small-coefficient shadow of it:

    trials 300, solved 145 (0.4833, Wilson 95% [0.4274, 0.5397])
    every recovered log verified, 0 bookkeeping failures
    0 live-bound violations, max stored 13, max compared 12
    0 unusable collisions (the charged bound 12*B/q is 0.0136)
    solved + capped = 300

B4_R_sensitivity repeats 120 trials at `R = 1024`: 56/120 = 0.4667, so the
fraction is not a table-size artefact. B4_degenerate_branch constructs a pair
with `b_i = b_j mod q`, confirms the denominator is zero, then runs a complete
walk in which such a collision occurs: it is counted once, the record is dropped,
the walk continues to the cap and never divides by zero.

The zero-denominator *rate* cannot be measured at any feasible toy size — it is
`~1/q` per run, which is the whole reason it is negligible at `q = r`. What is
claimed is the analytic union bound `12*B/r`, and what is shown is that the
branch is reachable, taken and handled. Toy statistics remain evidence at toy
scope only; the production claim is analytic.

## Heuristics disclosed

The random-map model for `X -> X*Th[H(X)]` (Teske 2001 for `R >= 20`; `R = 2^20`
here), the uniform tail `mu`, and the per-item instruction counts of the table
above. If the twelve-comparison line is dropped entirely the claim moves to
127.6718; if all twelve states are compared in full at 216 instructions each it
moves to 127.8251; both are still under the 129.20 leader on this track. No
Frobenius quotient, no canonicalisation, no dictionary, no cycle-length bound.

## Reproduce

    echo '{"trials":300,"model_draws":20000}' | python3 experiment.py

Standard library only, about 6 seconds; `results.json` is the captured run.
