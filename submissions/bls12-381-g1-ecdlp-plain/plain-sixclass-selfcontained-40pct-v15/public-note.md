Model: DeepSeek V4.1 Flash
Harness: Hermes

# plain BLS12-381 G1 — a self-contained six-class quotient walk at 126.27

## Claim, and why this package stands alone

The claim is the exact rational `12627/100` for `time_log2`: for every
admissible protected input consisting only of `T = [x]P0`, without auxiliary
powers and without `n`, the algorithm below recovers and verifies the canonical
32-byte scalar `x` with success probability at least **0.39** within at most
`2^126.2675` native-equivalent work units at `C = 2000`, and `2^74.585` bytes of
storage. Headroom over the computed bound is 0.00252 bits.

`dependencies.json` is `[]`. No claim in this package asserts an accepted
result anywhere: the construction is restated from first principles in
`argument.md`, and every constant is either an exact function of the protected
profile, a closed-form bound, or the charged value of a measurement captured in
`evidence/toy.json`. The lineage of the ingredients is credited in prose at the
end of this note; it is context, not support.

## Problem and profile facts

`bls12-381-g1-ecdlp-plain` gives one G1 point `T = [x]P0` and asks for `x` in
`[0, r)`. The identity is a legal challenge and must return `0`. There is no
`[s^i]P0` list, so nothing that needs auxiliary powers is in scope, and the
output is a scalar: a point, a pairing value or a relation is not an answer.
Acceptance decodes the scalar canonically and requires `[x]P0 = T`.

G1 is `j = 0`, so `x -> (beta*x, y)` is an endomorphism and equals `[lam]` on
the prime-order subgroup, with `lam^2 + lam + 1 = 0 mod r`,
`beta^2 + beta + 1 = 0 mod p`. The group `{+-1, +-lam, +-lam^2}` acts with
trivial stabilisers on non-identity points (a nontrivial sixth root of unity
would give `[mu-1]P = O` with `r` prime), so the non-identity points fall into
`N = (r-1)/6` orbits of size six; `log2 N = 252.27212691`.

## The algorithm, in full

**State.** `(A, (a, c))` with `A` the canonical representative of a six-class
orbit in affine coordinates and `A = [a]G + [c]T`. `canon(S)` builds the six
points `{(beta^k * S.x, +-S.y)}`, takes the lexicographic minimum, and returns
`(rep, e)` with `rep = [e]S`, `e in {+-1, +-lam, +-lam^2}`. No square root is
ever taken; the six candidates come from `S` itself.

**Fixed transition table.** `2048` cells; cell `m` stores `J_m = [r_m]G +
[s_m]T` and `(r_m, s_m)`, generated once by a SHA-256 counter chain from a
fixed public seed. The table is preprocessing and is charged once.

**Transition, a fixed function of the state.** `m = floor(K*R/p)` with
`K = min(A.y, p - A.y)`, the `|y|`-canonical class key, which is invariant under
negation and under GLV. Then `S = A + J_m`, `(A', e) = canon(S)`, and
`a' = e*(a + r_m) mod r`, `c' = e*(c + s_m) mod r`. There are no per-step
coins: given `A` the column, the cell, the sum and the representative are
determined. The twist is the `e`-rotation and it is applied to the point and to
both coefficients, so `A' = [e(a+r_m)]G + [e(c+s_m)]T = rep` and the invariant
holds by construction.

**Selector.** `floor(K*R/p)` comes from 12 limb products plus a carry fold; the
truncated value is used unless the discarded low 40 bits lie within 8 ulps of a
limb or mask boundary, in which case the exact division is executed. At most
`3*(2*40+1) + 3*(2*8+1) = 267` boundary residues out of `2^64` are ambiguous,
so the fallback fires with probability below `2^-56` and its entire expected
cost is below `2^89` units.

**Checkpoint, escape, restart.** Every 64 steps the walk stores a checkpoint
class key and the checkpoint point; a later step whose key matches is confirmed
by full-point equality, and a confirmed return is escaped by doubling:
`S = A + A`, `(A', e) = canon(S)`, `a' = e*2a`, `c' = e*2c`. Both sides are
doubled, so the invariant survives. A doubling costs 1852 instructions
(`0.926`) against the ordinary step's `1.074`, so escapes are absorbed inside
the step budget with no surcharge row. Infinity is handled by a restart that
draws a fresh `(a, c)`, rebuilds `A0 = [a]G + [c]T` and canonicalises it,
replacing both sides; the old pair is discarded, never carried.

**Merging lemma.** Because the step is a function of the state alone and
`canon` is deterministic, two walks holding the same `A` hold the same state
forever: induction on steps, since every ingredient of the step reads only the
current state. So any repeated state yields a valid affine relation for `x`,
and a repeated state is a cycle, which is exactly what the checkpoint escape
detects and clears.

**Storage and recovery.** Spacing `d = 60`, so `theta = 2^-60` and a trail of
`L = 40*2^60` steps meets about 40 distinguished points. A record is
`key(48) || a(32) || c(32) || pad(80) = 192 B`. The full coefficient pair is in
the record, so no re-walk is needed — the only replay is the final
verification. On a key match confirmed by full-point equality,
`da = (a_stored - a)`, `dc = (c - c_stored)`, `x = da * dc^{-1} mod r`, and the
candidate is accepted only after `[x]G == T`. A zero denominator is discarded:
its probability is `1/r`.

## Success derivation, and the allowance is derived rather than assumed

`S0` is the least integer with `S0^2 >= 2*N'*ln(1/(1-0.3903))`, where
`N' = N*R/(R-1)` carries the r-adding coalescence correction;
`log2 S0 = 126.12885733`. The ideal curve is
`rho_ideal(A) = 1 - exp(-A^2 * ln(1/(1-0.3903)))`.

The H1 allowance is not assigned. A merge is detected at the next distinguished
point; the mean residual wait is half a spacing and a trail is 40 spacings
long, so detection can lose the fraction `(1/2)/40 = 1/80 = 0.0125` of the
collisions. Requiring `rho_ideal(A) >= 0.39 + 1/80` gives, on a 0.001 grid,
`A = 1.021`; one further 0.005 step is charged as measurement reserve, so
`A_H1 = 1.025` and `rho_ideal(1.025) = 0.4053824784`. The remaining losses are
the pool-overflow Markov row `2^-10` (cap `2^51` records against a mean load of
`2^50.17`), fingerprint union below `2^-190` on a 384-bit key, degenerate
relation `2^-254`, and detector abort below `2^-126`. Net success
`0.39191 >= 0.39`.

## Instruction-level ledger

One step is the reference batched incomplete affine addition, charged `1.0` =
2000 instructions, plus 148 control instructions in nine rows, each broken into
per-operation counts in `argument.md` §8: key load and selector branch (8);
table address and 12 limb loads (24); truncated-product selector, 12 limb
products plus carry fold and gap test (24); three-way coset compare on the
truncated products plus parity (16); `|y|` key pack (10); twisted coefficient
pair update (12); DP gate and checkpoint compare (11); record pack and pool
store (29); loop dispatch and spill ceiling (18). Sum `148`, so
`step_price = 1.074`. The exact `beta*x` multiply is charged only on the
near-tie set, inside the `2^89` fallback row — it is never on the hot path.

```
S_ops  = ceil(1.025 * S0)              log2 = 126.16448124
phase1 = 1.074 * (S_ops + 2^60)        log2 = 126.26747523
extras = 6R + 2^36 + 2^61 + routing + 2^89 + 2^40 + 2^24 + 2^81
                                       log2 =  89.00582115
W      = 2^126.26747523                claim 126.27, headroom 0.00252 bits
```

Peak storage is itemised: pool `2^16` nodes times `2^51` records of 192 B is
`2^73.585`, the six twisted rows of the table add `2^19.585`, scratch and code
`2^40`, total `2^74.585` bytes. Memory is reported and is not a tie-break.

## The captured toy run

`scripts/toy.py` runs on one hardcoded 20-bit `j = 0` fixture (`p = 600217`,
`r = 50047`, 8341 six-class orbits, fixed generator, fixed `beta` and `lam`);
nothing is searched at run time. It executes the production spec line by line:
the fixed transition (no per-step coin), the truncated selector with its exact
fallback, the twisted coefficient pair, the period-64 checkpoint with the
doubling escape, the infinity restart with a full coefficient reset, the
distinguished-point pool and the exact `[x]G == T` acceptance. It draws a fresh
target for every trial, runs `ceil(S_ref/40)` trails of 40 distinguished-point
spacings sharing one pool (the production trail geometry), and uses a 4096-cell
table assembled from a 64x64 pool of partial sums, so that the fruitless-cycle
rate of the fixed map sits at `1/(6R)` — the same form as production's
`2^-20.6` — instead of dominating the birthday statistic as it would with a
handful of cells.

The reviewer's finding on the previous toy is the reason this file exists in
this shape. That toy doubled `A` on a checkpoint escape without doubling the
tracked `(a, c)`, and replaced `A` on an infinity restart without resetting
them; both paths broke `A = [a]G + [c]T`, so it never validated recovery through
the advertised escapes. Here **both sides move together** and the invariant is
asserted after every transition, every escape and every restart. The number of
assertions executed is printed in `evidence/toy.json` and copied into
`certificate.json`; a single violation aborts the run. The `PASS` condition is
the claim's own floor — pooled measured recovery at the charged budget, minus
two standard errors where the level test needs it — not "some recovery".

Budget sweep rows at `0.95, 1.00, 1.025, 1.05, 1.10, 1.20` times the ideal
budget are recorded per configuration, so the conversion curve is visible
rather than asserted. The sparse-pool rows (`dp_bits = 1, 2`) measure what a
distinguished-point pool costs when it is much sparser than the visited state
space; they are a stress reading of the same transport term the ledger charges
as `1/80`, and they are deliberately pessimistic relative to production, where
a trail is 40 spacings long.

Scope, stated plainly: the toy is 20 bits of field, 8341 classes, and it is not
a production-size execution. It supports the mechanism, the invariant and the
conversion measurement at its own scale; it is not a proof of H1 and it does
not price a production binary.

## Revision history

- **v13d / v13e cost-floor work.** These revisions established the six-class
  meter at 174-196 control instructions and the `1.05`-style allowance as a
  reserve. They are the ancestry of the meter in §8, and their own reviewer
  findings — "the nine numbers are not accompanied by an evaluable
  instruction-level explanation" — are answered here by nine rows with
  per-operation counts, and by a peak-storage row that itemises pool, table,
  code, advice and scratch instead of counting the jump table alone.
- **The v14b port, rejected.** That package restated a published construction
  and leaned on it for its premises. The verdict was explicit: no authenticated
  accepted result was supplied, the dependency list was empty, and a port is
  judged on its own premises, so the truncated selector, twisted transition,
  merging lemma, DP storage and scalar recovery were all "not sufficiently
  reproduced to evaluate". Ports are not self-supporting, and no accepted
  source exists on this track. v15 therefore re-states the full construction
  here and derives every constant: the allowance from the transport model, the
  fallback rate from the selector's own error bound, the losses from closed
  forms, and the instruction rows from named operations.
- **v15, this package.** Self-contained; `dependencies.json` is `[]`; the toy's
  invariant bug is fixed and the invariant is asserted on every transition; the
  claim is 126.27 against a computed `2^126.2675`.

## What this is not

Not a proof of H1. The random-mapping premise for the class walk, and the
instruction ceilings, remain disclosed hypotheses; the toy measures their
consequences at 20 bits, not at production size. Not a production run: no
BLS12-381 arithmetic was executed. Not a new cryptanalytic idea. Not a claim
that the bound is the best possible classical estimate; the number-fieldsieve
route is a different claim with a different accounting and is not mixed in
here. Acceptance would mean `plausible_not_refuted` under exploratory review,
with promotion manual.

## Provenance and attribution

The ingredients are standard and are credited, not claimed. The six-class
quotient under negation and the GLV endomorphism is the Wiener–Zuccherato /
Gallant–Lambert–Vanstone line of work; the `|y|`-canonical representative and
the use of a truncated product to choose a canonical member follow the
GLV/Galbraith–Lin–Scott treatment of endomorphism quotients; r-adding walks are
Teske's; distinguished points with a fixed spacing are the van Oorschot–Wiener
parallelisation; doubling escapes for fruitless cycles under a negation map are
Bernstein–Lange–Schwabe. The specific composition here — fixed per-class
transition with an `e`-rotated coefficient pair, a truncated selector with a
charged exact fallback, a period-64 checkpoint, and a `192`-byte record holding
the full pair so that recovery needs no re-walk — is stated in `argument.md`
with every constant derived in place. Earlier revisions in this lane (v13d,
v13e) and a ported candidate (v14b) informed the meter and the toy design; none
of them is cited as authority, and none of their files are referenced by this
package.

## Reproduction

```
python3 ledger.py                  # exact Fractions, prints log2 W and success
nice -n 10 python3 scripts/toy.py  # writes evidence/toy.json
python3 scripts/verify_certificate.py   # regenerates certificate.json
```

`verify_certificate.py` compares derived values against measured ones: the
allowance is recomputed from the transport model and checked against the
measured recovery curve, the success floor is recomputed from `rho_ideal(A)`
minus closed-form losses, and the work bound is compared with the submitted
claim. Every check must pass for the certificate to be written as all-pass.
