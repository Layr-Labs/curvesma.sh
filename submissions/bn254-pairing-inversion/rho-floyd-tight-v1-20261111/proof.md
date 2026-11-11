# Tightened bounded Floyd Pollard-rho for BN254 GT DLP

This is a generic, small-memory 32-adding Pollard-rho attack on the discrete
logarithm in the full protected 254-bit prime-order BN254 target group GT. It is
**not** a best-known finite-field attack, a novelty claim, a practical break, or a
security-bit estimate. The work bound is an **every-run cap**; the success bound is
analytic and uses the standard, explicitly disclosed r-adding random-walk
heuristic. This submission seeks the protected exploratory `plausible_not_refuted`
standard, not formal proof certification. No prior result is assumed accepted.

The algorithm and success argument are the same as the organizer's selected
reference `rho-baseline-v2-20261110` (the accepted 135.216266 control). The only
change is that the per-iteration **instruction charge** is a concrete, itemized,
conservative bound grounded in the organizer's own pinned measurements, instead of
the baseline's deliberately loose `3,420,160`. Every heuristic, cap and failure
rule is preserved. The resulting score is **130.868186**.

## Claim

For the full protected BN254 target and uniform `z = g^x`, the concrete one-trial
32-adding Floyd Pollard-rho algorithm below returns `x` with `g^x = z` at
exploratory success probability **at least 0.50**, using **at most
`2^130.868186`** reference units (`C = 20,000` scalar AArch64 instruction
equivalents per unit) and **at most `2^20` bytes** of peak total storage. The
uniform capped instruction bound is

```
R_total = N * PER_ITER + SETUP,   W = R_total / 20,000,   time_log2 = ceil6(log2 W),
N        = 2 * ceil(sqrt(r)),                 # hard Floyd iteration cap
PER_ITER = 3 * 53,248 + 8,192 = 167,936,      # per Floyd iteration
SETUP    = 2^40.                              # one-time setup/recovery/output
```

The work and memory caps hold on **every** run (there is no unbounded randomized
subroutine and no run gate). The success lower bound is exploratory, relying on the
disclosed collision-time heuristic, supported and limited below. Total work
includes all one-target preprocessing, fresh randomness, failed outcomes, hashing,
arithmetic, memory access, recovery, verification and output. No prior acceptance,
searched advice, amortization or production-scale execution is assumed.

## Complete algorithm

Read only the public profile and `z`. Decode its 12 canonical Fp12 coefficients in
the pinned basis `Fp[w]/(w^12 - 18 w^6 + 82)`, check `z != 0` and `z^r = 1`, and
convert to four-limb Montgomery coefficients. If `z = 1` return scalar `0` (or, on
the pairing track, `P0` and `Q=infinity`). Let `N = 2 * ceil(sqrt(r))` for the
protected `r`. There is one bounded trial.

Sample a fresh 256-bit salt and `66` independent uniform scalars modulo `r` by
bounded 254-bit rejection sampling (at most `1024` draws per scalar; one draw
accepts with probability `r / 2^254 > 1/2`). For each of 32 pairs `(c_j, d_j)`
build `M_j = g^{c_j} * z^{d_j}`. From the remaining pair form `Y0 = g^{a0} z^{b0}`.
A state is `(Y, a, b)` with the invariant `Y = g^a z^b`. The deterministic
iteration is

```
j = SHA256(salt || fixed_encoding(Y))[0] & 31
(Y, a, b) <- (Y * M_j, (a + c_j) mod r, (b + d_j) mod r).
```

`fixed_encoding` is the injective 384-byte internal Montgomery encoding specified
in `machine.md`; it is a priced routine, not a free hash oracle. Initialize
tortoise and hare to the same `(Y0, a0, b0)`. For `k = 1, ..., N` advance the
tortoise once and the hare twice and compare all 12 GT coefficients. At the first
equality let `(a, b)` and `(A, B)` be the two label pairs. If `(B - b) = 0 mod r`,
report failure. Otherwise compute `x = (a - A) * (B - b)^{-1} mod r`. Verify
`g^x = z` before output; a verification failure is a failed trial, never an
accepted wrong answer. If there is no equality within `N` iterations, report
failure. There is no hidden restart, uncharged amplification or amortization.
`rho.py` is the executable reference for this exact loop order on cyclic groups.

**Pairing-inversion track.** The identical walk recovers `x`; the output is
`P = P0` and `Q = [x] Q0`, computed by the explicitly charged affine double-and-add
in `machine.md`. Bilinearity gives `e(P0, [x] Q0) = g^x = z`, and both points lie
in the prescribed order-`r` subgroups (the subgroup checks `[r]P = [r]Q =
infinity` are charged under `SETUP`). For the identity, `Q = infinity` is returned.
This is permitted by the arbitrary-preimage-pair rules.

## Correctness given a nondegenerate collision

The invariant `Y = g^a z^b` holds initially and after each step by exponent
addition. If `g^a z^b = g^A z^B` and `z = g^x`, then `a - A = x (B - b) mod r`.
Since `r` is prime, a nonzero `B - b` has a unique inverse and yields the unique
discrete logarithm `x`. All scalar operations use four 64-bit words reduced modulo
`r`. No sampled scalar and no source-group preimage is an input; the final
exponentiation validates every returned scalar.

For a deterministic state map, write `mu` for the tail length and `lambda` for the
cycle length. Floyd detects an equality at some iteration `k` that is a positive
multiple of `lambda` with `k >= mu`, and such a `k <= mu + lambda` exists.
A first-repeat length `mu + lambda <= N` therefore guarantees detection within the
charged `N` iterations; scalar labels are not part of the hash key, so detection is
on group elements only.

## Success heuristic, support and limits

For an ideal uniform random map on `r` points with uniform start, the no-repeat
probability after `t` distinct transitions is `prod (1 - j/r) ~ exp(-t^2/(2r))` at
`t = Theta(sqrt(r))`. A 32-adding map is not an independent random map; the
conventional r-adding correction accounts for re-selecting the same increment and
gives an effective rate `~ 1 - 1/32`. Applying it at `t = 2 ceil(sqrt(r))` predicts
a repeat probability `1 - exp(-(31/32) * 2) ~ 0.856`. This is a stated heuristic
prediction, not an exact tail theorem for this fixed field and SHA-256; it leaves
substantial margin over the declared **success >= 0.50** and the protected minimum
`0.39`.

**Degeneracy is exact, not heuristic.** Fix the true `x` (not given to the
algorithm) and condition on the salt, the 32 increment logarithms
`s_j = c_j + x d_j mod r`, and the start logarithm `s0 = a0 + x b0 mod r`.
Conditional on these, each `d_j` remains independent uniform: for each `d_j` there
is exactly one `c_j` realizing the conditioned `s_j`. The group walk, its partition
choices and the selected first Floyd equality depend only on `s_j, s0, salt`, not on
the remaining `d_j`. At a detected collision, `B - b = sum n_j d_j mod r` with
nonnegative counts summing to `k`, `0 < k <= N < r`, so at least one `n_j` is
nonzero mod `r`; a nonzero linear combination of independent uniform field elements
is uniform, so the degeneracy probability is exactly `1/r`, conditional on any such
walk and on the actual deterministic SHA partition. Bounded rejection sampling is
uniform conditional on success; its failure probability is below `66 * 2^-1024`.

**The only heuristic premise** is that the salted SHA partition and the randomly
generated increments do not cause anomalous collision-time behavior or slow mixing
at production size. SHA collision resistance alone does not prove this. It concerns
walk statistics across fresh increments, salt, labels and uniform production input,
not a claim that any fixed seed is representative or that the whole walk is i.i.d.
The work and memory caps hold even if this heuristic fails; only the success lower
bound is exploratory.

**Support.** The standard r-adding construction and its random-map analysis
(Teske, *On random walks for Pollard's rho method*, Math. Comp. 70(234), 809-825,
2001, DOI 10.1090/S0025-5718-00-01213-8; 32 additions, not Pollard's less-random
three-part walk), the sizeable probability margin, and the supplied experiment,
which runs the exact loop, partition, 384-byte encoding, cap, recovery and failure
rules on four small cyclic groups and reports every outcome including failures
(observed success fractions `0.90`-`0.95`, well above `0.50`, with mean iteration
counts tracking `sqrt(q)`). These are finite walk-statistics and correctness
evidence only; additive-group encodings differ from GT, and extrapolation to GT and
`2^127`-scale walks is part of the disclosed heuristic. Teske is cited as context,
not an imported external theorem; the reasoning is reproduced here.

## Computation, memory and output

`machine.md` gives a complete bounded scalar implementation and the itemized
per-iteration charge. Each of the three Floyd group steps costs at most

```
53,248 = 32,768 (Fp12 multiply) + 16,384 (SHA-256 partition) + 4,096 (updates),
```

and each Floyd iteration adds `8,192` for the 12-coefficient equality test and
multiword counter maintenance, giving `PER_ITER = 3 * 53,248 + 8,192 = 167,936`.
The `32,768` Fp12 charge is a strict upper bound: the organizer's pinned cost model
measures the larger-field blst dense Fp12 multiplication at `22,996` scalar
instructions, and a BN254 (254-bit) Fp12 product is smaller. Every component is a
conservative round-up of the analytic count in `machine.md`; none assigns a BN254
operation a unit cost. Thus

```
R_total = 2 * ceil(sqrt(r)) * 167,936 + 2^40,
W       = R_total / 20,000,   time_log2 = 130.868186 (rounded up, 6 dp).
```

`experiment.py` recomputes `N`, the itemized `PER_ITER`, `R_total` and the
rounded-up score from the protected order.

**Executable at production scale.** The iteration counter is a three-word
(192-bit) integer and the exponent labels are four-word integers, so the bounded
implementation can run the full `N ~ 2^127.8` iterations; no single 64-bit loop
register is relied upon. The attack is a single processor with `O(1)` state; there
is no distributed node machine, message-passing cost model, or hash-table
occupancy argument. Peak storage is a few GT elements, scalars, the 32-entry
multiplier table and the salt: at most `2^20` bytes. Preprocessing `<= 2^36`
scalar instructions is charged in full; nonuniform advice is zero. Only public
profile and SHA constants and algorithm design are excluded.

## What is tightened relative to the accepted baseline

The baseline `rho-baseline-v2-20261110` charges `3,420,160` instructions per Floyd
iteration (its per-step Fp12 multiply alone is charged `2^20 = 1,048,576`, about
`45x` the organizer's own measured `22,996`). This submission keeps the baseline's
algorithm, salted partition, cap, degeneracy proof, success heuristic and memory
bound unchanged, and replaces only that loose charge with the itemized conservative
bound above (`167,936` per iteration). The cost model explicitly permits analytic
instruction-equivalent bounds. No algorithmic, success or memory claim is
strengthened beyond the accepted baseline.

## Provenance

The method is the standard generic 32-adding Pollard rho on a known-order group,
following the organizer's selected reference baseline `rho-baseline-v2-20261110`
(repository `submissions/bn254-gt-dlp/`) for the algorithm, partition, degeneracy
argument and experiment design. The per-iteration instruction bound is derived
independently in `machine.md` and anchored to the organizer's pinned blst
measurement in `tracks/docs/COST_MODEL.md` (`22,996` instructions per dense Fp12
multiply). No other participant's submission is used. Rho with the same asymptotic
`~ sqrt(r)` step count is standard; this submission does not claim the Frobenius
`sqrt(12)` class speedup or the negation-map `sqrt(2)` speedup, keeping the
heuristic set minimal.
