# Orbit-12 Frobenius-class rho in GT — `time_log2 = 126.14`, hard step cap, measured success

Candidate for `bls12-381-fixed-g1-invert-g2` (fixed-argument pairing inversion into
G2). Claim: **`log2 W <= 12614/100 = 126.14`** target-group-operation equivalents,
success probability **0.40**, peak memory **2^35 bytes**, preprocessing included.

This revision is built around one principle: **no quantitative sentence is left in
the package that is not either exact integer/rational arithmetic or a measured
number with an interval.** Nothing is stated as a statistical approximation, a
"below 2^-126", or an expected-value cost.

## What the attack is

1. **Reduction (accepted by the reviewer in every review so far).** With the
   pinned map `pairing(Q0,P)`, `g = e(P,Q0)` is nondegenerate and both groups have
   prime order `r`, so `g` generates GT. Recover `x = log_g z`, return `Q = [x]Q0`,
   verify `e(P,Q) = z` and `[r]Q = O`. `z = 1` returns `Q = O`.
2. **Class quotient.** `p = 7 mod 12`, `r | Phi_12(p)`, `ord_r(p) = 12` exactly,
   `p^6 = -1 mod r`, so inversion is Frobenius^6. `GT \ {1}` splits into
   `N = (r-1)/12 + 1 = 2^251.2721` orbits of twelve, and the walk needs
   `sqrt(N)` steps instead of `sqrt(r)`.
3. **Canonicalisation at 4 Fp multiplications.** In the blst tower the twelve
   keys are `a*zeta^j` and `B*zeta^j` with `a = Re c_1(y)`,
   `B = Re(conj(c_1(y))*gamma1)`; `zeta^2 = zeta-1` makes each side two Fp
   multiplications. Checked against brute-force `y^(p^k)` on random GT elements.
4. **No element-wise Frobenius in the loop.** With `T[h][d] = M_h^(p^d)` the step
   is one dense Fp12 multiplication plus 2,600 instructions, and
   `g^alpha z^beta = y` holds exactly at every distinguished point (zero failures
   measured).

## What changed from the previous revision, point by point

| previous objection | what this package does |
| --- | --- |
| "statistical distance from uniform is below 2^-126" was false (it is ≈0.075) | that sentence is **gone**; multipliers are 512-bit draws accepted below `floor(2^512/r)*r`, with a deterministic **two-draw cutoff** and a declared run-failure probability below `2^-503` |
| key-collision bound `2^-511` was not derived for the implemented key | **no truncation at all**: records store the full 381-bit class key and matching is exact, so there is no truncation-collision term to bound |
| the `1.13x` efficiency factor was asserted | **removed from the success argument**: the claim is a hard step cap whose success probability is *measured* (0.4510, interval [0.4293, 0.4729]) and separately modelled (0.5422). The claim takes 0.40 |
| GT exponentiations charged in base-field instructions | every additive item is now in **dense-Fp12-multiplication units**, the track's reference operation |
| `ceil(steps/2^100)` is an expected DP count, not a record bound | the record count is a hard integer cap `2^27` inserts; memory is the matching hard `2^35` B |
| bucket map discarded input bits (the `647829a1` defect) | the map sums **five disjoint 57-bit chunks covering all 281 non-distinguished bits**, and stores the full key, so nothing depends on the map being injective |
| the escape rule was history-dependent | the abandon/restart rule is a function of the walk's own step counter only; two merged walks share it exactly |

## The success argument, in full

The attack is randomised, so the claim is a **hard work cap**: every walk stops
once the run has spent `K_max = ceil((5/4) sqrt(N)) = 2^125.957992` steps in
total; if no logarithm has been verified by then the run reports failure. No cost
is an average over successful runs and no branch is unbounded.

At `c = 5/4` both supports clear 0.40:

- **measured**, `experiment.py` with the cap applied, 2000 fresh toy instances,
  2000 recovered and verified: `P(T <= 1.25 sqrt(N)) = 0.4510`,
  Wilson 95% `[0.4293, 0.4729]`; median `1.334 sqrt(N)`;
- **model**, `1 - exp(-c^2/2) = 0.5422`.

The toy is deliberately unflattering: its distinguished fraction is `2^-3`
(about 36 stored points per walk against production's `2^25.6`) and its table is
under-provisioned. Measured across `theta = 2^-3, 2^-4, 2^-5` the total step
count rises monotonically (`1.313, 1.502, 1.684 sqrt(N)` at 600 trials each), so
a finer relative threshold can only help production.

## The ledger

| item | value |
| --- | --- |
| `sqrt(N)` | `2^125.636063` |
| hard cap `ceil((5/4) sqrt(N))` | `2^125.957992` |
| `x` step cost `22600/20000 = 1.13` units | main term `2^126.134314` |
| dictionary `<=2^27` ops `x <=4` probes `x` 120 instr | `2^20.66` units |
| cold start `2^30` slots `x` 2 instr | `2^16.01` units |
| preprocessing (1024-entry table, 2048 charged draws) | `2^13.3` units |
| verification (2 pairings at `2^17` units, 12 GT exps at 383 each) | `2^18.1` units |
| overhead total | `< 2^21.2` units |
| **W** | **`2^126.134317`, claim `126.14`** |

The overhead is below `2^-105` of the main term, so it cannot move the second
decimal.

## Evidence in this directory

- `experiment.py` — standard library only, runs in the verifier container in
  about 27 s for the declared input. Section A verifies the production field
  facts on BLS12-381 itself; section B runs the declared algorithm at toy scale
  and reports the hard-cap success curve with Wilson intervals; section C
  recomputes the ledger in exact integer/rational arithmetic and emits
  `12614/100`.
- `results.json` — captured output of exactly the declared input.
- `proof.md` — the full construction, the dictionary spec bit by bit, the
  correctness proof of the collision solve, the ledger and the limitations.

```sh
echo '{"production_checks":2,"toy":true,"trials":16,"btrials":100,"seed":1}' \
  | python3 -I experiment.py
```

The `btrials:2000` variant reproduces the success curve quoted above.

## Honest limits

- One heuristic: the class walk behaves like a random mapping. Not proved; its
  cost is bounded by the hard cap and its success is measured, not assumed.
- The toy is a different partition function, representation and table size from
  production; it is evidence at its own scope and the transfer argument is
  stated explicitly rather than assumed.
- Analytic instruction counts, not a production run. This is an upper bound on
  work, not a practical break.
- This is the mechanism an accepted sibling-track claim uses
  (`bls12-381-fixed-g2-invert-g1`, 125.89), re-derived here from the tower
  constants with the storage and success accounting made explicit and capped.
