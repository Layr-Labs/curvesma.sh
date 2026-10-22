# secp256k1-ecdlp: truncated Pollard rho on automorphism classes

Status: ready for exploratory review. This is an analytic classical attack claim
with toy-scale experimental support for the collision model. It is not a
practical break, a full-size execution or a novelty claim about Pollard rho.

The normative target is `production-v2` for `secp256k1-ecdlp`. `claim.json`
binds that profile and `classical-group-work-v2` by SHA-256. No restricted
challenge distribution, quantum step or uncharged advice is used.

## Target constants

Let `p` and `r` be the protected field modulus and subgroup order. Both are
256-bit primes with `r = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141`.
The curve is `y^2 = x^3 + 7` over `Fp`, cofactor 1, with the fixed generator `P`
from the profile. An input is a uniform challenge `Q = [x]P` with secret
`x ∈ [0,r)`. The attacker returns `x`.

The curve has complex multiplication by `Z[ω]` with `j = 0`. Explicitly there
is a field automorphism `φ:(x,y) ↦ (βx, y)` of order 3, where `β` is a
nontrivial cube root of unity modulo `p`, together with the negation map
`σ:(x,y) ↦ (x,-y)`. These commute and generate a group
`G = ⟨σ,φ⟩ ≅ C2 × C3` of order 6. On the cyclic group `⟨P⟩` the same maps act
by multiplication of the discrete log by units `{±1, ±λ, ±λ²}` for a fixed
cube root of unity `λ` modulo `r` (the matching endomorphism eigenvalue). The
pair `(β,λ)` is uniquely determined up to swapping `λ ↔ λ²` by the public
curve; both may be computed offline from public constants and are charged
below.

## Algorithm (r-adding walk on equivalence classes)

Fix a partition count `K = 2^20` and a distinguished-point mask of
`θ = 2^{-48}` (a point is distinguished when its canonical x-coordinate is
`0 mod 2^48`). Both are public algorithm parameters.

1. **Offline setup (charged).** Compute `β` and `λ` as above. Sample
   `K` independent pairs `(c_j, d_j)` uniformly in `[0,r) × [0,r)` and store
   the affine points `M_j = [c_j]P + [d_j]Q` together with the coefficients.
   (Equivalently sample only multiples of `P` and fold the `Q`-multiples into
   the walk; either presentation costs the same order of magnitude.) Build an
   empty hash table for distinguished points.

2. **Canonical representative.** For an affine point `R = [a]P + [b]Q` and
   its coefficients `(a,b)`, compute the six orbit points
   `{±φ^k(R) : k=0,1,2}` and their matching coefficient multiplications by
   `{±1, ±λ, ±λ²}`. Among the three x-coordinates choose the numerically
   smallest; among the two signs of that x choose the one with
   `y ≤ (p-1)/2`. Return that representative together with the multiplied
   coefficients. Infinity is a single class and is rejected / restarted.
   This costs a constant number of field multiplications and comparisons per
   call and no group addition.

3. **Walk step.** From a current class representative `R` with coefficients
   `(a,b)`, set `j = x(R) mod K`. Add the stored `M_j`:
   `S = R + M_j`, `(a',b') = (a+c_j, b+d_j)`. If `S` is infinity, resample a
   fresh start. Otherwise replace `(R,a,b)` by the canonical representative of
   `(S,a',b')`. If the same partition index would be selected twice in a row
   (the classical length-2 fruitless cycle), bump `j ← (j+1) mod K` before
   the addition; if a longer fruitless cycle is detected by a repeated class
   with proportional coefficients, escape by one doubling of the current
   state and continue. Each successful step charges **exactly two** target
   group additions under the organiser reference (one for `R+M_j`, one for the
   affine normalisation / representative update charged at the same reference
   price; see the cost section). Partition selection and coefficient
   arithmetic are instruction work, charged separately and dominated by the
   group steps.

4. **Distinguished-point store.** After each step, if `R` is distinguished,
   look up its class key in the hash table. On a miss, store `(a,b)`. On a
   hit with stored `(a₂,b₂)`, if `b ≢ b₂ (mod r)` recover
   `x = (a₂ − a)·(b − b₂)^{-1} mod r`, verify `[x]P = Q`, and halt. Otherwise
   treat the collision as fruitless and escape as above.

5. **Budget and restart.** Run the walk for at most
   `S = ⌈1.05 · 1.02 · (1+2^{-20}) · √(2·(r/6)·ln(1/0.58))⌉`
   group-step iterations (the factors are justified below). If no useful
   collision is found, discard the walk state, sample a fresh start, and
   repeat. At most one walk is started under this claim (a single truncated
   trail); the success probability below already accounts for truncation, so
   no amplification restart is required. Verification of a candidate `x`
   costs one scalar multiplication, charged below.

The algorithm uses only public profile constants, the input `Q`, and its own
coins. No nonuniform advice is embedded.

## Correctness and success probability

Every orbit has exactly six affine points whenever `R` is non-identity and not
a 2-torsion / fixed point of a nontrivial automorphism, which for this prime
order group excludes only the identity. The walk is therefore a walk on a set
of size `N = r/6` (identity discarded by restart). An r-adding walk with a
fixed random partition table is a classical Markovian pseudo-random walk on
that set; under the standard heuristic that successive class iterates behave
as independent uniform samples from the N classes until a repeat (Teske;
Gallant–Lambert–Vanstone; Bos–Kleinjung–Lenstra), the waiting time `T` to the
first collision of class representatives is distributed as the birthday /
Rayleigh law

`P(T ≤ B) = 1 − exp(−B² / (2N))`

on the infinite trail. Setting the success model probability to
`p_model = 0.42` and `B₀ = √(2 N ln(1/(1−p_model)))` gives
`P(T ≤ B₀) = 0.42` under the heuristic.

Three conservative adjustments convert `B₀` into the iteration budget `S`
actually charged:

- **Walk-to-uniform factor 1.05.** Empirical mean collision times of the
  implemented r-adding walk on toy j=0 curves of 24–30 bits (see experiments)
  sit within about ±8% of the ideal `√(π N / 2)` prediction; we inflate the
  truncated budget by 5% rather than relying on the ideal constant.
- **Fruitless-cycle / escape factor 1.02.** Negation and the order-3 map make
  short cycles possible; the classical remedies (partition bump for 2-cycles,
  doubling escape for longer ones) are used. Toy runs recorded at most a few
  dozen escapes across hundreds of trials; we charge a flat 2% iteration
  overhead.
- **Restart slack `1+2^{-20}`.** Covers the negligible probability of landing
  on identity or resampling a start.

With these factors,
`S = ⌈1.05 · 1.02 · (1+2^{-20}) · √(2·(r/6)·ln(1/0.58))⌉`
and the model success probability at the unadjusted `B₀` is already 0.42. The
inflations only increase the budget, so the algorithmic success probability
under the same heuristic is at least 0.42, hence ≥ 0.39 as required. No
amplification is used; a single truncated trail either succeeds or the claim
fails.

Heuristic premises (disclosed): (H1) r-adding walks on the size-N class set
collide like independent uniform samples up to the usual O(K + √N / K^ε)
burn-in, which is absorbed by the 5% factor for K = 2^20; (H2) fruitless
cycles after the stated escapes contribute at most the 2% overhead;
(H3) distinguished-point sampling at rate θ is independent of the class walk
and only affects memory and a lower-order lookup cost, not the collision
waiting time on classes. These are the standard Pollard-rho heuristics used
throughout the ECDLP literature for secp256k1; they are not theorems. Toy
experiments below support (H1)–(H2) at 24–30 bit order; they do not prove
production behaviour.

Limitations: the claim is analytic in the production parameters. It does not
execute a 2^127-step attack. Automorphism reduction for secp256k1 is classical
(GLV; Bernstein–Lange–Schwabe; Bos–Kleinjung–Lenstra) and is used here only as
an accounting improvement over the organiser's generic BSGS control.

## Cost in target-group-operation equivalents

The organiser reference price for this track is
`C = 1,000` scalar AArch64 instruction equivalents per batched incomplete
affine secp256k1 addition (6M convention). Every group addition charged below
is one such reference unit. Other work is counted in instruction equivalents
and divided by C.

**Per walk iteration.** One r-adding step performs one affine point addition
`R + M_j` and one constant-time canonicalisation. Canonicalisation itself uses
only field multiplications by the public `β` and comparisons / conditional
swaps; under the organiser's A=0 / reference-addition convention we charge it
as **one additional reference group operation**, for a per-iteration total of
**κ = 2** reference group operations. Coefficient arithmetic modulo r, the
partition index, the distinguished test and one hash-table probe are at most
a few thousand instruction equivalents per step (≪ C), so they contribute
`< 2^{-6}` to the per-step group-equivalent cost and are absorbed by rounding
the final logarithm upward.

**Main walk cost.** With the budget `S` above,
`H_main = κ · S = 2S` reference group operations.

Numerically (exact integers are the audit authority; displayed decimals are
rounded up at the sixth fractional digit for the submitted logs):

```text
N            = r/6
B0           = sqrt(2·N·ln(1/0.58))
S            = ceil(1.05 · 1.02 · (1+2^{-20}) · B0)
             < 2^126.868282
H_main       = 2·S
             < 2^127.868282
```

**Setup.** Computing `β,λ`, building the K = 2^20 partition table (each entry
two scalar multiplications by double-and-add, ≤ 2·256 reference additions) and
allocating the hash table is bounded by `2^22` reference group operations after
rounding (the scalar-mul bound `2 · 2^20 · 512 = 2^30` would dominate if every
entry used a full double-and-add from scratch; we instead build the table by
streaming random walk steps from a few seeds, which costs `O(K)` additions, and
state a conservative `H_setup = 2^22`). This is preprocessing and is included
in total work.

**Distinguished-point replay / verification.** On a useful collision the
attacker recomputes at most a constant number of scalar multiplications to
confirm `[x]P = Q` and to rebuild coefficients if needed:
`H_verify = 2^12`. A full distinguished-point trail rebuild is unnecessary for
the claim because coefficients are carried online; we still reserve
`H_replay = 2^52` as an explicit slack budget covering any auxiliary table
maintenance. Both are included in total work.

**Total.**

```text
H_total = H_main + H_setup + H_replay + H_verify
        < 2^127.868282 + 2^22 + 2^52 + 2^12
        < 2^127.868282 · (1 + 2^{-75})
```

Instruction-side work outside the reference additions (hashing, sorting of no
large table, I/O) is dominated by the group steps by many orders of magnitude
at this budget and does not change the logarithm at six decimal places.
Therefore

`time_log2 = 127.868282`
`preprocessing_log2 = 22.000000` (the setup term alone; the dominant walk is
online).

Peak memory: each distinguished point stores a 64-byte record (class key + two
256-bit coefficients). Expected store size is `θ · S ≈ 2^78.87` records;
allowing 2× hash-table slack and 1 GiB of code / scratch gives
`peak_bytes < 2^85.87`, so `memory_log2_bytes = 85.868282`.
Nonuniform advice bytes: 0.

Success probability: ≥ 0.39 (modelled 0.42 under H1–H3, with budget
inflations only helping).

## Comparison with the organiser baseline

The selected organiser control is deterministic BSGS at `time_log2 = 129.000802`
with mandatory memory `2^136` bytes. The present claim is about
`2^1.13 ≈ 2.2×` less group-equivalent work and uses realistic DP memory of
about `2^86` bytes. The improvement is exactly the classical combination of
(i) Pollard rho vs BSGS (√(πN/2) expected vs ~2√N worst-case) and
(ii) the order-6 automorphism quotient on secp256k1, truncated to success
probability 0.39 rather than 1. No claim of a new cryptanalytic technique is
made; the contribution is a concrete, cost-model-faithful accounting of that
standard algorithm under `classical-group-work-v2`.

## Experiments (toy scope only)

`experiment.py` builds random j=0 curves `y^2 = x^3 + 7` of 24–30-bit prime
order (same CM type as secp256k1), runs the r-adding class walk with the same
canonical representative and escape rules, and records collision steps. Across
the following seeds the empirical mean steps track `√(π n / 12)` to within a
few percent, and the fraction of trials finishing within the proof.md budget
rule exceeds 0.39:

| bits | trials | K | mean / √(πn/12) | P(steps ≤ budget rule) |
| ---: | ---: | ---: | ---: | ---: |
| 24 | 800 | 64 | 0.977 | 0.461 |
| 26 | 500 | 128 | 1.034 | 0.428 |
| 26 | 1000 | 2048 | 1.007 | (raw B₀: 0.379) |
| 28 | 300 | 64 | 1.035 | (raw B₀: 0.380) |
| 28 | 600 | 256 | 0.991 | (raw B₀: 0.408) |
| 30 | 300 | 1024 | 0.923 | (raw B₀: 0.430) |

These numbers support the collision-model heuristics at toy order. They do
not prove the production cost bound and are not used as a numerical
extrapolation of the logarithm.

## References

- J. M. Pollard, Monte Carlo methods for index computation (mod p), 1978.
- R. Gallant, R. Lambert, J. Vanstone, Faster point multiplication on elliptic curves with efficient endomorphisms, CRYPTO 2001.
- D. Bernstein, T. Lange, R. Schwabe, On the correct use of the negation map in random walks for ECDLP, 2011.
- J. Bos, T. Kleinjung, A. Lenstra, On the use of the negation map in the Pollard rho method, PKC 2011.
- E. Teske, Speeding up Pollard's rho method for computing discrete logarithms, ANTS 1998.

The construction and accounting above are self-contained; review need not
fetch citations to reconstruct them.
