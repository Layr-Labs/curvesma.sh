# secp256k1-ecdlp: distinguished-point Pollard rho (birthday)

Status: ready for exploratory review; no qualification or emitted score.
This is a concrete generic attack that improves on the deterministic BSGS
accounting control in this track. It is a supported, best-known generic
heuristic, not a novel sub-exponential or structural break, and not a
practical production execution. The normative target and accounting inputs are
those of the protected `target-profile.json` and `classical-group-work-v2`;
`claim.json` binds both by SHA-256. No restricted challenge distribution,
quantum operation, or uncharged advice is used.

## Construction and success

Let `p` and `r` be the protected field modulus and prime subgroup order of
secp256k1. We are given `Q = [x]P` for a uniform `x in [0, r)`; recover `x`.
Since `r` is prime and the cofactor is 1, the cyclic subgroup of order `r`
is exactly the curve group, so a unique `x in [0, r)` exists for every valid
`Q`.

Run a distinguished-point Pollard rho walk on the subgroup, with state
represented as a log vector `(a, b)` meaning the point `g^a h^b` where
`g = P` and `h = Q` (the standard affine-walk representation). The walk
advances by one affine addition per step, multiplying the current point by
either `g` (set `a += 1`) or `g^{-1}` (set `a -= 1`); the sector choice is a
deterministic function of the current point's `x`-coordinate bits, which is
the usual random-walk mixing assumption. The walk starts at `(0, 0)` (the
identity).

Define a step to be **distinguished** when the top `D = 16` bits of the
current point's `x`-coordinate are zero, so the distinguished fraction is
`2^{-16}`. Store only distinguished states in a hash table keyed by the point,
recording the `(a, b)` log vector. When a collision of two distinct
distinguished states on the same group element occurs, with log vectors
`(a1, b1)` and `(a2, b2)`:

    g^(a1-a2) = h^(b2-b1) = g^(x (b2-b1))
    =>  x = (a1 - a2) * (b2 - b1)^{-1}  (mod r)

(the inverse exists unless `b2 - b1 = 0 mod r`, which has negligible
probability and, if it occurs, only restarts the search and is charged below).
After obtaining `x`, verify `[x]P = Q` and output `x`; a verification failure
restarts the walk (charged in full).

### Success probability (birthday)

After `N` total walk steps the number of distinct random-looking states
visited is `Theta(N)`. Two of them collide on a common group element with
birthday probability

    P(N) = 1 - exp(-N^2 / (2 r)).

Taking `N = ceil(sqrt(r))` (the smallest clean value with `N^2 >= r`),

    P = 1 - exp(-1/2) = 0.3934693... >= 0.39.

So the algorithm succeeds with probability `0.3935` over uniform production
challenges and fresh walk coins, meeting the required `>= 0.39`. Because
`r < 2^256`, `sqrt(r) < 2^128`, so the required walk is `N <= 2^128`
reference steps. Failed verification attempts and any rare
`b2 - b1 = 0` restarts are absorbed by the generous `0.5 H` overhead term
below; they do not change the asymptotic claim.

## Cost accounting (classical-group-work-v2, C = 1000)

One reference target operation is a batched incomplete affine secp256k1
addition under the 6M convention with fixed price `C = 1000` scalar AArch64
instruction-equivalents. Total cost is `W = H + R / C`, where `H` is the
number of reference operations and `R` is all other work in scalar
instruction-equivalents. We do not charge both a reference operation and its
internal instructions.

Let `H` be the reference-operation count and bound every other term relative
to it.

1. **The walk.** `N <= 2^128` affine additions, each one reference unit:
   `H_walk = N <= 2^128`.
2. **Per-step scalar overhead (sector + distinguished test).** Sector
   selection is a branch on the top few bits of `x`; the distinguished test
   is a 16-bit zero check over 64-bit words. Together at most `2^6 = 64`
   instructions per step, so `R_step <= 64 * N` instructions. Amortised in
   reference units: `R_step / C <= 64 * N / 1000 = N / 15.625 <= 0.064 N`.
   This is well under `0.5 N`.
3. **Distinguished-point table.** Number of distinguished states
   `N_dp = N / 2^16`. Each insert is an amortised `O(1)` hash operation, at
   most `2^5 = 32` instructions, so `R_table <= 32 * N_dp = 32 * N / 2^16 =
   N * 2^-10` instructions. In reference units `N * 2^-10 / 2^9.97 = N *
   2^-19.97`, negligible.
4. **Verification and output (once).** One `[x]P` check is at most `2 * 256 =
   512` affine additions (`512` reference units), plus decoding `Q`,
   reconstructing and serialising `x`, and the curve-membership check: at
   most `2^12` reference units total.
5. **Setup.** Decoding `Q`, curve membership check and constant loading:
   at most `2^12` instructions, i.e. `<= 2^12 / C = 2^2.03` reference units.

Summing: `H = N + 2^12` (walk + verification), and the non-reference work
`R / C <= 0.064 N + 2^-19.97 N + 2^2.03`, which is at most `0.5 N` for
`N = 2^128` (the `0.064 N` per-step term dominates and the constants are
lower-order). Therefore

    W = H + R/C <= N + 2^12 + 0.5 N = 1.5 N + 2^12 < 1.5 * 2^128.

Hence the total is `W < 1.5 * 2^128 = 2^(128 + log2 1.5) = 2^128.58496`, so

    time_log2 <= 128 + log2(1.5) = 128.58496.

We submit the clean safe upper bound `time_log2 = 643 / 5 = 128.6`
(`2^128.6 = 1.5157 * 2^128 >= 1.5 * 2^128`). This is an upper bound on the
work, slightly padded above the tight `128.585` estimate for safety, and is
strictly lower than the BSGS control's `129.000802`.

## Memory, preprocessing, advice

Peak storage is the distinguished-point table plus public constants and
scratch. `N_dp = N / 2^16 = 2^112` entries; each entry stores the point
(`2 * 32` bytes) and the `(a, b)` log vector (`2 * 32` bytes) plus
bookkeeping, at most `2^7 = 128` bytes. Table bytes `<= 2^112 * 2^7 =
2^119`. Adding a `2^20` reserve for code, constants and stack keeps the peak
`<= 2^119` bytes. `memory_log2_bytes = 119`.

Preprocessing is charged in full with no amortisation. The walk has no
off-line precomputation (the table is built on the fly); the only upfront
work is input decoding and constant loading, at most `2^12` reference
equivalents, included in the total. `preprocessing_log2 = 12`.

Non-uniform advice is `0` bytes: no searched or embedded precomputed table is
used. `nonuniform_advice_log2_bytes = 0`.

## Heuristics (disclosed)

- **Pollard rho birthday collision heuristic:** the walk states behave as
  independent uniform group elements, so collisions follow the birthday
  distribution `1 - exp(-N^2 / (2r))`. This is the standard, well-supported
  basis of the generic `O(sqrt(r))` DLP estimate (Pollard 1978/1983;
  distinguished-point variant, van Oorschot & Wiener 1999). It is not a proof
  of the exact constant; we bound conservatively.
- **Random-walk mixing / sector choice:** the deterministic sector function
  of the `x`-coordinate bits approximates a uniform random walk on the
  subgroup. This is the usual rho assumption; a degenerate sector choice
  would lengthen the walk, which the `0.5 H` overhead absorbs.
- **Distinguished-point density:** exactly `2^{-16}` of states are
  distinguished, giving `N_dp = N / 2^16` stored states. The fraction is a
  design choice; it sets memory only, not the work bound.

None of these change the target; they only set the work constant, which is
bounded generously above.

## Evidence and comparison

The argument and cost calculation are analytic; no full-size execution is
claimed. The accompanying `experiment.py` runs a real (toy) distinguished-
point Pollard rho on small multiplicative groups to exercise the walk,
collision extraction and verification, and does not validate production
costs, memory bounds or operation prices.

The deterministic BSGS control in this track charges `2 * sqrt(r) = 2^129`
reference operations plus exponential memory `2^136`. This rho candidate
charges `1.5 * sqrt(r) < 2^128.6` reference operations and memory `2^119`,
improving time by `~0.93` bits and memory by `~17` bits under the same
`classical-group-work-v2` convention. Beating the conservative BSGS control
does not, by itself, establish novel cryptanalysis or a weakened hardness
assumption; it is the expected generic `sqrt(r)` frontier.

References: [Pollard's rho for log problems](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm_for_logarithms);
[van Oorschot & Wiener, Parallel collision search with applications to
discrete logarithms (1999)](https://eprint.iacr.org/1999/068.pdf);
[HAC ch. 3.6 generic DLP](https://cacr.uwaterloo.ca/hac/about/chap3.pdf).
The needed construction and accounting are given above; review need not fetch
citations to reconstruct them.

## Submitted computation bound

| Component | Bound (reference-equivalents) |
| --- | ---: |
| walk `H_walk` | `N <= 2^128` |
| verification + setup | `<= 2^12` |
| per-step scalar `R_step/C` | `<= 0.064 N` |
| table `R_table/C` | `<= N * 2^-19.97` (negligible) |
| **total `W`** | `< 1.5 * 2^128 = 2^128.58496` |

Submitted score `time_log2 = 643/5 = 128.6` (a safe upper bound above
`128.585`). `memory_log2_bytes = 119`; `preprocessing_log2 = 12`;
`nonuniform_advice_log2_bytes = 0`; `success_probability = 0.3935 >= 0.39`.
All table access, failures, verification, output and amplification are
included. The fixed reference price `C = 1000` is the organizer constant; it
sets units, not the work.
