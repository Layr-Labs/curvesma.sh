# secp256k1-ecdlp: distinguished-point Pollard rho — tight accounting

Status: ready for exploratory review; no qualification or emitted score.
Same distinguished-point Pollard rho generic attack as the sibling candidate
`pollard-rho-birthday`, but with a tighter, less-padded cost bound. The
previous candidate bounded all non-group-op work at a generous `0.5 H`; this
one bounds it at `0.064 H` (at most 64 scalar instructions per walk step),
giving a strictly lower `time_log2`. The algorithm, success model and target
are unchanged; only the overhead bound is tightened. The normative target and
accounting inputs are those of the protected `target-profile.json` and
`classical-group-work-v2`; `claim.json` binds both by SHA-256. No restricted
challenge distribution, quantum operation, or uncharged advice is used.

## Construction and success

Let `p` and `r` be the protected field modulus and prime subgroup order of
secp256k1. Given `Q = [x]P` for uniform `x in [0, r)`; recover `x`. Since `r`
is prime and the cofactor is 1, the cyclic subgroup of order `r` is the whole
curve group, so a unique `x in [0, r)` exists.

Run a distinguished-point Pollard rho walk. Represent each state as the log
vector `(a, b)` with group element `g^a h^b`, where `g = P`, `h = Q`. The
walk starts at `(0,0)` (the identity). Each step multiplies the current point
by `g` (set `a += 1`), by `h` (set `b += 1`), or by `g*h` (set `a += 1,
b += 1`), chosen by a deterministic sector function of the point's
`x`-coordinate (the standard random-walk mixing). A step is distinguished when
the top `D = 16` bits of `x` are zero (density `2^-16`); distinguished states
are stored in a hash table keyed by the point, recording `(a, b)`. On a
collision of two distinguished states `R1 = R2` with log vectors `(a1,b1)`,
`(a2,b2)`:

    g^(a1-a2) = h^(b2-b1) = g^(x (b2-b1))
    =>  x = (a1 - a2) * (b2 - b1)^(-1)  (mod r),

then verify `[x]P = Q` and output `x`; a failed verification or a zero
denominator restarts the walk (charged in full).

After `N` total steps the birthday collision probability is
`P(N) = 1 - exp(-N^2 / (2r))`. Taking `N = ceil(sqrt(r))` gives
`P = 1 - exp(-1/2) = 0.39347 >= 0.39`. Since `r < 2^256`, `sqrt(r) < 2^128`,
so `N <= 2^128`.

## Cost accounting (classical-group-work-v2, C = 1000)

One reference target operation is a batched incomplete affine secp256k1
addition (6M convention), fixed price `C = 1000` scalar AArch64
instruction-equivalents. Total `W = H + R/C`. We do not charge both a
reference operation and its internal instructions.

**Reference operations (`H`).**
- The walk: `N <= 2^128` affine additions, one reference unit each.
- Verification + setup: one `[x]P` check is at most `2*256 = 512` affine
  additions, plus decoding `Q` / serialising `x` and the curve-membership
  check: at most `2^12` reference units total.
- So `H = N + 2^12 < 2^128 + 2^12`.

**Scalar overhead (`R`, in instructions) — the tightened term.**
- Per-walk-step log-vector update: up to two modular additions mod `r`
  (a `256`-bit add + compare + conditional subtract is ~10 instructions;
  the worst sector updates both `a` and `b`), amortised over the three
  sectors ~19 instructions/step.
- Per-walk-step sector selection: a cheap hash of the `x`-coordinate (a few
  limb operations), ~8 instructions.
- Per-walk-step distinguished test: a 16-bit zero check over the top limbs
  of `x`, ~8 instructions.
- Distinguished-point table: `N/2^16` inserts at ~32 instructions each,
  amortised to `N * 2^-10` instructions, negligible.
- Randomness/amplification: the walk is driven deterministically by its
  sector function (no external coins); failed attempts are absorbed by the
  restart charge above.
- Verification/output/setup instructions: at most `2^14`.

Bounding every per-walk-step scalar cost generously at **64 instructions/step**
(round number, well above the ~19+8+8 ~= 35 actual), the walk overhead is
`R_walk <= 64 * N` instructions, i.e. `R_walk/C = 64*N/1000 = 0.064*N`
reference-equivalents. Adding the negligible table term and the `2^14`
verification/output/setup instructions (which are `< 2^12/C`-order and
dominated), the total non-reference work satisfies `R/C <= 0.064*N + 2^12`.

**Total.**
    W = H + R/C = (N + 2^12) + (0.064*N + 2^12)
        = 1.064*N + 2^13
        < 1.064*2^128 + 2^13
        < 1.064*2^128          (2^13 is lower-order vs 2^128).

So `W < 1.064 * 2^128 = 2^(128 + log2 1.064) = 2^128.0895`, hence

    time_log2 < 128.0895.

We submit the clean safe upper bound `time_log2 = 1281/10 = 128.1`
(`2^128.1 = 1.069 * 2^128 >= 1.064 * 2^128`). This is strictly lower than the
sibling candidate's `128.6` (which padded the overhead at `0.5 H`) and lower
than the BSGS control's `129.000802`.

## Memory, preprocessing, advice

Peak storage is the distinguished-point table: `N_dp = N/2^16 = 2^112`
entries, each storing the point (`2*32` bytes), the `(a,b)` log vector
(`2*32` bytes) and bookkeeping, at most `2^7` bytes. Table bytes `<= 2^112 *
2^7 = 2^119`; a `2^20` reserve keeps the peak `<= 2^119` bytes.
`memory_log2_bytes = 119`.

Preprocessing is charged in full: the walk has no off-line precomputation;
only input decoding and constant loading, at most `2^12` reference
equivalents, included in the total. `preprocessing_log2 = 12`.

Non-uniform advice is `0` bytes. `nonuniform_advice_log2_bytes = 0`.

## Heuristics (disclosed)

- **Pollard rho birthday collision:** walk states behave as independent uniform
  group elements; collisions follow `1 - exp(-N^2/(2r))`. Standard,
  well-supported basis of the generic `O(sqrt(r))` DLP estimate (Pollard
  1978/1983; van Oorschot & Wiener 1999). Not an exact-constant proof;
  bounded conservatively.
- **Random-walk mixing / sector choice:** the deterministic sector function
  approximates a uniform random walk; a degenerate choice only lengthens the
  walk, absorbed by the bound.
- **Distinguished-point density `2^-16`:** sets memory only, not the work
  bound.
- **Per-step instruction bound (64/step):** a round, generous upper bound on
  the log-vector update + sector + distinguished test; the tight estimate is
  ~35 instructions/step, so 64 carries margin. This is the term tightened
  relative to the sibling candidate.

## Comparison

| candidate | overhead bound | time_log2 | memory_log2 |
| --- | ---: | ---: | ---: |
| BSGS control | `2*sqrt(r)` + `2^136` mem | `129.000802` | `136` |
| `pollard-rho-birthday` | `0.5 H` | `128.6` | `119` |
| **`pollard-rho-tight`** | `0.064 H` | **`128.1`** | `119` |

A tighter, still-supported bound on the same attack. Beating the BSGS
control (and the sibling candidate) does not by itself establish novel
cryptanalysis; this is the expected `sqrt(r)` generic frontier. A GLV
endomorphism refinement (secp256k1 has CM) could tighten the constant
further but is not claimed here.

References: [Pollard's rho for log problems](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm_for_logarithms);
[van Oorschot & Wiener 1999](https://eprint.iacr.org/1999/068.pdf);
[HAC ch. 3.6](https://cacr.uwaterloo.ca/hac/about/chap3.pdf).

## Submitted computation bound

| Component | Bound (reference-equivalents) |
| --- | ---: |
| walk `H_walk` | `N <= 2^128` |
| verification + setup | `<= 2^12` |
| per-step scalar `R_walk/C` | `<= 0.064 N` (64 instr/step) |
| table `R_table/C` | `<= N * 2^-19.97` (negligible) |
| **total `W`** | `< 1.064 * 2^128 = 2^128.0895` |

Submitted `time_log2 = 1281/10 = 128.1` (safe upper bound above `128.0895`).
`memory_log2_bytes = 119`; `preprocessing_log2 = 12`;
`nonuniform_advice_log2_bytes = 0`; `success_probability = 0.3935 >= 0.39`.
All table access, failures, verification, output and amplification are
included; the fixed reference price `C = 1000` sets units, not the work.
