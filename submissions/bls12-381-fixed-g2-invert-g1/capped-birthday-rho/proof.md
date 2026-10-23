# Capped birthday rho in GT for fixed-G2 pairing inversion into G1

## Instance

The challenge supplies a nonidentity `Q` in G2 and `z` in GT. Under the pinned
profile, `Q = [b]Q0` with `b` uniform in `1..r-1`, and `z = e(P0,Q)^a` with `a`
uniform in `0..r-1`. The required output is `P` in G1 with `e(P,Q)=z`.

## Exact reduction to a target-group discrete logarithm

Let `g = e(P0,Q0)`, a generator of GT of prime order `r`. Since `Q = [b]Q0`
and the pinned pairing is bilinear,

    h := e(P0,Q) = e(P0,Q0)^b = g^b.

Because `b` is nonzero modulo the prime `r`, `h` is again a generator of GT.
The challenge value satisfies

    z = e(P0,Q)^a = h^a.

Therefore `a = log_h z`, and `P = [a]P0` gives
`e(P,Q) = e(P0,Q)^a = h^a = z`. The reduction is exact; it uses only
bilinearity and nondegeneracy of the pinned pairing, and it works for every
nonidentity supplied `Q`, not only for `Q0`.

The reduction costs one pairing to form `h`. A BLS12-381 pairing (Miller loop
plus final exponentiation) is conservatively charged at `2^12` dense
Fp12-multiplication equivalents, which is well above the pinned measurement.

## Group and reference unit

GT is cyclic of prime order

    r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001,
    log2(r) = 254.8570894130472.

The track reference operation is the dense Fp12 multiplication in blst, with
`C = 20000` scalar AArch64 instruction equivalents. Every GT walk step below is
a general Fp12 multiplication, so it costs exactly one reference operation. No
cheaper cyclotomic multiplication or squaring is claimed anywhere in the
ledger.

## Walk

Use Teske's `r`-adding walk with `L = 32` classes. Write every walk point as
`x_i = h^{α_i} z^{β_i}`. A fixed partition function reads the low five bits of
the first canonical coefficient of `x_i` and selects a precomputed multiplier
`m_j = h^{a_j} z^{b_j}` with small `a_j, b_j`. Then

    x_{i+1} = x_i * m_j,   α_{i+1} = α_i + a_j,   β_{i+1} = β_i + b_j  (mod r).

The starting pair `(α_0, β_0)` is drawn uniformly. The partition is a fixed
hash, so the walk is deterministic given the instance and the starting pair.

Each step is one dense Fp12 multiplication, i.e. one reference operation.
Bookkeeping is two 256-bit additions with conditional reduction plus the
partition hash; it is charged at 64 scalar instructions per step, which is
deliberately conservative.

## Distinguished points and relation extraction

A point is distinguished when the low 84 bits of its first coefficient vanish,
so the distinguishing probability is `θ = 2^-84`. Every distinguished point is
stored as `(x, α, β)` in a hash table. If a distinguished point repeats, then
`x_i = x_j` with `i < j`, hence

    h^{α_i-α_j} = z^{β_j-β_i} = h^{a(β_j-β_i)},

so

    a ≡ (α_i - α_j) * (β_j - β_i)^{-1}  (mod r)

whenever `β_i ≠ β_j`. The computed `a` is checked by `h^a = z` and by the final
pairing. The probability that a collision has `β_i = β_j` is `~1/r`, far below
the success slack; a detected failure costs one restart, which the bound below
absorbs.

## Success probability

Until its first repeat the walk visits distinct group elements. Under the
standard random-map model the probability of no repeat in `t` steps is

    Π_{i<t} (1 - i/r) ≈ exp(-t^2 / (2r)),

which is the usual Pollard-rho birthday heuristic (Teske 2001; Brent 1980;
van Oorschot and Wiener 1999). With the cap

    B = ceil(sqrt(r)),

the success probability is

    P(T ≤ B) ≈ 1 - exp(-1/2) = 0.3934693403... ≥ 0.39.

The algorithm runs exactly `B + 2^84` steps and declares success if and only if
the first repeat occurs by step `B`; the extra `2^84` steps only allow the
repeated distinguished point to be observed. Since `2^84` is negligible
relative to `B`, the cap is effectively `B`.

At the first repeat the walk has entered a cycle of length `λ ≈ sqrt(πr/8) =
2^126.75`, which contains `θλ ≈ 2^42.75` distinguished points in expectation,
so the repeated distinguished point is found within `2^84` further steps with
overwhelming probability. Both quantities are far below `B`.

## Cost ledger

- Walk steps: `B + 2^84 = 2^127.428545 * (1 + 2^-43.4)`.
- Per-step reference operations: `1` dense multiplication plus `64/20000 =
  0.0032` for bookkeeping, i.e. `1.0032`.
- Walk total: `2^127.428545 * 1.0032 = 2^127.433154` reference operations.
- Preprocessing: one pairing (`2^12`) plus the 32-multiplier table (`2^6`);
  reported as `2^13.0` and included in total work.
- Distinguished-point operations: `θB = 2^43.43` stores/lookups at at most 100
  scalar instructions each, i.e. `2^50.1` instructions `= 2^35.8` reference
  operations.
- Output: one G1 scalar multiplication (`2^6`) and one verification pairing
  (`2^12`).
- Sum: `W = 2^127.433154 + 2^35.8 + 2^13 + 2^12.1 < 2^127.434`.

The declared score is `time_log2 = 127.434`, rounded up from `127.433154`, so
it is an upper bound on the charged total work.

## Memory

Peak storage is the distinguished-point table: `θB = 2^43.43` entries, each
640 bytes (12 Fp coefficients of 48 bytes plus two 32-byte exponents), with a
2× hash-table allowance, giving `2^43.43 * 1280 = 2^53.75` bytes. The declared
peak is `2^54.0` bytes, including code, walker state and the multiplier table.
Memory is supplementary and has no scalar contribution or tie-break.

## Heuristics, uncertainty and scope

The only heuristic is the standard random-map/birthday model for the rho walk.
It is the accepted exploratory standard and is supported by the cited
literature; no proof of ideal randomness is claimed. The bound is analytic; no
full-size run was performed and none is feasible. The attack is classical, uses
no uncharged advice, no auxiliary powers and no quantum steps. It does not
claim novelty, a practical break or a security estimate. The declared memory is
required by the distinguished-point method.
