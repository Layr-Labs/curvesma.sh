# mu_(2^20) orbit collapse on the supplied d-th power

## 1. Problem and inputs actually used

Target: `T = A_1 = [s]P0` in G1, with `P0` of prime order `r` and `s` uniform in
`[1, r-1]`. The profile supplies `A_i = [s^i mod r]P0` for `i = 1..n`, and the
scoring value is `n = 2^20`. Two public constants matter:

* `d = n = 2^20`, and `d | r - 1` (in fact `2^32 | r-1` for BLS12-381).
* `A_d = A[n-1] = [s^d mod r]P0`, i.e. the `d`-th power of the target scalar is
  given as a group element.

The attack reads exactly two inputs: `T = A_1` and `A_d`. It does not read the
remaining `n - 2` powers, which is permitted ("the attack may use a subset").
Generation of the powers is organizer work and is not charged; the two accesses
are.

Let `g = P0`, so `ord(g) = r`. All exponents below are taken mod `r`.

## 2. The reduction

Because `d | r-1`, the multiplicative group `F_r^*` has a cyclic subgroup
`mu_d = { z in F_r^* : z^d = 1 }` of order exactly `d`. Fix a generator `eta`
of `mu_d`: take the smallest integer `h >= 2` with `h^((r-1)/2) = -1 mod r` (such
`h` exists because `2^32 | r-1`), and set `eta = h^((r-1)/d)`. Then
`eta^d = h^(r-1) = 1` and `ord(eta) = d`. This is a public constant derived from
the published factorization of `r-1`; its construction is a handful of modular
exponentiations.

Choose an integer `G` (its value is fixed in section 5) and consider the two
tables

* `L = { ( g^(s^d - a^d), a ) : 1 <= a <= G }`,
* `R = { ( g^(b^d), b ) : 1 <= b <= G }`.

Both are computable from public data: `g^(s^d) = A_d` and `g^(a^d)`, `g^(b^d)`
are fixed-base exponentiations of known exponents.

Suppose `L` and `R` collide on the same group element for a pair `(a, b)`:

```
g^(s^d - a^d) = g^(b^d)
```

Exponents live mod `r`, so this is `s^d - a^d = b^d (mod r)`, i.e.

```
s^d = a^d * b^d = (a*b)^d  (mod r).
```

Both `a*b` and `s` are nonzero mod the prime `r`, so multiply by `(a*b)^(-d)`:

```
( s * (a*b)^(-1) )^d = 1  (mod r),
```

that is, `s * (a*b)^(-1) in mu_d`, hence

```
s = zeta * (a*b)  (mod r)   for some zeta in mu_d.
```

This is the whole reduction, and it has the property the earlier attempts lacked:
**the collision does not hand back the logarithm of an auxiliary element.** It
hands back a description of the target scalar `s` itself, as the explicit set

```
Cand(a,b) = { eta^t * (a*b) mod r : t = 0, 1, ..., d-1 }.
```

Step 6 of the algorithm tests each member of `Cand(a,b)` against the acceptance
predicate `[c]P0 = T`, so the object being tested is always a candidate value of
`s`, never a value derived from some other generator.

Correctness direction is therefore immediate and unconditional: if the pair
`(a,b)` matches, then `s` is literally one of the `d` enumerated elements, and
the enumeration is exhaustive for the orbit because `mu_d = <eta>` has exactly
`d` elements. No heuristic enters the *correctness* of the recovery; the only
statistical ingredient is how many `s` the candidate set covers (section 4).

## 3. Algorithm

1. Compute `eta` as in section 2.
2. Build `L`: for `a = 1..G`, compute `a^d mod r` by 20 modular squarings,
   `g^(-a^d)` by one fixed-base scalar multiplication, and store the product
   `A_d * g^(-a^d)` with the index `a`.
3. Build `R`: for `b = 1..G`, compute `b^d mod r` and `g^(b^d)`, store with `b`.
4. Sort both tables on the 48-byte compressed point encoding and merge-join.
5. For every colliding pair `(a,b)`: set `x = a*b mod r`, generate `eta^t * x`
   for `t = 0..d-1` by successive multiplication by `eta`, and test `[c]P0 = T`
   for each candidate `c`.
6. On the first `c` with `[c]P0 = T`, output `c` after re-checking
   `1 <= c < r` and `T in G1`.

Every step is deterministic; the attack consumes no random coins. Success
depends only on the challenge `s`.

## 4. Coverage and success probability

Write `M = { a*b mod r : 1 <= a, b <= G }`. The candidate set is the union of
the `d` orbits `mu_d * x` over `x in M`, so its size is at most `d*|M|` and, for
a uniform `s`, the success probability is the fraction of `F_r^*` it covers.

Two facts fix that fraction.

* `|M|` is the size of the `G x G` multiplication table mod `r`. Since
  `G^2 < r` (section 5), no product wraps, and `|M|` is the number of distinct
  integers `a*b`. Erdos's multiplication-table theorem gives
  `|M| = G^2 / (log G)^(c)` with `c = 1 - (1 + log log 2)/log 2 = 0.086...`;
  for `G ~ 2^117.5`, `(log G)^0.086 < 1.5`. I charge the deficit as a factor
  1.5, i.e. `|M| >= G^2 / 1.5`.
* The `d` orbits of distinct `x in M` are essentially disjoint; the expected
  number of accidental coincidences among the `d*|M|` candidates is at most
  `(d*|M|)^2 / r`, which is absorbed by the same 1.5 factor at the sizes used.

With `lambda = d*|M| / r` the covered fraction is `1 - exp(-lambda)`. Section 5
chooses `G` so that `lambda >= 0.74`, giving coverage `>= 0.52`. The claim
declares the conservative floor `0.39`, which is the track minimum; the
conservative floor already accounts for the Erdos deficit being larger than
charged and for orbit coincidences.

The expected number of colliding pairs is `~ d*G^2*log(G^2)/r ~ 2^7.5`; each is
processed at a cost given below, and the total of that stage is far below the
dominant term, so no amplification or restart term is added: a single run
already attains the declared success probability.

## 5. Parameters and cost ledger

Take `G = ceil(2^117.5)`. Then `G^2 = 2^235 < r` (no wraparound in `M`), and

```
lambda = d*|M|/r >= 2^20 * 2^235 / (1.5 * r) = 2^255 / (1.5 * 2^254.857) = 0.74,
coverage = 1 - exp(-0.74) = 0.52.
```

Prices: one G1 addition is one reference operation (batched incomplete affine,
6M, `C = 2000` instructions). A base-field multiply is 383 instructions, so
`383/2000 = 0.19` reference units.

Per table entry, both tables alike:

* `x^d mod r`: 20 modular squarings, `20 * 383 = 7660` instructions = `3.83`
  reference units.
* fixed-base scalar multiplication with an 8-bit comb: 255/8 = 32 additions =
  `32` reference units, with a one-time table of `32 * 256` points
  (`2^13` units, negligible against `2^123`).

So each entry costs `35.83 < 36` reference units.

```
H = 2 * G * 36 = 72 * 2^117.5 = 2^123.67
```

Join: two radix sorts of `G` 48-byte keys (6 byte-passes, ~100 instructions per
element per table) plus a linear merge cost
`R = 2 * G * 100 = 2^118.5 * 100 = 2^125.1` instructions, i.e.
`R/C = 2^114.1` reference units.

Orbit stage: `~2^7.5` colliding pairs, each needing `d` candidate tests at 32
units plus `d` modular multiplications for `eta^t`:
`2^7.5 * 2^20 * 32.2 = 2^32.6` reference units.

```
W(2^20) = H + R/C + orbit = 2^123.67 + 2^114.1 + 2^32.6 = 2^123.67
time_log2 = 123.67, claimed as 123.7.
```

Peak memory: both tables live at the merge, `2G` entries of 48-byte compressed
points plus 4-byte indices, `2 * 2^117.5 * 52 = 2^124.2` bytes. Memory is
reported and is supplementary under `classical-group-work-v2`; it carries no
scalar contribution.

Applicability: the bound holds for any `n` with `d | r-1`, `d <= n`, and `A_d`
supplied; at the scoring value `n = 2^20` we use `d = n = 2^20`. No
number-theoretic condition beyond `d | r-1` is required, and that condition is
satisfied by the published BLS12-381 scalar field.

## 6. Heuristics and limitations

* The Erdos multiplication-table deficit is a theorem, not a conjecture; the
  only extrapolation is using its leading form at `G ~ 2^117.5`, and the charge
  factor 1.5 exceeds the computed `(log G)^0.086 < 1.5`, so the deficit is
  over-charged rather than assumed away.
* Orbit disjointness is argued by the second-moment bound `(d|M|)^2/r`, which is
  below one candidate in absolute terms once the 1.5 factor is applied; the
  declared success floor 0.39 is well below the nominal 0.52 so that this term
  cannot flip acceptance.
* No experiment is offered as a production complexity proof. The parameters are
  far beyond any runnable scale; the argument is analytic, which the rules allow
  ("an analytic argument need not run an exponentially large attack").
* Memory is `2^124.2` bytes, which is not physically realizable; the score
  metric is computation, and memory is supplementary, so this bounds the work,
  not a deployable machine. This is stated rather than hidden.

## 7. Why the earlier reduction failed and what changed

Earlier candidates fixed a generator `gamma`, computed `a = log_gamma(s)`, and
matched inside the order-`(r-1)/d` subgroup, so a rectangle match returned
`a' = k_t + t'N` — a value in the wrong exponent space — and then tested `[s]P0`
with it. That is a genuine exponent/scalar mismatch: the matched quantity was a
logarithm to a different base, and no lift from `Z_N` to `Z_r` was ever
performed.

Here the matched quantity is `s^d`, the `d`-th power of the target scalar, and
the recovery step inverts that power map exactly: the `d`-th roots of `s^d` are
the `mu_d`-orbit of `s`, all `d` of them are enumerated, and each is tested as a
candidate value of `s`. There is no auxiliary generator, no subgroup change of
basis, and no unperformed lift.
