# bls12-381-fixed-g1-invert-g2: special tower NFS discrete logarithm in GT

Status: ready for independent review; no qualification or emitted score.
This is an analytic, heuristic best-known-attack estimate. It is not a
practical break, and nothing at production size has been run. `claim.json`
binds the protected target profile and `classical-group-work-v2` by SHA-256.

## Reduction to a discrete logarithm in GT

Let p, r, P0, Q0 and e be as in the protected profile. The input is a
nonidentity P in G1 and a target z in GT. P has order r because r is prime,
and e is nondegenerate on G1 x G2. So z0 := e(P,Q0) has order exactly r and
generates GT, the order-r subgroup of F_{p^12}^*. Every z in GT has a unique
a in [0,r) with z = z0^a. Then Q := [a]Q0 satisfies e(P,Q) = e(P,Q0)^a = z by
bilinearity. When z = 1 this gives a = 0 and Q = infinity, which the profile
accepts exactly when z = 1. The rules name this reduction as admissible.

```text
1. Validate P and z as the profile requires.
2. Compute z0 = e(P,Q0) with one pairing evaluation.
3. Shared precomputation for F_{p^12}, k = 12 = 2*6 with p = p(x) from the
   BLS12 family (x = -0xd201000000010000): run the special extended TNFS
   (SexTNFS). Represent F_{p^12} as a tower R/(h) with deg h = 2 or 6 over Q.
   Use the special-form polynomial f from p(x) (the Joux-Pierrot / SNFS
   selection) with conjugate g, sieve over the tower ring, and filter.
   Solve the sparse linear system modulo r with block Wiedemann to get
   virtual logarithms of the factor base modulo r.
4. Individual logarithms: for t in {z0, z}, apply an initial-split
   (boot) step with fresh random exponents of t (randomised by multiplying
   by known factor-base elements), then a special-q descent down to the
   factor base. This gives log(t) mod r relative to an unknown common base.
5. If log(z0) = 0 mod r, which cannot happen because z0 has order r, or if a
   descent fails, rerun step 4 with fresh coins.
   Set a = log(z) * log(z0)^{-1} mod r.
6. Compute Q = [a]Q0 and check e(P,Q) = z. Output Q, or null when a = 0.
```

The precomputation in step 3 depends only on the field F_{p^12} and the
modulus r. It is the same for every instance, but under the rules it is
charged in full to this one cold attack. Both individual logarithms reuse
the same precomputation. A ratio of logarithms cancels the unknown base,
so the factor-base logarithms only need to be consistent up to a common
scalar modulo r.

## Cost accounting

Unit: one dense F_{p^12} multiplication equals C = 20,000 scalar AArch64
instruction equivalents (profile and cost model). log2 C = 14.287712.

Precomputation (step 3). Guillevic (PKC 2020, "A short-list of
pairing-friendly curves resistant to special TNFS at the 128-bit security
level") estimates the full SexTNFS discrete-logarithm cost for BLS12-381's
4572-bit F_{p^12}. The method follows Barbulescu-Duquesne (J. Cryptology
2019, "Updating key size estimations for pairings"). It optimises the
smoothness bound, the sieving region and the tower choice, and balances
relation collection against linear algebra. The resulting estimate is about
2^126 elementary operations. Here "elementary operation" means one sieve-cell
update during relation collection, or one sparse row-times-vector term
(a multiply-add of a residue modulo the 255-bit r by a small coefficient)
during linear algebra. Later follow-up runs of that estimator, reported in
the survey literature (e.g. arXiv:2212.01855, section on Guillevic's
short-list), drop to about 2^120. This claim uses the higher 2^126 figure.

Converting to instructions:

- A sieve-cell update is an index computation, a load, an add of a scaled
  logarithm, a store and loop control. Even with lattice-sieve bucket
  bookkeeping that is at most about 16 instructions.
- A linear-algebra term reads a column index and a small signed coefficient,
  then multiply-adds into a 4-limb (256-bit) accumulator modulo r, amortising
  the reduction. That is at most about 64 instructions including memory
  traffic.
- Cofactorisation of survivors (ECM/batch smoothness) and filtering are
  already inside the estimator's balanced total as lower-order terms.

To leave margin for the estimator's unknown o(1) and for cache and transfer
overhead, this claim charges 2^8 = 256 instructions per elementary operation.
That is 4x the worst per-operation figure above. The total is
2^(126+8) = 2^134 instructions, which is 2^(134 - 14.287712) = 2^119.712288
F_{p^12} multiplication equivalents.

Individual logarithms (step 4), two targets, each with at most 4 expected
restarts. Boot and descent cost L_{p^12}[1/3, c] with c well below the
precomputation constant. Published descent estimates for this field size are
below 2^100 elementary operations. Charging 2^8 instructions each and 8 runs
gives at most 2^(100+8+3-14.29) < 2^97 units. That changes log2 W by less than
2^-21.

Pairing evaluations and the scalar multiplication [a]Q0, plus the
verification pairing, cost fewer than 2^14 units in total. This is
negligible.

Total: W <= 2^119.712288 + 2^97 + 2^14 < 2^119.7123 units. The submitted
score rounds up to time_log2 = 120 (exact rational 120/1). That leaves about
0.29 bits for any additive constant missed above.

Memory. The factor base at the estimator's optimum is about 2^60 ideals.
The sparse matrix holds about 2^60 rows of weight about 2^7, with 8-byte
entries. That is about 2^70 bytes, and the Wiedemann state is smaller.
Declared peak memory: at most 2^71 bytes. No nonuniform advice is used.
Preprocessing equals step 3 and is included in W.

## Success probability

Under the standard NFS heuristics (smoothness probabilities of the norms
follow the Canfield-Erdos-Pomerance estimate; the relation matrix has
corank small enough that the virtual logarithms modulo r are determined up
to the common scalar), the precomputation succeeds. Its failure modes are
corank defects that are detected and fixed by collecting more relations,
which is already inside the balanced estimate. Each descent is Las Vegas:
a failed boot or special-q step is detected and restarted with fresh coins,
and step 6 checks the output. With the restart budget charged above, the
success probability over uniform P, uniform z and the coins is at least
0.9, above the 0.39 minimum. The z = 1 case is handled exactly.

## Heuristics and limitations

- The NFS smoothness and linear-algebra heuristics are standard and
  supported by record computations (e.g. the 2016-2021 TNFS/NFS records in
  F_{p^6} and F_{p^12} at small size). Extrapolating them to a 4572-bit
  F_{p^12} is not experimentally validated.
- The 2^126 figure is a published estimator output with an unknown o(1);
  the 2^8-instructions-per-operation charge and the rounding to 120 are the
  margin for that. The estimate is not an upper bound proved for this
  field.
- No production run, no machine code. This is an analytic estimate under
  the exploratory review standard.
