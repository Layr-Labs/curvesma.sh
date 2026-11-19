# Fixed-G1 pairing inversion into G2: one GT logarithm by a distinguished-point rho walk on reference Fp12 multiplications

Problem: given the instance's nonidentity P in G1 and z in GT, return Q in G2
with e(P,Q)=z, where e is the protected function
`py_ecc.optimized_bls12_381.pairing(Q,P)` (library takes G2 first). We use the
admissible route named in rules.md: solve z=e(P,Q0)^x in GT and return [x]Q0.

This is a generic Pollard-rho attack. It is not a best-known finite-field
attack, a novelty claim, a practical break or a security-bit estimate. The work
bound is analytical and capped. Success relies on the standard, explicitly
disclosed r-adding random-walk heuristic. We seek the protected exploratory
`plausible_not_refuted` standard; no prior acceptance is assumed.

Relative to a Floyd baseline, two accounting choices are made, both explicit:

* Every walk step performs exactly one **reference operation**: the pinned blst
  dense Fp12 multiplication `blst_fp12_mul`, in blst's own tower
  representation. By the protected cost model a reference operation costs one
  unit. All other per-step work (canonicalization, partition and DP tests,
  scalar labels, control) is charged separately as at most E = 8,192 scalar
  instructions (machine.md), i.e. 0.4096 units.
* One walk with distinguished points (DPs) replaces Floyd's three
  multiplications per iteration. The walk is capped at
  S = floor(3*(isqrt(r)+1)/2) + 2^115 steps.

## Complete algorithm

1. Read the profile and the instance (P, z). Decode P and the 12 canonical
   coefficients of z; check P on the G1 curve, P != infinity, z != 0 and
   z^r = 1. If z = 1 return Q = infinity, which the profile makes valid
   exactly for z = 1.
2. Compute g = e(P,Q0) **once** by evaluating the protected pairing function
   on the actual P and the public Q0. Cost is in machine.md. g has order r
   (below).
3. Convert g and z into the blst tower representation. The tower is
   Fp2 = Fp[u]/(u^2+1), Fp6 = Fp2[v]/(v^3-(1+u)), Fp12 = Fp6[w]/(w^2-v).
   Its w is the protected w (w^2=v, w^6=1+u, so u=w^6-1). The change of basis
   is the fixed linear map `poly_to_tower` in rho.py: coefficient
   c0 = f_n + f_{n+6}, c1 = f_{n+6} for the term v^k w^e with n=2k+e. Every
   coefficient is then put in canonical Montgomery form (R = 2^384).
4. Sample 66 independent uniform scalars mod r by bounded 255-bit rejection
   sampling (at most 1024 draws each). Build 32 records M_j = g^c_j z^d_j and
   a start Y0 = g^a0 z^b0. The state is (Y,a,b) with Y = g^a z^b. All of
   these elements are canonical.
5. Walk. Let t_0..t_11 be the canonical Montgomery residues of the 12 Fp
   coefficients of Y in tower order. Then:
   - The index is j = (word 0 of t_0) & 31.
   - Y is a **DP** iff the low D = 110 bits of t_1 are zero.
   - Its **fingerprint** is the low 128 bits of t_2.

   One step is: Y <- canonicalize(blst_fp12_mul(Y, M_j)), a <- a+c_j mod r,
   b <- b+d_j mod r. After each step, if Y is a DP, look up its fingerprint
   in a DP table of at most 2^21 entries:
   - **Absent:** store (fingerprint, a, b). If the table is full, report
     failure.
   - **Present** with stored (A,B): if b = B mod r report failure
     (degenerate). Otherwise set x = (A-a)/(b-B) mod r and test g^x = z.
     If it holds, stop with x. If not (a fingerprint false match), replace
     the entry and continue; after 64 failed tests report failure.

   Stop with failure after S steps.
6. Return Q = [x]Q0 by affine Fp2 double-and-add (machine.md). There is one
   trial, with no hidden restart, amplification or amortization. `rho.py`
   (`rho_dp`) is the executable reference, with the exact loop order.

## Correctness of the output

g = e(P,Q0) with e exactly the protected function. P is a nonidentity point of
the order-r group G1 and Q0 generates G2. Nondegeneracy of the pairing
therefore gives that g has order r, so g generates GT. The challenge itself is
generated as z = e(P,Q0)^a with the same function, so z = g^x for a unique x
mod r. Once g^x = z has been verified, bilinearity in the G2 argument gives
e(P,[x]Q0) = g^x = z. [x]Q0 lies on the G2 curve and satisfies
[r]Q = infinity. The instance's own P is used; no scalar is returned instead
of Q.

The walk invariant Y = g^a z^b holds initially and after each step, by
exponent addition; the basis change and the canonicalization do not change
the element. If two visits have equal Y with labels (a,b) and (A,B), then
a + x b = A + x B mod r, so x = (A-a)/(b-B) whenever b != B. The final test
g^x = z rules out any wrong output, including one caused by a fingerprint
false match. Canonical residues make equality of elements equivalent to
equality of coefficient words, so the partition, DP flag and fingerprint
are functions of the group element.

## Detection, success heuristic, support and limits

The walk Y_0, Y_1, ... is a deterministic function iteration on GT. Let mu be
its tail length and lambda its cycle length, so Y_mu is the first repeated
state, revisited at step mu+lambda. Let delta be the number of steps from Y_mu
to the first DP at or after it. If lambda > delta, that DP lies on the cycle.
It is first stored at step mu+delta and revisited at step mu+lambda+delta,
where its fingerprint is found with labels from two visits lambda steps
apart. So the algorithm reaches the verification step within the cap whenever

    (i) mu+lambda <= 1.5*sqrt(r),   (ii) delta <= 2^115 < lambda.

Earlier DP fingerprint matches are either true collisions or false matches;
neither can produce a wrong output.

Heuristic for (i). For an ideal random map on r points, the first repeat
satisfies P(mu+lambda > t) ~= exp(-t^2/(2r)). The conventional r-adding
correction for 32 increments gives the rate factor 31/32. At t = 1.5*sqrt(r)
this gives

    P(i) ~= 1 - exp(-(31/32)*1.125) ~= 0.664.

This is a stated heuristic, not an exact theorem for this group and
partition. We declare **success >= 0.50**, leaving margin over the protected
minimum of 0.39. A stricter reading of the r-adding penalty would still need
an effective factor below 0.44 to drop under 0.39.

Heuristic for (ii). A DP occurs with probability 2^-110 per fresh state under
the uniform-coefficient-bit heuristic. The chance of 2^115 consecutive
distinct states with no DP is about exp(-32). For a random map, lambda is
roughly uniform on [1, mu+lambda]; given (i), P(lambda <= 2^115) is about
2^115/2^127 = 2^-12. Further losses:
- The DP table overflows only if more than 2^21 DPs occur, against a mean of
  S*2^-110 ~= 2^18. By a Chernoff bound this is below 2^-1000.
- A fingerprint false match needs two distinct DPs sharing 128 further bits:
  at most 2^42 pairs, so below 2^-86, and 64 such events are far less likely.
- Degeneracy has probability exactly 1/r (next paragraph).
- Rejection sampling fails with probability below 66*2^-1024.

The total predicted success is about 0.66, against 0.50 declared.

Degeneracy needs no random-relation assumption. Fix the true x and condition
on the 32 increment logarithms s_j = c_j + x d_j and s0 = a0 + x b0.
Conditional on these, each d_j remains uniform (one c_j realizes each s_j).
The group walk, its indices, the DP flags and fingerprints, and the first
matching pair depend only on the s_j and s0. For a true collision between
visits lambda steps apart, b - B = -sum n_j d_j with sum n_j = lambda,
where 0 < lambda < r. So some n_j is nonzero mod r, and the combination is
exactly uniform. Degeneracy therefore has probability exactly 1/r for the
actual partition.

The material heuristic is that the partition bits taken from canonical
Montgomery coefficient words behave like random bits. This applies both to
the index from t_0 and to the DP and fingerprint bits from t_1 and t_2. It
also requires that the random 32-adding walk mixes like the analyzed model at
production size. The bits are low-order bits of x*2^384 mod p for
coefficients of pseudo-random GT elements. Partitions on coordinate bits are
the standard practice for rho/DP walks (Teske, "On random walks for
Pollard's rho method", Math. Comp. 70 (2001) 809-825; van Oorschot-Wiener,
"Parallel collision search with cryptanalytic applications", J. Cryptology 12
(1999) 1-28). These are context; the reasoning used is reproduced here.

The work and memory caps hold even if the heuristic fails: the cap is
deterministic. The success bound is exploratory, not an unconditional theorem.

Evidence (`experiment.py`, development output in `development-evidence.json`):
- The exact `rho_dp` loop runs on four toy prime additive groups (orders
  ~1e4 to ~1e7), 128 trials each. The coefficient bits are modelled by
  SHA-256 of the toy element, with DP rate 2^-2 and gap 32/theta as in
  production. It reports every outcome, plus the fraction that succeeds
  within 1.5*sqrt(q) steps without the gap.
- The tower multiplication is checked equal to the protected-basis
  multiplication through the basis map, on random full-size elements.
- A real full-size GT walk of 2048 steps uses the production key on an
  order-r element h^((p^12-1)/r). It reports the partition histogram (and
  chi-square), the DP-word low-bit frequency, and the final label invariant.

These are correctness and finite-statistics evidence only; extrapolation to
2^128-step walks is part of the disclosed heuristic. The experiment does not
evaluate the py_ecc pairing or blst itself.

Finite-field (Sex)TNFS may be stronger; nothing here claims tightness.

## Work, memory and score

    S = floor(3*(isqrt(r)+1)/2) + 2^115
    W = S*(1 + 8192/20000) + 2^40/20000     (reference units)
    submitted score = ceil(10^6*log2(W))/10^6

- The S*1 term is the S reference multiplications.
- The S*8192 term is the per-step non-reference work.
- The 2^40 term is all fixed work, with every item bounded in machine.md:
  the pairing, conversions, table and start construction, randomness, DP
  table handling and verifications, output and checks.

`experiment.py` recomputes the score from the protected r; the value is
stated in submission.json.

Peak memory is below 2^28 bytes, mostly a DP table of 2^21 slots of 96
bytes. Advice is zero. Preprocessing is one-target only and charged in full.
