# Argument: class walk for fixed-argument pairing inversion into G2, claim 2^(12589/100)

Track `bls12-381-fixed-g1-invert-g2`. Score is log2 of total classical work in units of
one dense Fp12 multiplication at the protected C = 20000. Success must be at least 0.39 over
uniform challenges and fresh coins. Submitted claim **125.89 = 12589/100**; itemised work
**2^125.78971**; every judgement-bounded row carried at 13/10 gives **2^125.85097**.

## 1. Reduction

P is a nonidentity point of G1 of order r and z is in the order-r subgroup of GT. Let
g = e(P, Q0). Reduced pairing non-degeneracy at prime order makes g generate GT, so
z = g^a for a unique a mod r and Q = [a]Q0 satisfies e(P,Q) = e(P,Q0)^a = z. Finding a is a
DLOG in GT. `part_a.py` checks the two properties this uses, on BLS12-381 itself: the map is
bilinear in both slots, g^r = 1, g != 1, and for a planted instance (P = [b]P0, z = e(P,Q0)^a)
the returned Q = [a]Q0 satisfies e(P,Q) = z exactly, while Q = [a+1]Q0 does not. The pairing
is evaluated by a Miller loop over L = |x| on the untwisted point followed by the final
exponentiation, in the tower fixed by `target-profile.json`.

## 2. Classes

ord_r(p) = 12 (`part_a` checks p^6 = -1, p^4+p^8+1 = 0, and no smaller power is 1), so
raising to p partitions GT into orbits of twelve and the walk lives on N = (r-1)/12 =
2^251.272 classes. Inversion is the sixth power of Frobenius, so folding loses nothing.

A class is named by the w-coefficient of a representative: with A and B the two base-field
reads of that coefficient, the twelve conjugate keys are {A z^j} U {B z^j}, j < 6, where
z = |gamma|^2 is a primitive sixth root. Because z^3 = -1 the twelve values are the signed
orbits of six absolute values, so the class key is the minimum of six |.| tests and the
position of the minimum fixes the canonical representative (kappa). The fold therefore costs
six base-field multiplications, four base-field adds, and eleven six-limb compare/moves:
**2850 instructions**, itemised in `ledger.py` and summed again in `verify_certificate.py`.

`part_a.py` checks on real BLS12-381 GT elements that this rule reproduces the keys of the
twelve explicitly powered conjugates as a set, that the twelve keys are distinct, that the
minimum key is invariant and its position covariant over the orbit, and that the min-of-six
identity holds. The keys therefore name classes, not representatives.

## 3. Step and mixture

From the canonical representative a step is one dense Fp12 multiply by a tabulated multiplier
(the table stores 2^14 multipliers with all twelve twists, so one multiplication suffices) or,
on a fraction f = 3/8 of classes selected by the routing word, one Granger-Scott cyclotomic
squaring. Squaring is legitimate on classes because it commutes with Frobenius; it is also
**injective on classes**: c^2 = c'^2 forces c = +-c', and -1 = p^6 lies in the orbit, so c and
c' are the same class. A pair of squaring steps therefore cannot merge, so the merge rate is
reduced by (1 - f^2) and the budget carries the factor 1/(1-f^2) = 1.0792 on sqrt(N). This is
charged, not assumed away: `sim.py` measures the effect at 1.0295, i.e. the charge is
conservative by 0.02 bits.

The squaring row is the exact operation count of the implemented formula (two Fp6 squarings
of six Fp2 squarings each, one v-shift, four Fp6 adds): 24 base-field multiplies and 122
base-field adds = **11632 instructions**, 0.582 of a reference multiplication. Both rows are
paid on every step in addition to the fold, because the class key of the current state decides
the next step type and cannot be deferred.

## 4. Detection, recovery, store

A two-cycle rule (key against the previous key, plus a running minimum checkpoint) costs two
compares per step and rejects the closed short patterns a deterministic walk enters; the
escape is applied to the minimum-key class, which is what makes the merging lemma hold: two
trails at a common class follow the same deterministic map, detect the same cycle and escape
to the same class. A distinguished-point density of 2^-80 gives 2^45.7 reports at the budget,
each an 80-byte record for 2^52.8 bytes with the full 48-byte class key, so a key match is a
class match. Reports are ordered; on a repeat both trails are replayed from their starts with
(a,b) exponent bookkeeping (2^80 steps, negligible against 2^125.8), aligned by the Frobenius
power that maps one endpoint to the other, and solved as
a1 + b1*x = p^j (a2 + b2*x), i.e. x = (a2 p^j - a1)/(b1 - b2 p^j) mod r, then verified against
z. Nothing is trusted: the exponent solve is rejected unless g^x = z.

Routing uses one 256-bit mixed word per class key with bit-disjoint fields: 14 bits table
index, 32 bits twist (mod 12), 80 bits distinguished test, 3 bits step type. `routing.py`
histograms 120000 synthetic production-width keys: 16369 of 16384 indices hit, 120000
distinct distinguished words, zero collapse. A 64-bit mixer would place the 80-bit
distinguished field and the index in the same 64 bits and collapse one of them.

## 5. Budget and success

S is the least integer with S^2 (1-f^2) >= 2 mu N ln(1/(1-0.392)), mu = 11/10, computed once in
`constants.py` and recomputed as an integer in `verify_certificate.py`: log2 S = 125.81064.
The random-map model then gives exactly 0.392 at the allowance, and the disclosed post-event
losses (cycle escapes after meeting, cut trails, degeneracies, sampling distance) are carried
at 3e-4, so success >= 0.392 - 0.0003 = 0.3917 > 0.39.

`sim.py` runs the same walk (injective doubling step for squaring, table-additive step for
multiplication, same routing fields) at the production constants for 4500 groups over three
moduli and two table sizes: pooled P(T_c <= S) = 0.4073 +- 0.0073 against the model's 0.3961
at the integer budget and the 0.39 requirement. A f = 0 control measures the injectivity
factor quoted above.

## 6. Ledger

    dense Fp12 multiply (reference operation)           20000   1.0000
    class fold                                          2850   0.1425
      A, B reads                                766
      four zeta products                       1532
      adds / negations                           80
      eleven compare+cmov                       264
      three hash words                           72
      distinguished test                          8
      two-cycle detector + checkpoint            64
      48-byte previous-key slot                  32
      addressing, counters, loop                 32
    cyclotomic squaring (squaring steps only)          11632   0.5816

    c_mul = 1.1425, c_sqr = 0.7241, c = (5/8) c_mul + (3/8) c_sqr = 0.98560
    W_raw = S c / C = 2^125.78971
    W_stress (every judgement-bounded fold row at 13/10) = 2^125.85097
    claim 125.89: pad over W_raw 0.11131 (>= 0.10), pad over W_stress 0.04959 (0.005..0.06)

Prices are the pinned ones: 383 instructions per base-field multiply (measured), 20000 per
dense Fp12 multiply (reference), with 20 per base-field add, 24 per six-limb compare/move and
24 per hash word derived as instruction bounds with a margin over their derivation.

## 7. What is not claimed

No production execution; the accumulation of the measured walk quality from N <= 2^14 to
2^251 is the principal uncertainty and is disclosed. Random-map behaviour, pseudorandom hash
fields and ideal coins are hypotheses, not theorems. The toy verifies the algorithm, not the
production complexity; the simulation measures the walk statistic, not the field arithmetic.
