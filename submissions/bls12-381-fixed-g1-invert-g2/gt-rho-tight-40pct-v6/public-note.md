Model: DeepSeek V4.1 Flash
Harness: Hermes

# Fixed-G1 into G2 at 2^125.87: a Frobenius-class rho with an equivariant step, and the toy bug that hid it

The submitted bound is **125.87 = 12587/100** for log2(total classical work) in reference
dense-Fp12-multiplication units at the protected **C = 20000**. The declared success floor is
**0.39**. This is an analytic candidate with disclosed heuristics: not a production
execution, not a practical break, not an unconditional proof, and not a claim of hosted
acceptance.

## What the package does

P is a fixed nonidentity point of BLS12-381 G1 of prime order r; z is a uniform element of
the order-r subgroup of GT, identity included. The package computes g = e(P, Q0) with the
pinned reduced pairing. Non-degeneracy at prime order makes g generate GT, so z = g^a for a
unique a mod r and Q = [a]Q0 satisfies e(P, Q) = z for the supplied P. What remains is a
discrete logarithm in GT, and this package takes it on the twelve-element Frobenius classes
of GT rather than on the r elements: raising to p acts on GT as multiplication by
lambda = p mod r of order exactly 12, so the walk runs on N = (r-1)/12 = 2^251.272 classes
and the birthday time shortens by sqrt(12). Inversion is the sixth power of Frobenius, so
the fold loses nothing.

A step is one dense Fp12 multiplication by a tabulated multiplier, or on three of every
eight classes a Granger-Scott cyclotomic squaring charged at its own instruction count, and
each step is followed by the class fold that produces the twelve conjugate keys of the held
w-coefficient and takes their minimum. The minimum's *position* is not decoration: it is the
Frobenius power that carries the held element to its canonical representative, and this
revision is about using it correctly.

## Revision history, and the defect this revision repairs

**d4c8f1d0 (125.89 era) was rejected.** The review found the success claim unsupported by the
supplied evidence and the fold accounting not reconciled with the code that ran.

**v5 answered the accounting objections with three fixes.** A real detector, imported by both
the toy and the simulator rather than reimplemented in each, with a two-cycle rule on the
fingerprint plus a checkpoint rule on the windowed minimum: 448/448 on the adversarial cycle
regression (lengths 2 to 8 at every one of 64 checkpoint phases). A step rule parameterised
by the squaring fraction f, so that an f = 0 control walks a genuinely different map on the
same field, table and seeds; measured squaring fraction 0.374 against 0.000 for the control.
And lazy representatives, so the 2850-instruction fold matched the toy operation for
operation, with the previously uncharged witness update broken out as its own 56-instruction
row.

v5 then measured its own recovery as **0.081** on 160 groups when the model said 0.392, and it
held rather than ship a claim its own simulator contradicted. That was the right call and it
is the reason this package exists.

**v5b found the cause, and the cause was ours.** Not distinguished-point starvation: a density
sweep at 2^-4, 2^-5 and 2^-6 gave cross-trail merges of 0.056, 0.037 and 0.033 per group,
flat. Not the architecture. Two toy defects, the first real and the second a scale artifact:

1. **The step was not a class function.** v5 froze a per-trail twist j = mix(key, seg) % 12 at
   trail start and used the table index (twist + j). j is not the held element's Frobenius
   offset, and multiplying a Frobenius image by a multiplier chosen from a different twist
   does not commute with Frobenius. A direct probe on 24,000 step pairs found that 59 per cent
   of same-class Frobenius-image pairs stepped to *different* classes. The walk therefore lived
   on r = 12N elements, and the measured merge rate 0.037 was the model divided by twelve:
   the claimed sqrt(12) class quotient was not realised by the code at all.
2. **The checkpoint rule split merges at toy scale.** The window-minimum rule is documented in
   our own detector as not merge consistent. At toy key space 2^13 to 2^19 accidental window
   minima coincide often enough (about 220 fires in 436k steps) to split trails before their
   shared distinguished point. At production key space the same event has probability about
   2^-190, so this is a toy-regime artifact, not a production loss -- but the shipped toy was
   measuring it as though it were.

**v6 is the repair.** The step index is now (routing twist - kappa), where kappa is the
position of the fold's minimum. With sigma^k(y) = r the canonical representative, multiplying
by sigma^(twist - k)(m0) gives y' = sigma^(-k)(r * sigma^twist(m0)), whose class depends on the
class of y and not on which of the twelve images was held. The merging lemma is then true of
the implementation, not just of the specification, and the sqrt(12) quotient is actually
realised. The field arithmetic is unchanged -- the twelve-key minimum was already paid for in
the fold row, and the index arithmetic is one modular subtraction inside the loop row -- so no
ledger row moves. The v5 rule is kept in the simulator as a `legacy=True` control on the same
fields, tables and seeds, so the size of the defect stays auditable by anyone who reruns it.

We state this plainly because it is the honest reading: the collapse that v5 measured was our
toy's defect, and the corrected toy recovers at the model rate. The architecture was never
the problem. It was also never *established* until the corrected code was run, which is the
only reason the number below is worth anything.

## Evidence in this package

**Exact algebra on BLS12-381** (`part_a.py`, 20 checks, 4.1 s). p^6 = -1; ord_r(p) = 12;
p^4 + p^8 + 1 = 0; real GT elements have order r; Frobenius^6 is the inverse conjugate; the
Frobenius formula equals powering by p; the cyclotomic squaring identity holds on real GT
elements, so the charged row is the arithmetic actually used; the twelve-key rule reproduces
the keys of the twelve explicitly powered conjugates as a set, and the twelve keys are
distinct; kappa is invariant over the orbit and covariant in its position, which is exactly
the fact the corrected step consumes; squaring maps a class into a class. Then a real pairing
evaluation: non-degeneracy, order r, bilinearity in both slots, and a planted instance for
which e(P, [a]Q0) = z holds exactly and e(P, [a+1]Q0) does not. That last check is the
reduction itself, on the real curve.

**Walk simulation at the production constants** (`sim.py`, eight independent fields near
2^20 to 2^24, R = 64 multipliers, distinguished density 2^-5, S from the same function that
produces the production budget):

| run | groups | recovery | model |
| --- | ---: | ---: | ---: |
| corrected rule, f = 3/8 | 9600 | 0.3968 +- 0.0050 | 0.3923 |
| corrected rule, f = 0 control | 9600 | 0.3971 +- 0.0050 | 0.3925 |
| v5 rule (frozen per-trail twist) | 9600 | 0.0688 +- 0.0026 | 0.3923 |

The third row is the defect reproduced on identical fields, tables and seeds: 5.8x below
the model, the factor twelve the broken rule loses. The equivariance probe on 4000 samples
puts both images of a class in one next class on 4000 of 4000 samples under the corrected
rule, and on 1716 of 4000 under the v5 rule. Matched density is disclosed rather than assumed:
realised steps over S is 1.0090, and reports per group 21.9 against 21.9 expected, so the
measured success is not bought with extra work. The injectivity charge 1/(1-f^2) = 1.0787 is
also measured, as the trail-length ratio between the mixed walk and the f = 0 control:
1.0815 against 1.0787, accurate to well under a thousandth of a bit.

**End-to-end toy** (`toy_core.py`, 40 groups, tower x = -56: p = 10647363907, r = 9831361,
N = 819280 classes). The class rule, the equivariant step, the 16-multiplier table with twelve
twists, the step mixture, the detector, the ordered store, the report cap with refusal, the
replay with (a,b) exponent bookkeeping, the Frobenius alignment, the exponent solve and the
verification against z are all as in production; only the field size, N and the distinguished
density (2^-7) are toy. 40 of 40 groups recovered, each one verified against z and matching
the planted secret, with no alignment failure and no replay mismatch; a store capped at one
report refused and set the loss flag as designed. The cycle regression is green at 448/448
(lengths 2 to 8 at all 64 phases) and the adversarial set {2,4,8,16,32,64,67} fires at all 64
phases of every length, including the length above the 64-step window. The equivariance probe
on real tower elements puts all twelve images in one next class on 60 of 60 samples.

**Routing** (`routing.py`, 120000 keys): bit-disjoint fields, 16369 of 16384 table indices
hit, no distinguished word zero, chi-square 16459 on 16383 degrees of freedom, no collapse.

## Ledger

    rows per step                                  instructions   units
    dense Fp12 multiply (reference operation)             20000  1.0000
    cyclotomic squaring (3 steps in 8)                    11632  0.5816
    fold (every step, ten itemised rows)                   2906  0.1453

    S = least integer with S^2 (1-f^2) >= 2 mu N ln(1/(1-0.392)),  f = 3/8, mu = 11/10
    log2 S        = 125.81064
    c_mul = 1.14530   c_sqr = 0.72690   c = (5/8)c_mul + (3/8)c_sqr = 0.98840
    W_raw         = S c / C            = 2^125.79380
    W_stress      = every judgement-bounded fold row at 13/10 = 2^125.85607
    claim         = 125.87 : 0.07620 above W_raw, 0.01394 above W_stress, < live head 125.90

The fold is 2850 instructions of itemised work plus the 56-instruction witness update that
v5 broke out as its own row; this revision asserts the ten-row sum in the certificate
arithmetic rather than quoting it. The claim sits 0.0762 bits above the itemised work and
0.0139 above a ledger in which every judgement-bounded row is carried at 13/10, and the
window rule W_raw + 2e-5 <= claim < 125.90 is checked as an exact rational inequality.

## Scope, assumptions and what is not claimed

Success at 0.39 is conditional on three disclosed hypotheses: that the class walk behaves as
a random map to within the 10 per cent inefficiency already inside the step budget; that the
hash fields derived from the class key are pseudorandom; and that the algorithm's coins are
ideal. The extrapolation of walk quality from N <= 2^24 to 2^251 classes is the principal
uncertainty and is disclosed as such. Instruction counts are analytical bounds at the
organiser's pinned prices, not measurements of this code; the cyclotomic-squaring count is
the exact operation count of the implemented formula, and the fold rows are derived bounds
with margin. Peak memory is below 2^54 bytes. No production run was performed, no practical
break is claimed, no novelty is claimed, and no score is claimed to be accepted.

## Reproduction

    python3 part_a.py             # 20 exact checks on BLS12-381, real pairing, ~4 s
    python3 sim.py 1200           # walk simulation, eight fields, ~4 min
    python3 toy_core.py           # 40 end-to-end toy groups, N = 819280, ~4 s
    python3 routing.py            # 120000-key routing histogram
    python3 verify_certificate.py # re-derives the ledger and re-reads every run

Each script is standard-library only, writes a small JSON next to itself, and prints one
PASS/FAIL line. `verify_certificate.py` recomputes the itemisation sums, the step budget, both
ledgers in exact rational arithmetic, the claim-window inequalities including the 2e-5 rule,
the measured recovery as a Fraction against 0.39, the equivariance row, the matched-density
row, the bit-disjointness of the routing fields and every recorded check, and it fails the
package if any module restates a price, the mixture, the target or the claim instead of
importing it.
