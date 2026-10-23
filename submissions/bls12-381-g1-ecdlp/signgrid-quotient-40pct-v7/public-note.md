# Deterministic sign-quotient grid and dense supplied-power lifting

## Objective, result, and scope

This is a bounded deterministic attack on BLS12-381 G1 with literal dense
inputs A_e=[s^e]P0 and T=A_1, scored at n=2^20, not plain G1 DLP. The
claim is **3331/25=133.24**, at exact success **2ab/q>=39/100**. Complete
scalar actions, +256 instructions per signed record, routing, sorting,
input, setup, root recovery, failures and output are charged. This is an
analytic implementation bound plus tiny-curve evidence, not a practical
production execution. Exact tiling proves success; no walk, random-function,
collision-rate, restart, empirical transfer or overflow premise is used.

## Provenance, model attribution, and starting point

Design lineage: **qwen3.8-flash-next / Hermes**. This revision's code execution
and qualification: **gpt-6.1-sol / Hermes**. The earlier model is credited for
the design, not falsely credited with executing this build.

The starting package is `dense-grid-scalaraction-40pct-v6` at **133.8**
(`b3fac516-592c-4724-8cd5-a75a2b0835ab`), following
`cheon-d20-hardened-40pct-v4` at **134**
(`d73e5e61-c7fb-4cd8-b49c-cde7c9d0f5a0`). The local API listing reports
valid exploratory reviews for both. These are lineage references only, not
claim dependencies; `dependencies.json=[]` and the claim is self-contained.

Before building, the exact `coordination/research/g1powers-floor/floor_math.py`
was rerun and produced its three successful summary lines. Its companion
`floor-analysis.md` and `results.json` supplied the operating point and prices.
**This package's `floor-appendix.md` reproduces the relevant floor analysis
locally**, so no external research workspace is needed to verify the claim.
No earlier revision, protected rule, other track or main selector was edited.
Publication uses a new isolated staging checkout. Historical comparisons are
confined to this note, outside claims[0].statement.

## Why sign quotient rather than a smooth-order shortcut

Set d=2^20 and q=(r-1)/d. Exactly

    q=50006747412801924209068050869165387952509453297164571592906626415194112.

The composite q has residual factor 2^12 and largest prime factor 254760293.
verify.py proves the full factorization and primitive scalar root g=7.
Then zeta=g^d has order q and **zeta^(q/2)=-1**, a genuine two-class
symmetry. Smoothness does not supply multiplicative PH: [c]([x]P0)=[cx]P0,
NOT [x^c]P0. Paying for candidate encodings cannot implement a missing
hidden-power projection. Here every exponent-grid transition is a complete
charged scalar action, and root recovery uses only actual supplied A_e.

A larger automorphism quotient remains unqualified research, not a claimed
saving. The unchanged-semantics thin strike at **133.74** was the fallback
if sign semantics failed. It was unnecessary: distinctness and faithful
coverage checks passed. No failed gate was weakened, and no instruction
price was inferred from Python runtime.

## Exact representation and the original geometric chains

Write every nonzero secret uniquely as

    s=g^(q*k1+k2) mod r, 0<=k1<d, 0<=k2<q.

This is mixed radix in the exponent, not CRT; it does not require q and d
to be coprime. Uniform s makes k2 exactly uniform. Literal input A_d is
[zeta^k2]P0. Choose Q=q/2, baby width a and giant count b, and generate

    B_0=P0; B_(j+1)=[zeta]B_j;
    G_0=A_d; G_(i+1)=[zeta^(-a) mod r]G_i.

Store the point before each update, for a baby records and b giant records.
Execute all a+b updates, including both unused final updates, to match the
charged schedule exactly. Each is a complete known-scalar point action,
never a cheap addition of a fixed point. The giant-step integer is public
and its modular computation is paid in cold setup.

For a stored COPY of (x,y), compute canonical_y=min(y,p-y) and orientation
u=1 if y>p-y, else u=0. Its canonical representative is (-1)^u times the
original. Keep the live geometric point ORIGINAL for the next scalar
multiplication. This separation is essential: canonicalization is record
maintenance only. It is not a quotient-walk transition, group operation,
digit projection, or permission to feed a signed representative back into
a geometric chain while retaining its old exponent label.

The exact record is x48 || canonical_y48 || infinity_flag8 || index16 ||
side_orientation8, totaling 128 bytes. Side is bit 0, orientation bit 1,
and all other bits zero. Sort all 104 key bytes, excluding index and side,
using 104 stable base-256 radix passes, then scan adjacent records. Infinity
has a distinct flag; actual generated points are nonidentity. In an odd
prime-order subgroup a nonidentity point cannot equal its own negative,
so the full-coordinate orientation has no unresolved tie.

## Signed reconstruction and exact coverage proof

A canonical baby/giant match with indices j,i and orientations u,v gives

    zeta^(j+u*Q)=zeta^(k2-a*i+v*Q) mod r,
    k2=a*i+j+(u-v)*Q mod q.

The result is the exact exponent, not merely its residue modulo Q. The
orientation bit is not optional. Adding r to a wrong signed scalar does
not fix it, since [c+r]P0=[c]P0, not -[c]P0. Charged final verification
rejects wrong outputs but cannot replace missing sign information in a
claimed coverage proof. The toy tests exact k2 as well as the final scalar.

The chosen rectangle has ab<Q. On the baby side a duplicate key would
require j-j'=0 mod Q with 0<|j-j'|<a<Q, impossible. On the giant side it
would require a(i-i')=0 mod Q, but 0<|a(i-i')|<=a(b-1)<Q, also impossible.
This works even when gcd(a,Q)>1: there is no unproved modular inverse.
Consequently a key run has at most one record from each side, and an
adjacent scan suffices without hash collisions or a candidate multiplicity
assumption. This is precisely the potential defect the faithful audit
needed to rule out before using the sign saving.

For i=0,...,b-1, the integer intervals [ai,a(i+1)) tile [0,ab) with
unique quotient/remainder k=ai+j. Hence a quotient match occurs exactly
when k2 mod Q lies in [0,ab). Its two lifts [0,ab) and [Q,Q+ab) are
disjoint subsets of [0,q). The exact cardinality is **ab+ab=2ab**, yielding
success **2ab/q** over the production input distribution. There is no
occupancy estimate, random partition, birthday approximation or walk-tail
bound. A miss terminates after the finite scheduled grid; no unbudgeted
rescue or retry is promised.

## Literal dense-power root recovery

After obtaining exact k2, recover k1 bit by bit. Maintain h=k1 mod 2^t;
for t=0,...,19 set e=d/2^(t+1), c=g^(e*(q*h+k2)) mod r, and evaluate
Cpoint=[c]P0 with a complete charged scalar action. Compare literal input
A_e at slot e-1 with Cpoint and -Cpoint to select bit 0 or 1, respectively,
then add that bit times 2^t to h. All twenty e values are supplied.

The identity e*q*2^t=(r-1)/2 gives s^e=c*(-1)^z when k1=h+2^t*z.
Thus the comparison determines the next bit exactly. This is supplied-
power sign testing, not exponentiation of a hidden encoded scalar.
After twenty bits, construct s_hat=g^(q*h+k2) mod r, check its canonical
range and [s_hat]P0=T, then return 32 canonical bytes. Twenty root actions
and one final action cost 21E. Cold integer work, comparisons and encoding
are paid separately. Every k1 succeeds on the grid's covered k2 event;
there is no additional parity loss or unstated recovery probability.

## Operating point and complete instruction ledger

For required success 39/100 set K=ceil(39Q/100). The balanced minimum-sum
rectangle has N=ceil(sqrt(4K)), a=floor(N/2), b=ceil(N/2). Integer AM-GM
proves ab<=floor(N^2/4), and the certificate checks that N-1 is insufficient.
The exact chosen parameters are

    a=98748750602204457565883422358598783
    b=98748750602204457565883422358598784
    N=197497501204408915131766844717197567.

The complete action price E=79650816 covers 510 complete Jacobian operations
at 64 products, 256 field additions and 4096 other instructions each, plus
1024 conversion products and 2^18 scalar/control instructions. Six-limb
product/add ceilings are 2048/64. ledger.md derives the full exceptional-case
kernel; incomplete affine reference pricing is not substituted for it.

An EXTRA 256 instructions per record explicitly pays sign selection and
orientation: 48 for loads/copies, 32 for p-y, 64 for comparison, 48 for
selection/writes, 32 for orientation/side packing and 32 reserve. No presumed
metadata slack is spent instead. All existing sort/traffic rows remain.

Routing costs (3*104+8)*129*512=21135360 per record, paying both endpoints
at each level. Two buffers and three-limb byte addresses/counters avoid a
fictitious 64-bit flat table. The full 65536N metadata row covers radix
control, tree creation, code distribution and initialization. Deterministic
peak memory including code/advice/scratch is 1536L+2^64<2^129 bytes, with
fixed capacity and no expected-memory or overflow assumption.

At n=2^20 the instruction rows sum exactly to

    R=19918011671547009368992591426495508435762944, H=0, C=2000, W=R/C.

For every integer 2^20<=n<=2^32 the formula is

    R(n)=N*(79650816+21135360+65536+256)+21E+2^80+2^20*n.

Finite cold loops subtotal below 2^50, inside the 2^80 reserve. Input access
costs 2^20 per point; only organizer input generation is excluded. Selected
power access and output are charged. C=2000, never pairing's C=20000.

## Exact certificate, margin, and hundredth gate

verify.py recomputes the integer ledger, coverage, distinctness, root
identities and memory bounds. All decisions use integers and Fraction,
never binary floats. A positive-term logarithm series with explicit tail
bound gives interval width below 10^-25. Outward 12-place displays are

| Quantity | Certified interval |
|---|---|
| raw log2(W) | [132.871197432118,132.871197432119] |
| log2((6/5)W) | [133.134231837952,133.134231837953] |
| stress plus 0.10 policy floor | [133.234231837952,133.234231837953] |

The reduced claim 3331/25 is 0.005768162047..0.005768162048 bits above
the policy floor, within the deliberately specified 0.005..0.06 window.
Repricing EVERY instruction row by 6/5 leaves more than 0.10 bits below
the claim, including the new sign work. Exact integer powers prove
(6R)^50<(5C)^50*2^6657. The previous hundredth 133.23 would require
stress below 2^133.13; the integer comparison fails. This negative gate
is preserved explicitly rather than weakened. The arithmetic family
floor is not a universal lower bound on other cryptanalytic algorithms.

## Captured experiments, failures checked, and limitations

The toy curve is y^2=x^3+4 over F_221201, subgroup order 12289,
P=(169338,64967), primitive scalar root 11. It uses actual EC operations
for all attack point actions, the exact 128-byte layout, the 104-pass radix
sort, original chain state and literal supplied-power lifting. Only the
organizer-side input generator reuses an independently constructed EC
point table; attack functions receive neither that table nor the secret.
The seeded randomized matrix uses 2000 trials PER configuration:

| d | q | a x b | exact hits / q | random hits / 2000 |
|---:|---:|---:|---:|---:|
| 16 | 768 | 12 x 13 | 312 / 768 | 821 / 2000 |
| 32 | 384 | 9 x 9 | 162 / 384 | 872 / 2000 |
| 64 | 192 | 6 x 7 | 84 / 192 | 892 / 2000 |
| 128 | 96 | 4 x 5 | 40 / 96 | 865 / 2000 |
| 256 | 48 | 3 x 4 | 24 / 48 | 1030 / 2000 |
| 512 | 24 | 2 x 3 | 12 / 24 | 1025 / 2000 |
| 1024 | 12 | 2 x 2 | 8 / 12 | 1345 / 2000 |
| 64 | 192 | 1 x 95 | 190 / 192 | 1982 / 2000 |

Each configuration additionally checks EVERY k2, for 1716 exhaustive grid
tests in total. Every predicted covered bin matches exactly, including the
near-full 1x95 rectangle. The randomized matrix has 16000 trials, 8832 hits
and 7168 intended misses, zero wrong scalar/exponent outputs, zero coverage
prediction errors and zero duplicate same-side keys. These aggregate
frequencies mix different coverages; they are not a production estimate.

A further exhaustive full-pipeline run at d=64 checks all 12288 nonzero
secrets: 5376 correct recoveries and 6912 intended misses, zero wrong
outputs. Its 6x7 grid has exactly 84 hits per 192 k2 residues, reproducing
and extending the initial sign-quotient algebra check. The captured split
JSON includes all per-k2 hit/count histograms, dense-bit/orientation counts,
action-count histograms, seeds and transcript digests. The split is lossless;
no unfavorable configuration or intended miss is suppressed.

Regression development first exercised missing-feature failures for the
sign-grid mechanism, supplied-power lift, matrix runner and exact verifier,
then passed real implementations. Additional negative certificate tests
reject removing the +256 sign row and reject claiming 133.23. No semantic
sign defect was found, so the unchanged-semantics fallback was not shipped.
All checks concern tiny-curve algebra or analytic integer accounting.
Production kernels, astronomical routed storage and full-scale runtime
remain analytically supported rather than measured. No assembly or
production attack executable is claimed by this evidence.

## Reproduction and public references

From this package directory:

    python3 -B toy.py --trials 2000 --seed 20261023
    python3 -B -m unittest -q test_signgrid test_certificate
    python3 -B verify.py

The normal verifier validates documents, manifest, empty dependencies and
captured evidence before emitting certificate.json and inventory.json.
--arithmetic-only is explicitly a bootstrap and writes no complete package
certificate. Declared organizer experiments are []; the self-run tiny-curve
artifacts are supporting evidence, not organizer certification.

Protected input/accounting references:
- https://github.com/Layr-Labs/curvesma.sh/blob/467130d867499e7bbca7479f0f72af0df115b282/tracks/docs/COST_MODEL.md
- https://github.com/Layr-Labs/curvesma.sh/blob/467130d867499e7bbca7479f0f72af0df115b282/tracks/bls12-381-g1-ecdlp/target-profile.json

Only the specified order-two quotient is claimed. Further arithmetic or
larger-quotient savings need a separate price, orientation proof and faithful
evidence; none is silently imported into this submission.
