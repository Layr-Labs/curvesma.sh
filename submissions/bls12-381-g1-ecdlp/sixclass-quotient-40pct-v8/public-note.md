# Six-class exponent quotient with exact orbit digits

## Objective, outcome and scope

This deterministic auxiliary-powers G1 attack uses literal A_e=[s^e]P0
and T=A_1, scored at n=2^20; it is not plain G1 DLP. The reduced claim
is **2649/20=132.45**, at exact success **6ab/q>=39/100**. Complete point
actions, canonical records, sorting/routing, input, setup, root recovery,
failures and output are charged. This is an analytic implementation with
faithful toy evidence, not an executed production attack. Success is exact
integer tiling, with no walk or probabilistic collision premise.

## Provenance and the prior accepted mechanism

Design lineage: **qwen3.8-flash-next / Hermes** (floor analysis and brief).
Execution and qualification: **gpt-6.1-sol / Hermes**, the actual runtime
identifier, presented as GPT 6.1 Sol for CLI attribution.

The starting package is `signgrid-quotient-40pct-v7`, submission
`037d0cbc-67e0-4380-95db-71799c304d6e`. The retrieved API record reports
valid exploratory review at 133.24, explicitly supporting original chains,
stored-copy canonicalization, exact signed recovery and deterministic
coverage. Earlier lineage is `dense-grid-scalaraction-40pct-v6` and
`cheon-d20-hardened-40pct-v4` (cheon-d20). This history is context only,
not an authenticated dependency or part of claims[0].statement.

The existing floor_math.py was rerun: its six-class 8192-extra-instruction
model predicts [132.441864108206,132.441864108207]. This package reproduces
it independently; floor-appendix.md embeds the relevant arithmetic.
Dependencies are empty. No prior revision, protected rule, other track or
main selector is edited. Publication uses a new isolated checkout.

## The crucial order and root-alignment distinction

Set d=2^20, q=(r-1)/d and Q=q/6. The full factorization, factor primality
and Lucas/order tests prove r prime and g=7 primitive. Exactly

    q=50006747412801924209068050869165387952509453297164571592906626415194112.

Six divides q. Let zeta=g^d, of order q in F_r*. The EC subgroup has
prime order r, NOT q; q is the multiplicative hidden-scalar order. Every
nonzero s has a unique mixed-radix form g^(q*k1+k2), 0<=k1<d, 0<=k2<q.
No coprimality of d and q is needed. Uniform s gives uniform k2, and the
actual supplied A_d is [zeta^k2]P0.

On y^2=x^3+4, tau(x,y)=(beta*x,y) is an automorphism for beta^3=1,
beta!=1 in F_p. Set lambda=zeta^(q/3). Two unrelated cube roots are not
an alignment proof: verify.py performs full-size EC multiplication to
check [lambda]P0=(beta*x0,y0), proving tau=[lambda] throughout G1.
Both constants and the protected generator are embedded. The public
alignment action fits charged cold setup. Also zeta^(q/2)=-1.

The reconstruction correction lies in F_r*, NOT modulo composite q or
in a fictional q-order EC subgroup. Its known zeta exponent supplies an
ADDITIVE correction to k2. Exhaustive real-EC checks verify this direction.

## Exact record semantics and chain schedule

Store each point before updating these ORIGINAL geometric chains:

    B_0=P0, B_(j+1)=[zeta]B_j,             0<=j<a;
    G_0=A_d, G_(i+1)=[zeta^(-a) mod r]G_i, 0<=i<b.

Execute a+b complete actions, including both unused final updates. Each
costs E=79650816, never a fixed-point addition. The known giant scalar is
computed in cold setup. A scalar action [c]([x]P0)=[cx]P0 is NOT [x^c]P0;
no PH projection or encoded-scalar exponentiation is invented.

For a stored copy enumerate:

    digit:       0       1        2        3        4         5
    point:       P     tau(P)   tau^2(P)   -P     -tau(P)   -tau^2(P)
    multiplier:  1     lambda   lambda^2   -1     -lambda   -lambda^2
    n_digit:     0       2        4        3        5         1.

Here m_u=zeta^(n_u*Q). Choose the least full (x,y) lexicographically and
retain its digit. Two Fp products produce the three x values; reuse y
and p-y. Canonicalization is **RECORD MAINTENANCE ONLY**: keep each live
point original for its next complete action, without quotient-walk
transitions, label drift or a cheaper geometric update.

The record is x48||y48||infinity8||index16||side_digit8: 128 bytes.
Side is bit 0, digit bits 1..3, remaining bits zero. All 104 key bytes
are sorted with 104 stable base-256 LSD radix passes, then an adjacent
cross-side scan. There are no fingerprints, hash lookups or x-only matches.

Distinct multipliers cannot fix a nonidentity point in prime order r.
Thus every nonidentity orbit has size six; full-curve x=0 or y=0 torsion
is outside G1. Infinity has a distinct flag and digit zero, but every
actual chain point is nonidentity. No tie or exception is silently lost.

## Complete digit-pair table and exact exponent recovery

A canonical baby (j,u) and giant (i,v) match proves

    m_u*zeta^j = m_v*zeta^(k2-ai) mod r,
    c(u,v)=m_u/m_v=zeta^(Delta[u][v]*Q) mod r,
    k2=ai+j+Delta[u][v]*Q mod q.

Rows are baby u, columns giant v; Delta=(n_u-n_v) mod 6:

| baby u \ giant v | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 4 | 2 | 3 | 1 | 5 |
| 1 | 2 | 0 | 4 | 5 | 3 | 1 |
| 2 | 4 | 2 | 0 | 1 | 5 | 3 |
| 3 | 3 | 1 | 5 | 0 | 4 | 2 |
| 4 | 5 | 3 | 1 | 2 | 0 | 4 |
| 5 | 1 | 5 | 3 | 4 | 2 | 0 |

For Delta=0,...,5 the correction c is
1,-lambda^2,lambda,-1,lambda^2,-lambda. Giant equals [c] times baby;
the exact k2 in [0,q) is recovered, not merely k2 mod Q. Adding r cannot
repair a lost digit. Charged final verification does not replace this
reconstruction proof. The toy checks every ordered digit pair in every
nonzero orbit, including field/scalar alignment and all canonical keys;
73728 relations pass with none wrong. Grid checks also compare exact k2,
not only the final recovered scalar.

## Exact six-lift tiling and success

The chosen ab is strictly below Q. A baby same-side duplicate would require
j-j'=0 mod Q with 0<abs(j-j')<a<Q. A giant same-side duplicate would require
a(i-i')=0 mod Q with 0<abs(a(i-i'))<=a(b-1)<Q. Both are impossible, even
when a and Q are not coprime. Therefore each equal-key run has at most one
record per side, and the adjacent scan cannot miss or mispair it.

The intervals [ai,a(i+1)) for i=0,...,b-1 uniquely tile [0,ab) as ai+j.
A match occurs exactly when k2 mod Q is in this interval. Its six lifts
[tQ,tQ+ab), t=0,...,5, are disjoint subsets of [0,q), each with exactly ab
integers. Hence |covered|=sum_t ab=6ab and success=6ab/q>=39/100. This is
an exact integer identity, not a random occupancy, mixing or birthday-rate
claim. A miss terminates after the bounded schedule, with no free retry.
Every k1 succeeds on the same covered k2 event once dense lifting follows.

## Literal supplied-power recovery of the remaining root

Maintain h=k1 mod 2^t, initially zero. At t=0,...,19 let e=d/2^(t+1),
c=g^(e*(q*h+k2)) mod r and Cpoint=[c]P0. Compare literal A_e at array
slot e-1 with Cpoint and -Cpoint to recover bit 0 or 1; update h by that
bit times 2^t. Writing k1=h+2^t*z gives e*q*2^t=(r-1)/2, hence
s^e=c*(-1)^z mod r. All twenty e values are supplied. Since Cpoint is
nonidentity, exactly one comparison applies. This is actual input-power
sign testing, not exponentiating a hidden scalar encoding.

After twenty bits, form s_hat=g^(q*h+k2), check 1<=s_hat<r and
[s_hat]P0=T, then output its canonical 32 bytes. Twenty candidate actions
plus one final action cost 21E. Modular arithmetic, comparisons, selected
input decoding, output and every attack-side input access are charged.
Only organizer generation of the supplied sequence is excluded.

## Complete ledger and the reproduced floor

For K=ceil(39Q/100), integer AM-GM gives ab<=floor(N^2/4), so choose
N=ceil(sqrt(4K)), a=floor(N/2), b=ceil(N/2). The verifier checks N-1 fails.

    a=b=57012617742321962895928165416718183
    N=114025235484643925791856330833436366.

The complete action ceiling pays 510 complete Jacobian operations at
64 products/squares, 256 field additions and 4096 other instructions
each, plus 1024 normalization products and 2^18 scalar/control work.
Six-limb product/add ceilings are 2048/64. Exceptional-case candidates
are selected, not inverted; dependent chains do not get affine batching.
The reference calibration fixes C=2000, not these implementation prices.

The EXTRA 8192 record row pays 4096 for two conservative Fp products,
64 for p-y, 1280 for comparisons, 640 selections, 512 local traffic,
128 digit packing and 1472 reserve. The subtotal is 6720. A public fixed
beta*R_M Montgomery constant returns ordinary beta*x from ordinary x,
so no conversion is hidden. This row is ON TOP of the whole metadata
row, not presumed slack or a Python-derived 256-instruction product.

Routing costs (3*104+8)*129*512=21135360 per record, paying both
endpoints and all addressing. Two fixed buffers and binary routing trees
use THREE-limb byte addresses, bucket counts, prefixes and scatter cursors.
The 65536N metadata row includes radix control, allocation, initialization
and code distribution. Peak memory is 1536L+2^64<2^129 bytes, with
L=166153499473114484112975882535043072; this is finite but astronomical, not
a 64-bit flat-address oracle or expected-occupancy claim.

For every integer 2^20<=n<=2^32,

    H(n)=0,
    R(n)=N*(79650816+21135360+65536+8192)+21E+2^80+2^20*n,
    W(n)=R(n)/2000.

Cold finite subloops total <2^50 inside 2^80. Input scanning pays 2^20
per point; selected A_e are retained. At scoring,
R=11500574304558579830758679329268362223869952. ledger.md itemizes every instruction row;
verify.py recomputes the exact sum without double-counting H and R.

## Numerical qualification and the failed adjacent gate

Exact Fraction logarithm enclosures, outward-rounded:

| Quantity | Interval |
|---|---|
| log2 W | [132.078829702372,132.078829702373] |
| log2((6/5)W) | [132.341864108206,132.341864108207] |
| stress +0.10 policy floor | [132.441864108206,132.441864108207] |

Repricing EVERY instruction row by 6/5 leaves 0.108135891793..794
bits below the claim, exceeding 0.10. Slack above the policy floor is
0.008135891793..794, inside 0.005..0.06. The integer inequality
(6R)^20<(5C)^20*2^2647 proves the stress gate; 132.44 FAILS the same
headroom requirement. Negative tests preserve that failure and reject
removing the orbit row. This is a local conservative policy, not an
organizer rule or a generic lower bound. All decisions use integers and
Fractions, not binary floats or empirical hit rates.

## Captured experiments, course corrections and limitations

The old sign fixture p=221201 lacks a nontrivial cube root in F_p, so it
was not blindly reused. A small-prime search found p=257869,
#E=258069=21*12289, P=(225183,257160), primitive scalar root 11,
beta=234548 and lambda=6240. The verifier independently counts the curve,
checks primality, point order and alignment, and audits all nonzero subgroup
points. Attack operations use actual EC arithmetic, original live chains,
full records and radix sorting. Only organizer input generation reuses a
repeated-addition point table; no attack function sees it or the secret.

The seeded matrix uses 2000 randomized trials per configuration, seed
20261023 plus index, and enumerates every k2 in every configuration:

| d | q | a x b | exact hits / q | random hits / 2000 |
|---:|---:|---:|---:|---:|
| 16 | 768 | 7 x 8 | 336 / 768 | 880 / 2000 |
| 32 | 384 | 5 x 5 | 150 / 384 | 781 / 2000 |
| 64 | 192 | 4 x 4 | 96 / 192 | 1033 / 2000 |
| 128 | 96 | 3 x 3 | 54 / 96 | 1145 / 2000 |
| 256 | 48 | 2 x 2 | 24 / 48 | 1037 / 2000 |
| 512 | 24 | 1 x 2 | 12 / 24 | 994 / 2000 |
| 1024 | 12 | 1 x 1 | 6 / 12 | 966 / 2000 |
| 64 | 192 | 1 x 31 | 186 / 192 | 1936 / 2000 |

There are 16000 randomized trials: 8772 correct recoveries and 7228
intended misses, plus 1716 exhaustive k2 tests. An additional d=64,
a=b=4 pipeline run checks all 12288 nonzero secrets: 6144 correct
recoveries, 6144 intended misses. Wrong exponents, wrong scalars,
coverage-prediction errors and duplicate same-side keys are all zero.
Mixed-coverage aggregate hit frequencies are NOT a production estimate.
Exact 6ab/q is independently checked per configuration.

All 12288 nonzero subgroup points have six-member orbits, with zero
fixed-point exceptions. All 36 ordered digit pairs in each of 2048
orbits pass: 73728 relations, none wrong. The lossless split capture is
15283 bytes including summary, retaining per-k2 count/hit histograms,
dense bits, digit pairs, action counts, seeds and transcript digests.
Intended misses and the near-full unbalanced stress are not suppressed.

Test-first development exercised missing mechanism, grid, lift, matrix
and certificate failures before real implementations passed. Negative
pricing tests reject a missing orbit row and the preceding hundredth.
The sign-only fallback would not improve the named frontier and was not
submitted; the six-class mechanism qualified without it.

These tests establish tiny-curve algebra, not production instruction
measurements, practical memory or full-scale recovery. Production success
comes from exact tiling, while kernels and routed storage remain analytical
bounds. Exploratory acceptance is plausible-not-refuted, not formal proof
certification, a novelty guarantee or a demonstrated production break.

## Reproduction and public references

From this package directory:

    python3 -B -m unittest test_sixclass test_certificate -q
    python3 -B toy.py --trials 2000 --seed 20261023
    python3 -B verify.py

The verifier validates documentation, captured evidence and exact prices
before writing certificate.json and recursive inventory.json. Arithmetic-
only bootstrap explicitly leaves documentation pending and writes no final
certificate. Source is embedded; experiments=[]. Protected intake limits
are unchanged. Historical lineage is context, not a claim dependency.

https://github.com/Layr-Labs/curvesma.sh/blob/master/tracks/docs/COST_MODEL.md
https://github.com/Layr-Labs/curvesma.sh/blob/master/tracks/bls12-381-g1-ecdlp/target-profile.json
