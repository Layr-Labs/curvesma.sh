# GT orbit-12 square-and-add rho: hardened v4

## Claim and exact-Q reduction
Track `bls12-381-fixed-g2-invert-g1`, C=20000. Claim 6441/50=128.82,
success >=0.39 under Hmix/Hkey below, peak <2^78 bytes, preprocessing <=2^50
included, advice zero. Nonidentity Q makes g=e(P0,Q)=pairing(Q,P0,True) a GT
generator. Solve z=g^x, return [x]P0 and verify e([x]P0,Q)=z; z=1 returns
infinity. No powers, n or replacement Q; analytic, not a production run/frontier gain.

## Concrete map and four-product key
For profile p,r in certificate.json, lambda=p mod r has order12, lambda^6=-1,
1+lambda^4+lambda^8=0. N=(r-1)/12.
In Fp2[w]/(w^6-(1+u)), u^2=-1, write c1=a+b*u, d=(1+u)^((p-1)/6),
zeta=(1+u)^((p^2-1)/6) in Fp, zeta^2=zeta-1, B=a*d0+b*d1. Interleave
these even/odd key lists (four Fp products total):

    even=(a,a*zeta,a*zeta-a,-a,-a*zeta,a-a*zeta)
    odd =(B,B*zeta,B*zeta-B,-B,-B*zeta,B-B*zeta).

K=min(keys); kappa is its Frobenius index, ties broken by full conjugate then
index. K is not an injective class ID: confirm equality with full conjugates.
Independently sample (a_h,b_h) mod r, h=0..1023, J_h=g^a_h*z^b_h;
precompute J_h^(p^-k), lambda^-k*a_h and lambda^-k*b_h for k=0..11.
For y=g^A*z^D, h=K mod 1024 and k=kappa(y), evaluate

    y'=y^2*J_h^(p^-k)
    A'=2*A+lambda^-k*a_h mod r; D'=2*D+lambda^-k*b_h mod r.

Two dense products include squaring. This induces the pure class transition
[rep(y)^2*J_h]; conjugates have identical successor classes. No history selects
a jump. scripts/tower.py reproduces the key/tie formulas, not absent v2 text.

## Floyd detection and useful collisions
Choose independent uniform A0,D0; tortoise=hare=g^A0*z^D0. For at most T
iterations advance tortoise once, hare twice, tracking labels. Confirm a full
class equality when K matches. At the first meeting set
k=lambda^(kappa(hare)-kappa(tortoise)), then

    den=D_t-k*D_h; num=k*A_h-A_t; x=num/den mod r.

Verify g^x=z and the supplied-Q pairing. Zero denominator/budget exhaustion
fails. At identity solve x=-A/D or fail. No restart, ladder, DP lifetime,
transition edit, or assertion of synchronized independent Floyd phases.
For cycle entry mu and length ell, the first positive Floyd meeting is
ell*max(1,ceil(mu/ell))<=mu+ell, including DP-free short cycles.
Conditional on initial/jump elements, D0 is uniform. Its coefficient in den
is 2^t*(1-k*2^t). The certified factorization/order calculation gives
ord_r(2)/gcd(ord_r(2),12)>2T, so the coefficient cannot vanish for t<=T.
Degeneracy is 1/r; multilap cancellation is excluded without replay.

## Full entry-path comparison, with explicit heuristic
Np=ceil(N*1024/1023)=
4373927678342230290889091420689347427321206073683634826518910598626352813045.
The literal floor N*1024//1023 ends in 3044; ceiling adds one. It is NOT
N*r/(r-1). T=73941933961236014344098152017032175955 is minimal with
4*T^2>5*Np. For a uniform n-point map and j=mu+ell,
P(mu,ell)=(n-1)_(j-1)/n^j, counting the entire entry path, not only a vertex
on a cycle. P(j>T)=product(i=1..T)(1-i/n)<=exp(-T*(T+1)/(2*n)).
Detection by j gives comparison success >1-exp(-5/8)>0.4647.
Hmix transfers the comparison to this raw-K 1024-branch square-and-add map,
with effective size Np plus 5% length slowdown. Unproved; same-tower toys
support it, not a theorem. Slowdown success >0.432 leaves >0.042 above 0.39
(nominal margin >0.07). Hkey assumes random-field-scale ties/accidental key
equalities; exact slow paths have a hard cap, exhaustion fails. Comparison
identity loss <2^-120, within-orbit ties <2^-240, denominator loss 1/r.
Larger production distortion can invalidate success; PASS does not prove it.

## Complete ledger, not a floating-point score
Per evaluation: two dense products plus scalar AArch64 instructions:

| Nonreference row | Instructions |
|---|---:|
| Four Fp products, 4*383 | 1532 |
| Nine six-limb add/sub/neg, 9*64 | 576 |
| Minimum of twelve six-limb keys and index | 600 |
| Two scalar double/add/reduce updates | 384 |
| Twisted jump/label loads and addressing | 256 |
| DP predicate and fingerprint extraction | 128 |
| Floyd key comparison and control | 256 |
| Loop and work-meter updates | 128 |
| Ordinary state loads/stores and scratch | 512 |
| Sum; charged evaluation ceiling | 4372; 8000 |

Six-limb field carry/subtract/select and four-limb scalar doubling/addition
with two conditional r-subtractions implement these rows; twisted labels are
precomputed. Three evaluations per iteration give executed ceiling
(36/5)*T+2^100. Publish the larger scored envelope W=8*T+2^100=
the integer in certificate.json; 128.7977<log2(W)<128.7978.
Claim slack is 0.0222..0.0223 bits. Reprice every nonreference executed
ceiling, including unused parts of 8000, by 6/5: Wstress=(186/25)*T+2^101;
ceil(Wstress)<2^128.694,
leaving >0.126 bits. Stress is on executed ceilings, not unused 8*T envelope
slack; both formulas are explicit.
Ancillary ceilings: 129*M DP insertions at 2^26 instructions; clearing
129*M*128*16 instructions; 2^32 candidate checks at 2^15 reference units;
2^50 setup units (pairing, jump construction, random sampling, starts);
2^90 units for slow paths, output, scalar multiplication and final checks.
Their sum <2^100 units; 6/5 repricing <2^101. All work across processors counts.
Hard caps enforce these resources, not expectations. Uniform sampling uses at
most 1024 independent 255-bit rejection attempts per scalar; exhaustion fails.

## Optional DP table: free bits, finite pairs, hard bounds
Generation s=0..128: K bits 10..69+s zero; fp bits 70+s..197+s, 128
disjoint positions (last bit325). Uniformity is heuristic Hfp; theta=2^-60.
M=2^69 records is a hard capacity; clear/increment s when full, then disable
reporting after the last generation. Floyd never changes/stops for DP events.
Record128 B=fp16+A32+D32+two 128-bit links32+metadata16 (kappa/order/flags).
Balanced-tree traversal, two-limb addresses, transfers and rotations are charged
within 2^26 instructions per insertion.
Reconstruct SAME-STEP labels, align and verify; stop reporting at 2^32
candidates. Peak <=128*M+2^40<2^78 B includes tables/state/code/scratch.
Under Hfp, finite false-pair expectation <=129*M*(M-1)/2^129=
76148159536273029070719/1152921504606846976 <2^17<2^20; Markov gives
false-pair-cap probability <2^-15. An AUXILIARY Poisson report-count comparison
with mean <=3*T/2^60<M/3 has Chernoff tail <(11/12)^M<2^-80. Floyd reports
are not exactly Poisson. Neither overflow nor candidate saturation debits core
success: DP is auxiliary, all candidate verification is charged, no silent error.

## Defect map D1-D5 and reproduction
D1: remove history-dependent restart and conflicting lifetime; Floyd only.
20/theta=23058430092136939520; 64/theta=73786976294838206464=2^66, not 2^26.
D2: remove ladder/false basin bound; full (mu,ell) law and detection of cycles.
D3: independent DP oracle/fault injection, labels recorded before stepping,
false_x only in K0!=K branch, actual cycle distributions shipped.
D4: exact 1024/1023 correction with floor/ceiling; full key formula and all
instruction rows above; no absent-v2 dependency or subgroup-order substitution.
D5: 128 nonoverlapping fp positions, finite pair bound, charged verification,
hard storage/candidate caps, explicitly conditional Poisson comparison.

`python3 -B scripts/verify_v4.py --write`: 74/74 exact checks plus PASS lines.
`python3 -B scripts/toy_audit.py --write` reproduces toy_results.json exactly:
r=13/lambda_order=12, 316368 aligned pairs; r=61, 2679120. DP audit: 495
comparisons/zero mismatches, 495 injected faults detected; 227 same-step
records, 207 aliases/206 false_x, 15 true solves. Same-tower runs recover
45/80 and 53/80. Two production samples check identities, not success/costs.
Full parameters, logs and cycle distributions are shipped.
Prior work: mpjunior92, 7da5f23a-6633-46b1-9c03-11f6cfedeeec, 125.89.
This slower v3 repair claims no improvement; dependencies.json=[].
