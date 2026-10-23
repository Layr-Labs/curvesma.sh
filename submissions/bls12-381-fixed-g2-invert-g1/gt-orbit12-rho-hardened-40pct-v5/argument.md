# GT orbit-12 square-and-add rho: hardened v5

## Claim and unchanged exact-Q mechanism
Track `bls12-381-fixed-g2-invert-g1`, C=20000. Claim 2583/20=129.15,
success >=0.39 under Hmix-v5/Hkey, peak <2^78 B, preprocessing <=2^50
included, advice zero. For the supplied nonidentity Q, g=e(P0,Q) generates GT;
solve z=g^x, return [x]P0, verify e([x]P0,Q)=z. Library: pairing(Q,P0,True). z=1 returns infinity. No powers or replacement Q.

For profile p,r, lambda=p mod r has order12, lambda^6=-1,
1+lambda^4+lambda^8=0; N=(r-1)/12. In Fp2[w]/(w^6-(1+u)), u^2=-1,
write c1=a+b*u, d=(1+u)^((p-1)/6), zeta=(1+u)^((p^2-1)/6) in Fp,
zeta^2=zeta-1, B=a*d0+b*d1. Interleave:

    even=(a,a*zeta,a*zeta-a,-a,-a*zeta,a-a*zeta)
    odd =(B,B*zeta,B*zeta-B,-B,-B*zeta,B-B*zeta).

Four Fp products form the keys. K=min(keys); kappa is its index, ties broken
by full conjugate then index. K is not injective: confirm equalities fully.
Sample independent a_h,b_h mod r, h=0..1023, J_h=g^a_h*z^b_h; precompute
J_h^(p^-k), lambda^-k*a_h, lambda^-k*b_h for k=0..11. For y=g^A*z^D,
h=K mod1024, k=kappa(y), the UNCHANGED v4 step is

    y'=y^2*J_h^(p^-k)
    A'=2*A+lambda^-k*a_h; D'=2*D+lambda^-k*b_h (mod r).

It induces [rep(y)^2*J_h], independent of history and chosen conjugate.
One shared table; two dense products. No ladder, restart or escape rule.

## Retained Floyd and coefficient analysis
Initialize tortoise/hare using independent uniform A0,D0 at g^A0*z^D0.
Per iteration advance tortoise once, hare twice. On full class equality,
k=lambda^(kappa(hare)-kappa(tortoise)), den=D_t-k*D_h, num=k*A_h-A_t;
x=num/den, verify g^x=z and the supplied-Q output. Zero denominator,
exception-cap exhaustion or budget exhaustion fails. Identity gives x=-A/D
when D!=0, or conservatively fails. For entry mu and cycle ell, first meeting
ell*max(1,ceil(mu/ell))<=mu+ell, including DP-free cycles.
Condition on initial/jump elements: D0 remains uniform with coefficient
2^t*(1-k*2^t) in den. Exact factor/order calculation gives
ord_r(2)/gcd(ord_r(2),12)>2T. Thus no coefficient closure for t<=T,
conditional degeneracy is1/r, and no uncharged replay is required.

## Single v4 defect: faithful-branch shared-table measurements
Rejected v4 a002fa35-4e79-4cb3-bdb6-7eba68fc543c accepted the mechanism and
ledger but found Hmix unsupported: dead toy branches and no shared-table
correlation measurements. Old 45/80,53/80 runs remain algebra evidence ONLY.
`scripts/faithful_audit.py` exhaustively enumerates each small GT subgroup
in the SAME tower, raw K/tie rule, h=K mod B, and literal DP bits. An EXACT
exponent oracle accelerates class(2*canonical_exponent+j_h); literal tower
cross-checks cover metadata and32 steps/setting. No hash replaces K.
2500 independent table/map pairs per setting; each table is kept throughout
its trajectory. Controls are fresh uniform finite maps with the same start.
Identity is a separate state: controls have N+1 uniform states; shared outputs
have uniform-exponent quotient weights. Every NONIDENTITY branch is live.

|q|r|N|B|branch min..max|mean shared/control repeat ratio|
|---:|---:|---:|---:|---:|---:|
|14009539|12289|1024|64|9..25|99641/100472|
|2677963|49177|4098|256|3..29|199004/200635|
|159508651|196657|16388|1024|6..30|404169/396092|
|19066387|786433|65536|1024|36..89|808771/791560|
|19066387|786433|65536|8|8146..8258|816627/808090|
|19066387|786433|65536|32|1947..2140|799939/807169|
|19066387|786433|65536|128|456..561|66490/66933|

Largest ratio~1.02175 is NOT a CDF or confidence bound. Let
t0=isqrt(5*ceil(N*B/(B-1))//4)+1, t1=ceil(5*t0/4).
Shared repeats by t1:1658,1601,1536,1561,1665,1603,1582 of2500;
control repeats by t0:1238,1157,1172,1172,1281,1194,1165.
Paired one-sided Hoeffding lower CDF differences (radius2/25) are
11/125,61/625,41/625,189/2500,46/625,209/2500,217/2500.
Endpoint radius1/25 has union failure<=28*exp(-8)<1/100 over seven settings;
all shared endpoint lower bounds exceed0.57. Confidence concerns TOY coins,
not production. B=8 stress reuses entries on796639/816627 transitions (>97%).
Full branch histograms; exact quantiles/sums; binned mu,ell,mu+ell,Floyd and
DP first-hit/censor distributions at bit widths2,4,6; regenerable full-data
digests are captured. 

## Widened budget and explicit production extrapolation
Np=ceil(N*1024/1023)=
4373927678342230290889091420689347427321206073683634826518910598626352813045;
the literal floor ends3044. T0=73941933961236014344098152017032175955,
4*(T0-1)^2<=5Np<4*T0^2. Widen ONLY the charged length:
T=ceil(5*T0/4)=92427417451545017930122690021290219944.
For a uniform map, P(mu,ell)=(n-1)_(mu+ell-1)/n^(mu+ell), and
P(mu+ell>T0)<=exp(-T0*(T0+1)/(2n)). At n=Np success>1-exp(-5/8)>0.4647.
**Hmix-v5:** the production shared-table raw-coordinate map satisfies
P(mu+ell<=T)>=the Np comparison probability at T0. This is an explicit25%
length-stretch extrapolation, replacing the unsupported5% assertion.
The preselected5/4 reserve passes the paired CDF tests, including heavy reuse;
mean ratios are separate diagnostics. log2(N) spans approximately10,12,14,16,
not production251.27. Scaling variables are log2(N),sqrt(N)/B,q/B, branch
occupancies and observed reuse. Small-B stress increases reuse, not every production scale. Finite confidence does NOT prove Hmix-v5; production correlations can
defeat it. The comparison-to-0.39 gap is NOT a transfer-error bound. Under the
explicit exploratory hypothesis, success exceeds the floor after negligible
identity, sampling, denominator and Hkey losses. Hkey assumes random-field-
scale exceptional keys; full equality/tie slow paths have a hard work cap.
Comparison identity loss<2^-120 and within-orbit ties<2^-240 are conditional.

## Complete inherited ledger, repriced for larger T
Per evaluation: two dense products plus AArch64 instructions:
4 Fp products1532; nine six-limb add/sub/neg576; min/index600;
two scalar double/add/reduce384; twisted loads/address256; DP/fp128;
Floyd comparison/control256; loop/meter128; state/scratch512: sum4372,
ceiling8000. Field carry/subtract/select uses six limbs; scalar updates use
four limbs, a carry limb, at most two conditional r-subtractions; twists and
labels are precomputed. Wexec=(36/5)T+2^100. Scored envelope:
W=8T+2^100=739419340880010743669210921667024964928;
129.1196<log2(W)<129.1197. Claim129.15 slack0.0303..0.0304 bits.
Stress ALL executed nonreference ceilings by6/5:
Wstress=(186/25)T+2^101, ceil=687659988374796133856571616751805647136
<2^129.015, reserve>0.135 bits. This reprices executed ceilings, not unused
8T envelope slack. Hard ancillary charges:129*M insertions at2^26
instructions; clearing129*M*128*16 instructions;2^32 candidates at2^15 units;
2^50 setup;2^90 exceptions/output/pairings. Sum<2^100, repriced<2^101.
All work/randomness charged;1024 rejection attempts/scalar, exhaustion fails.
No cross-target amortization.

## Auxiliary DP: unchanged hard caps, updated comparison
s=0..128: DP bits10..69+s zero; fp bits70+s..197+s,128 disjoint positions
(last325). Hfp is conditional uniformity, not an entropy theorem. M=2^69
records; clear/advance when full, disable after final generation or2^32
candidates. Floyd never changes/stops for DP events. Record128 B =
fp16+A32+D32+links32+metadata16. Peak128M+2^40<2^78 B includes all scratch.
Tree traffic, addressing, reconstruction and verification are charged. Finite false-pair expectation129*M*(M-1)/2^129
<2^17, Markov cap loss<2^-15 under Hfp. Larger T changes auxiliary Poisson
mean to<5M/12; FULL Chernoff<=(3/4)^M<2^-80 since exp(7/12)<9/5.
Reports are NOT exactly Poisson. DP caps/heuristics never debit core success.

`python3 -B scripts/verify_v5.py --write`:151/151 integer/Fraction checks,
not Hmix certification. `python3 -B scripts/faithful_audit.py --write`
replays17500 pairs; `python3 -B scripts/toy_audit.py` regenerates retained
v4 algebra results. Prior: mpjunior92,7da5f23a-6633-46b1-9c03-11f6cfedeeec,
AI-valid125.89; its hashed toy did not validate raw K. dependencies.json=[].
