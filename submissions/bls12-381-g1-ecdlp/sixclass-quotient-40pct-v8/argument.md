# Six-class exponent quotient, with original geometric chains

## 1. Domain and complete claim

The protected inputs are T=A_1=[s]P0 and literal A_e=[s^e mod r]P0,
with uniform s in 1,...,r-1. For every integer 2^20<=n<=2^32 choose
fixed d=2^20, q=(r-1)/d and Q=q/6. The scalar order r and field p are
embedded in verify.py, together with the protected generator. Set

    a=b=57012617742321962895928165416718183, N=a+b=114025235484643925791856330833436366.
    H(n)=0,
    R(n)=N*(79650816+21135360+65536+8192)+21*79650816+2^80+2^20*n,
    W(n)=R(n)/2000.

The attack succeeds on exactly 6ab of the q equally likely exponents,
so success is exactly 6ab/q>=39/100. At n=2^20 the bound is
W<2^(2649/20), where 2649/20=132.45. Hard peak memory is at most
1536L+2^64<2^129 bytes, L the least power of two >=N. This is a finite
analytic implementation, not a measured production attack or a walk.

## 2. Which group has which order

verify.py proves the full factorization of r-1, factor primality and the
Lucas/primitive-root criterion with g=7. The composite q is

    50006747412801924209068050869165387952509453297164571592906626415194112.

In particular 6 divides q. Let zeta=g^d mod r, of exact order q in
F_r*. Every nonzero secret has a unique mixed-radix representation
s=g^(q*k1+k2), 0<=k1<d, 0<=k2<q. This is integer quotient/remainder,
not CRT; gcd(d,q)=1 is NOT needed. Then A_d=[zeta^k2]P0.

The EC subgroup has prime order r, NOT order q. Define
lambda=zeta^(q/3)=g^((r-1)/3) mod r and the nontrivial beta in F_p
embedded in verify.py. The verifier performs actual full-size EC scalar
multiplication to check [lambda]P0=(beta*x0,y0), and checks beta^3=1,
beta!=1, lambda^3=1, lambda!=1. Because y^2=x^3+4, tau(x,y)=(beta*x,y)
is an EC automorphism. Equality on P0 proves tau=[lambda] throughout G1.
Also zeta^(q/2)=-1. These facts establish scalar/field root alignment,
not merely the existence of two unrelated cube roots.

## 3. Exact orbit digit and exception handling

Enumerate each point's orbit in this order:

    digit u:    0      1       2       3       4        5
    m_u:        1    lambda  lambda^2  -1    -lambda  -lambda^2
    point:      P    tau(P)  tau^2(P)  -P    -tau(P)  -tau^2(P)
    n_u:        0      2       4       3       5        1
    m_u=zeta^(n_u*Q) mod r.

The six m_u are distinct. Every nonidentity subgroup point has six
orbit members: [m_u]P=[m_v]P would imply (m_u-m_v)P=O, impossible in
prime order r. Thus nonidentity fixed points of any nontrivial signed
tau power number ZERO. Potential full-curve x=0 or y=0 exceptions have
order dividing 3 or 2, respectively, and cannot occur here. O is given
a separate infinity key and digit zero; no actual chain point is O.

For a stored copy choose canon(P)=the lexicographically least (x,y)
among the six candidates, and retain its unique digit u. Compute only
two extra Fp products beta*x and beta*(beta*x), one y negation and
five fixed-width lex comparisons/selections. Canonicalization is RECORD
MAINTENANCE ONLY. Do not replace any live geometric chain point with
its canonical representative. Each subsequent grid update remains a
complete charged known-scalar action on the ORIGINAL point.

A record is x48||y48||infinity8||index16||side_digit8: 128 bytes. The
last word is side|(u<<1); three digit bits suffice. Sort the full 104
key bytes with 104 stable base-256 passes, then scan adjacent records.
There are no truncated keys, fingerprints or hash-table assumptions.

## 4. Grid and exact 6x6 reconstruction

Store before updating:

    B_0=P0, B_(j+1)=[zeta]B_j,             0<=j<a;
    G_0=A_d, G_(i+1)=[zeta^(-a) mod r]G_i, 0<=i<b.

Execute all a+b updates, including the two unused final updates. Every
update costs the complete action ceiling E=79650816; no fixed-point-add
substitution or encoded-scalar exponentiation is performed. For a
matching baby digit u and giant digit v, canonical equality means

    m_u*zeta^j = m_v*zeta^(k2-ai) mod r,
    G_i=[c(u,v)]B_j, c(u,v)=m_u/m_v in F_r*,
    k2=ai+j+Delta[u][v]*Q mod q,
    Delta[u][v]=(n_u-n_v) mod 6.

The correction c lies in the six roots of unity of F_r*, NOT in a
fictional q-order EC subgroup or an assumed field modulo composite q.
The correction to the exponent is ADDITIVE. The full table is:

| baby u \ giant v | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 4 | 2 | 3 | 1 | 5 |
| 1 | 2 | 0 | 4 | 5 | 3 | 1 |
| 2 | 4 | 2 | 0 | 1 | 5 | 3 |
| 3 | 3 | 1 | 5 | 0 | 4 | 2 |
| 4 | 5 | 3 | 1 | 2 | 0 | 4 |
| 5 | 1 | 5 | 3 | 4 | 2 | 0 |

For entry Delta=t, c=zeta^(tQ). In t=0,...,5 order these corrections
are 1,-lambda^2,lambda,-1,lambda^2,-lambda. Equality of points first
gives equality of scalar residues mod r, then exact order q of zeta
gives the exponent congruence. The stored pair selects the exact lift
k2 in [0,q), not just k2 mod Q. Adding r cannot repair a wrong digit.

## 5. Exact tiling, distinctness and success

The rectangle satisfies 0<ab<Q. Baby same-side duplicates require
j-j'=0 mod Q, impossible for 0<abs(j-j')<a<Q. Giant duplicates require
a(i-i')=0 mod Q, impossible for 0<abs(a(i-i'))<=a(b-1)<Q. No modular
inverse of a is assumed. An equal-key run therefore contains at most
one record of each side, and the adjacent scan finds its unique pair.

The adjacent intervals [ai,a(i+1)), 0<=i<b, tile [0,ab) uniquely as
k=ai+j. Matching occurs exactly for k2 mod Q in [0,ab). Its SIX lifts
[tQ,tQ+ab), t=0,...,5, are disjoint subsets of [0,q), each of length ab.
Thus the integer identity |covered|=sum_t(ab)=6ab proves exact success
6ab/q. The certificate checks the endpoints and 100*6ab>=39q. This is
counting, not occupancy, mixing, a birthday approximation or a walk.
A miss returns failure after its finite schedule, without free retries.

## 6. Recover all root bits from literal supplied powers

After exact k2, maintain h=k1 mod 2^t, initially h=0. For t=0,...,19,
let e=d/2^(t+1), c=g^(e*(q*h+k2)) mod r and evaluate Cpoint=[c]P0
with one complete action. Compare actual supplied A_e at slot e-1
with Cpoint and -Cpoint to determine bit 0 or 1, then h+=bit*2^t.
Indeed k1=h+2^t*z implies e*q*2^t=(r-1)/2 and s^e=c*(-1)^z mod r.
All twenty e values are supplied, and Cpoint is nonidentity. Return
s_hat=g^(q*h+k2) only after 1<=s_hat<r and [s_hat]P0=T, in 32 bytes.
The twenty candidates and final verification cost 21E. All remaining
known-integer arithmetic, input access, comparisons and output are paid.

A scalar action [c]([x]P0)=[cx]P0 does NOT compute [x^c]P0. No PH
projection is invented: the only unknown powers used are actual A_e.

## 7. Accounting, certificate and evidence scope

ledger.md derives the complete action, sort, routing, storage, cold and
input ceilings under C=2000. The EXTRA 8192 per record includes two
2048-instruction fixed-constant products, 64 for p-y, 1280 comparisons,
640 selections, 512 loads/scratch, 128 digit packing and 1472 reserve.
No metadata slack or reference-calibration interpolation is substituted.
At scoring R=11500574304558579830758679329268362223869952. Every instruction row
repriced by 6/5 still has log2 in [132.341864108206,132.341864108207],
leaving >0.10 bits below the claim. The policy floor after that reserve
is [132.441864108206,132.441864108207]. Claim slack is >0.005 and <0.06
above that floor; 132.44 FAILS the identical headroom test. All numeric
decisions use integers/Fractions with bounded-tail logarithm enclosures.

The real-EC fixture has p=257869, r=12289, #E=258069=21r and explicit
aligned beta/lambda. Exhaustive testing of all 12288 nonzero subgroup
points finds zero orbit-size or fixed-point exceptions. All 73728
within-orbit digit pairs check the complete table. Eight configurations
have 2000 randomized grid trials each, 1716 exhaustive k2 tests in all,
and an additional 12288 full-secret pipeline tests. Wrong exponents,
wrong scalars and coverage-prediction errors are all zero. The tests use
actual EC actions, original chains, full records and radix sorting.
These are tiny-curve algebra checks, NOT production instruction or
runtime measurements. Production success is proved by exact tiling.
Dependencies are empty; the evidence and verifier source are local.
