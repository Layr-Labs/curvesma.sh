# Exact finite public-parameter certificate

Read r as the integer with protected hexadecimal expansion
73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001.
Define all divisions and ceilings exactly as follows:

R=isqrt(r-1)+1
L=(1059*R+999)//1000
s=(L+61)//62
K=L+s=246410426563779385479884453130106700081
T=271; eta=1/65536; d=181; w=33554432
D=(L-T+1)*(L-T+2)/2
m=(1-eta)*D/r
U=(1+3eta)*D^2/r^2+(1+eta)*D/r
  +((724+722eta)/32400)*(L+1)^2/r
  +8*T*eta*(L+1)^2/r^2+12*eta*(L+1)/r

All quantities except R/L/s/K/T/d/w are exact rationals, not floating point.
The following are the finite inequalities used in the complete proof:

1. w >= 4*d^2*256; 4*r*3^256*2^100 < 8^256.
2. rho=1/d <=1/2; rho^T <= eta/r.
3. (3*m-U)/2 > 39015/100000.
4. (3*m-U)/2 -1/8192 -2^-100 -2^-980 -1/r >39/100.
5. K+1<2^129 and K<r.
6. (31008*K+2^90)^125 < 20000^125 * 2^16021.
7. 640*33554432+2^20=21475885056 <2^35.

Each is a finite rational/integer comparison; multiply by the positive
explicit denominators for a pure integer certificate. These are supporting
analytic checks, not executable experiment declarations or measured results.
Diagnostic floating approximations, not used as proofs: the net probability
lower bound is 0.3908488484572603, and the ledger logarithm is
128.16697138483252. The conservative first claim remains exactly 16021/125,
not either diagnostic. The bounds do not prove the fixed-SHA heuristic.
