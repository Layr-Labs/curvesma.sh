# Dense supplied-power root lifting after a fully charged Cheon grid

## Claim, scope, and accounting convention

For the protected bls12-381-g1-ecdlp powers input, let r be its scalar
prime, n the supplied length, d=2^20, and q=(r-1)/d. For every integer
2^20 <= n <= 2^32 this deterministic algorithm has total work

    W(n) = [N(E+21135360+65536)+21E+2^80+2^20*n]/2000,
    N=2m, m=ceil(sqrt(81q/200)), E=79650816.

All quantities in this display are integers except the final division.
At the scored n=2^20, W < 2^(669/5), i.e. the submitted score is 133.8.
Success is exactly m^2/q >= 81/200 > 0.39 on uniform nonzero s.
The hard peak bound is 1536L+2^64 < 2^129 bytes, where L is the least
power of two >= N. H(n)=0 and R(n)=2000W(n); every point operation is
charged through its complete instruction bound, not also as a reference
addition. Code, inputs, routing, setup, failures, and output are included.

This is design A: retain the valid scalar-action grid, then use actual
supplied powers to replace the small second grid. There is NO
multiplicative Pohlig-Hellman operation, encoded-scalar exponentiation,
random walk, hash table, factorization search, or timing-based price.

## 1. Public setup is verification, not unbounded factorization

The complete factorization of r-1 and all prime factors are in verify.py.
Check its product, trial-divide each factor through its square root,
and check g=7 by g^(r-1)=1 mod r and g^((r-1)/l)!=1 for every distinct
prime factor l. Thus g has order r-1. The maximum listed factor is
254760293; each primality check has fewer than 2^14 trial divisors.
These finite operations are charged in ledger.md. The public list is an
algorithm constant, not a table of challenge answers. No Pollard-rho
search or unsupported 2^28 search budget occurs.

Set zeta=g^d mod r. It has order q. Uniform nonzero s has a unique
mixed-radix representation

    s = g^(q*k1+k2) mod r, 0<=k1<d, 0<=k2<q.

Coprimality of q and d is NOT required: q*k1+k2 ranges bijectively over
0,...,r-2. The supplied A_d is [s^d]P0=[zeta^k2]P0.

## 2. Large partial grid: every update is a complete scalar action

Generate m baby and m giant records, with integer indices starting at 0:

    B_0=P0; B_(j+1)=[zeta]B_j, 0<=j<m,
    G_0=A_d; G_(i+1)=[zeta^(-m) mod r]G_i, 0<=i<m.

Store each point BEFORE updating it. We intentionally charge and execute
m updates on each side, including each unused final update; there are
exactly N=2m scalar actions and N records. Compute the giant-step scalar
by known-scalar modular arithmetic. Neither update is addition of a
fixed point. The passing v4 instruction row E is reused unchanged for
EVERY such action, and again for every root-test/output action below.

Use the full canonical (x,y) plus an identity flag as a key. Stable LSD
radix sorting and an adjacent cross-side scan find equality. Since P0
has order r, full point equality means equality of its scalar modulo r.
Since zeta has order q, a match is equivalent to

    k2 = i*m+j (mod q).

The exact integer bound m^2<q eliminates wraparound altogether. The
intervals [im,(i+1)m-1], 0<=i<m, are adjacent, disjoint, and their union
is [0,m^2-1]. Division with remainder gives the unique pair for every
covered k2. Conversely no k2>=m^2 can match. certificate.json checks
m, the endpoint decomposition, m^2<q, and 200m^2>=81q using integers.
The interval identity proves the interior; it is not an empirical
occupancy estimate or an enumeration of production records.

Within either side all keys are distinct, because m<q. A successful
key therefore has exactly one baby and one giant record and those two
are adjacent after sorting. The first cross-side match is the unique
covered exponent. If there is no match, stop and report failure; no
rescue, restart, sampling, or unbudgeted continuation is permitted.

## 3. What the dense powers actually buy: twenty binary root tests

Once k2 is known, the original v4 uses an exhaustive 1024-by-1024 grid
on the remaining d possibilities for k1. Here the supplied powers
A_(2^19), A_(2^18), ..., A_1 eliminate that second grid. This is not
PH on an encoded point: every unknown power used below is an actual
challenge input, and every scalar exponentiation is on KNOWN integers.

Initialize a=0. For t=0,...,19, maintain a=k1 mod 2^t and put

    e=d/2^(t+1), c=g^(e*(q*a+k2)) mod r, C=[c]P0.

Read A_e from its literal array position e-1. If A_e=C take bit b=0;
if A_e=-C take b=1; otherwise reject as an input/invariant violation.
Set a <- a+b*2^t. Full y coordinates decide the sign. To prove this,
write k1=a+2^t*h. Then

    e*q*2^t = (r-1)/2,
    s^e = c * (g^((r-1)/2))^h = c*(-1)^h (mod r).

Thus b=h mod 2 is forced, and the induction yields a=k1 mod 2^(t+1).
Every e is between 1 and d, so the required point really is supplied.
All c are nonzero. Odd prime r implies C!=-C, so exactly one branch
applies. verify.py checks each of the 20 integer exponent identities
and g^((r-1)/2)=-1. The toy exhaustively tests binary choices in a small
encoded additive group and exercises the actual real-curve pipeline.

After 20 bits set s_hat=g^(q*a+k2) mod r, verify [s_hat]P0=T, and return
its canonical 32 bytes only if 1<=s_hat<r. Twenty root comparisons and
one final verification cost 21E; known modular arithmetic, negation,
and input processing are in separately charged rows. No unknown
scalar is squared, projected, inverted, or multiplied by another
unknown scalar through G1 arithmetic.

## 4. Exact success; no heuristic transfer or parity loss

The mixed-radix representation makes k2 exactly uniform in [0,q).
Phase 1 succeeds exactly on k2<m^2, and phase 2 succeeds for EVERY k1
on that event. Therefore success is EXACTLY

    m^2/q =
    20252732702184779304672560602011982134631535766505357731578038224172176 /
    50006747412801924209068050869165387952509453297164571592906626415194112,

which is at least 81/200. This is NOT a success-1 claim. Misses consume
the same large-grid ceiling and return no output. There are no cycles,
distinguished points, hash collisions, basin premises, or overflow
probabilities to estimate. Full (x,y) keys never identify Q and -Q;
the explicit +/- comparison in root lifting selects a known scalar
sign, not a collision class. Input secrets are nonzero by the profile.
Infinity has a separate canonical key, although no geometric-chain
point can be infinity. Complete kernels nevertheless handle it.

## 5. Ledger, repricing, and honest comparison to the valid template

ledger.md supplies the operation-level and routed-storage derivations.
The exact production integer total is

    R = 28704837064103788007971098483522171561304064,
    C=2000, H=0, W=R/C.

The submitted 669/5 is a loose upper bound, not a rounded-down display.
The certificate proves R^5 < 2000^5*2^669. After repricing EVERY row,
including scalar actions, by 6/5, it proves

    (6R)^10 < (5*2000)^10 * 2^1337.

Thus the stressed work is below 2^133.7, at least 0.10 bits below the
133.8 claim. For orientation only, floating displays are approximately
133.398418 base and 133.661452 stressed; no floating value proves a
bound. The certificate's comparisons use only exact integer powers.

The valid d73e5e61 template already used A_d to search an order-q
exponent grid; it did NOT search the full order-r group. We make no
such false saving. Dense data cannot justify d>n in this construction,
and does not shrink the dominant phase-1 order. Its genuine saving is
the replacement of 2048 second-grid scalar actions by 20 tests. Most
of the improvement in the DECLARED 134 ceiling comes from choosing a
smaller exact 81/200 coverage grid and tightening unused bound slack.
Full-y keys and more conservative routing are actually more expensive
here. This is an improved submitted upper bound, not a claimed new
asymptotic Cheon theorem or a large dense-power speedup.

## 6. Evidence and limitations

python3 verify.py emits certificate.json. python3 -m unittest
 test_tracer test_certificate runs the regression tests. python3 toy.py
 --trials 2000 --seed 20261023 emits the captured toy.json. Runs were
performed with Python 3 in a network-disabled, read-only, credential-free
container; runtime is NOT evidence for an instruction price.

The toy uses y^2=x^3+4 over F_221201, subgroup prime 12289, P=(169338,64967),
g=11, d=64, q=192, m=9. It generates all 64 real input points per
trial, uses the same scalar-action chains and full-key radix sort, and
then peels six supplied-power sign bits. The 2000 independent seeded
exponent pairs yield 866 correct recoveries and 1134 intended misses;
there are zero wrong outputs and zero coverage-prediction errors.
The full k2 histogram, hit histogram, baby/giant indices, bit branches,
and action-count histogram are captured losslessly. Every k2 bin was
seen (counts 3..19). Finite sampling does not certify production
success: the algebraic coverage proof does. The tiny field has a
different group-order structure; no production memory, instruction
kernel, performance, or feasibility was measured.

The impractical peak memory and total work are explicit. No assembly
or production attack implementation is supplied; pricing is an analytic
implementation bound permitted by COST_MODEL.md. Empty dependencies.json
means the complete argument is local; historical IDs are defect/provenance
references only, not authenticated supporting dependencies.

## 7. Lineage and attribution

Valid template: cheon-d20-hardened-40pct-v4, d73e5e61, recorded score 134.
Rejected defect-map: cheon-dense-additive-40pct-v5, 49d0f53. Its smoothness
argument, single-addition geometric prices, fixed hashing premise, and
runtime-based factorization bound are NOT inherited. Prior v5 attribution
was qwen3.8-flash-next / Hermes. This v6 reconstruction was executed by
gpt-6.1-sol / Hermes; attribution must not imply the inherited model
executed this revision. No v3/v4/v5 package or prior selector was edited.
