# Deterministic sign-quotient grid with supplied-power root recovery

## 1. Claim and input domain

The input is T=A_1=[s]P0 and the literal sequence A_e=[s^e mod r]P0,
with s uniform in 1,...,r-1, on the protected BLS12-381 G1 subgroup.
For every integer 2^20 <= n <= 2^32 use d=2^20 and q=(r-1)/d.
The deterministic attack below has H(n)=0 and

    W(n) = [N*(79650816+21135360+65536+256)
            +21*79650816+2^80+2^20*n]/2000,
    a=98748750602204457565883422358598783,
    b=98748750602204457565883422358598784, N=a+b.

Success is exactly 2ab/q >= 39/100. At scoring n=2^20 the claimed
ceiling is W(n)<2^(3331/25), with 3331/25=133.24 reduced. Peak memory
is at most 1536L+2^64<2^129 bytes, L the least power of two >=N.
No probability model of a walk or hash table occurs. All failures,
preprocessing, addressing, code distribution and output are charged.

## 2. Public scalar setup and mixed radix

The scalar order r is the hexadecimal constant in verify.py. Its complete
factor list is supplied there; exact trial primality checks and the Lucas
criterion with witness g=7 prove r prime and g primitive. These are public
constants, not target-dependent advice; even trial-discovering the factors
fits the cold reserve as detailed in ledger.md. Set zeta=g^d mod r.

    q=50006747412801924209068050869165387952509453297164571592906626415194112.

The order q is EVEN and composite, with a residual factor 2^12. In
particular zeta^(q/2)=-1 mod r. Write Q=q/2. Every nonzero s has a UNIQUE
representation s=g^(q*k1+k2) mod r with 0<=k1<d and 0<=k2<q: this is
ordinary division of an exponent by q, not a Chinese-remainder claim.
No coprimality of q and d is needed. The actual supplied A_d equals
[zeta^k2]P0. Only known integers are modularly exponentiated.

## 3. Quotient keys; ORIGINAL geometric chains

Generate and store before updating:

    B_0=P0, B_(j+1)=[zeta]B_j,                 0<=j<a;
    G_0=A_d, G_(i+1)=[zeta^(-a) mod r]G_i,     0<=i<b.

Execute a updates and b updates, including the two unused final updates.
EVERY update is a complete known-scalar point action costing E=79650816
instructions. It is not a fixed-point addition and does not exponentiate
a hidden scalar. The same complete action price applies to every root
bit's candidate evaluation and to final verification.

For each stored COPY, set u=1 if its original y exceeds p-y, else u=0,
and key canon(x,y)=(x,min(y,p-y)) with a separate infinity flag. The
canonical representative is (-1)^u times the original point. Store the
orientation in the existing side word; keep the chain's live point
ORIGINAL for its next complete scalar action. Canonicalization is RECORD
MAINTENANCE ONLY, not a group operation or a quotient walk transition.
An extra 256 instructions is charged PER RECORD, without spending slack.

The 128-byte record is x48 || canonical_y48 || infinity8 || index16 ||
side_orientation8. Bit 0 of the last word is the side and bit 1 is u;
unused bits are zero. Stable 104-pass base-256 radix sorting compares all
104 key bytes, followed by an adjacent cross-side scan. Index and side
are payload, not key bytes. The infinity key is all-zero coordinates and
flag=1. Actual chain points are nonidentity because their scalar factors
are nonzero mod r. No nonidentity point has order two in this odd-order
subgroup, so the orientation is unique; no tie case loses a solution.

## 4. Exact signed exponent; distinctness and tiling

If a baby with index j and orientation u matches a giant with index i
and orientation v, equality of the canonical points implies

    [zeta^(j+u*Q)]P0 = [zeta^(k2-a*i+v*Q)]P0,
    k2 = a*i+j+(u-v)*Q mod q.

P0 has prime order r; equality first gives equality of scalar residues
mod r, and the exact order q of zeta then gives this exponent congruence.
The orientation fixes the unique lift from a quotient class to [0,q).
Dropping it gives k2 only mod Q. A wrong-sign candidate is NOT repaired
by adding r: [c+r]P0=[c]P0, not -[c]P0. The final scalar check is charged,
but is not used to pretend a missing orientation supplies full coverage.

Our rectangle satisfies 0<ab<Q. Babies have distinct canonical keys:
a duplicate requires j-j'=0 mod Q, impossible for 0<|j-j'|<a<=ab<Q.
Giants also have distinct keys: a duplicate requires a(i-i')=0 mod Q,
but 0<|a(i-i')|<=a(b-1)<Q. This argument needs NO assumption that a
is invertible modulo Q. Thus an equal-key run has at most two records,
one per side, and an adjacent scan cannot miss or mispair it.

For 0<=i<b the intervals [ai,a(i+1)) are adjacent and disjoint, with
union [0,ab). Division with remainder gives a unique (i,j) for every
integer in that union. The quotient match occurs exactly when k2 mod Q
is in [0,ab). Its two lifts [0,ab) and [Q,Q+ab) are disjoint subsets of
[0,q), each of cardinality ab. Hence the EXACT coverage identity is

    |covered exponents| = ab+ab = 2ab,
    Pr[grid succeeds] = 2ab/q >= 39/100.

This is an integer tiling theorem, not random occupancy or a birthday
estimate. The certificate checks ab<Q, both distinctness bounds, the
interval endpoints, the disjoint lifts, and 100*2ab>=39q. Tiny instances
also enumerate every tile and every k2. On a miss return failure after
the finite scheduled grid; there is no retry, rescue or free continuation.

## 5. Recover the complete scalar from literal supplied powers

After recovering EXACT k2, initialize h=0. For t=0,...,19 maintain
h=k1 mod 2^t and put

    e=d/2^(t+1), c=g^(e*(q*h+k2)) mod r, Cpoint=[c]P0.

Read literal A_e at array slot e-1. If A_e=Cpoint set bit=0; if
A_e=-Cpoint set bit=1; otherwise reject. Set h<-h+bit*2^t. Indeed,
writing k1=h+2^t*z gives e*q*2^t=(r-1)/2 and s^e=c*(-1)^z mod r.
All e are supplied; Cpoint is nonidentity, so exactly one sign applies.
The 20 integer identities are in the certificate. Return
s_hat=g^(q*h+k2) mod r only after 1<=s_hat<r and [s_hat]P0=T, serialized
canonically in 32 bytes. Twenty candidate actions plus one final action
cost 21E; comparisons, known modular work and encoding are charged too.

This is not multiplicative Pohlig-Hellman. Given [x]P0, a scalar action
[c]([x]P0)=[cx]P0 does NOT supply [x^c]P0. Smoothness of q does not
supply that missing operation. Here the large grid uses genuine charged
scalar actions, and each later unknown power is an actual input A_e.
No digit projection is invented; no canonicalized point becomes chain
state. This directly avoids both the encoded-power defect and an
uncharged geometric-update defect without a probabilistic replacement.

## 6. Scored certificate and evidence scope

ledger.md derives every analytic implementation ceiling under the
protected COST_MODEL.md: C=2000, H=0, R=19918011671547009368992591426495508435762944 at n=2^20.
All instruction rows, including complete actions and signs, are stressed
by 6/5. verify.py proves stressed work <2^(133.24-0.10) with exact integer
powers; 133.23 FAILS that same policy. A positive-term Fraction logarithm
enclosure gives the policy floor in [133.234231837952,133.234231837953].
The claim is between 0.005 and 0.06 bits above that floor. Both documents
and captured evidence are validated before the final certificate is made.

The real-EC toy executes the same original chains, exact 128-byte record
layout, 104 radix passes, signed reconstruction and literal-power lifting.
Eight seeded configurations have 2000 random trials each; 1716 exhaustive
k2 tests cover all configurations, and a further 12288 trials enumerate
all nonzero secrets at d=64. No wrong exponent or coverage error occurs.
The exact production success follows from the proof, not toy hit rates.
Production instruction kernels and astronomical memory are analytical
bounds, not a measured attack. dependencies.json=[]; the argument is local.
