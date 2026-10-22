# bls12-381-g1-ecdlp: Cheon's attack with a tame-wild class walk

Status: analytic cost bound for exploratory review. Heuristics are named in section 6.
No quantum step, restricted distribution or uncharged advice is used.

## 1. Setting and what is used

Challenge: s uniform in 1..r-1, input T=A_1=[s]P0 and A_i=[s^i]P0, i<=n, n=2^20. The attack reads only
B=A_d with d=2^20 (one point, 96 bytes) and T. Here r=z^4-z^2+1 with z=0xd201000000010000, and
2^32 divides r-1, so d=2^20 divides r-1. Put N=(r-1)/d (log2 N=234.857, N even and divisible by 3)
and let E=F_r^*{}^d be the subgroup of d-th powers; |E|=N. Cheon (Eurocrypt 2006; J. Cryptology 23, 2010)
reduces recovery of s to finding v=s^d in E from B=[v]P0, then extracting s from v (section 4).
Valid range: every n with 2^20<=n<r; W does not depend on n beyond needing A_{2^20}. For n<2^20 the bound is not claimed.

## 2. Phase 1: collision search for v

Let mu6 be the sixth roots of unity in F_r. They lie in E (-1 because N is even; the cube roots
because 3 divides N). On G1 the map phi(x,y)=(beta x,y), beta a cube root of unity in F_p, equals [lambda]
for a primitive cube root lambda in F_r (the matching pair is fixed by one scalar multiplication of P0).
The class of a point Q is {+-phi^i(Q)} (six distinct points). Let S={[e]P0 : e in E}. Both the tame points
[e]P0 and the wild points [e]B=[e v]P0 (e in E) lie in S, and S has M=N/6 classes (log2 M=232.272).

Walk. Fix 2048 jumps g_j=u_j^(2^20) for random u_j in F_r^*; each g_j is uniform in E. For a point Q define its
canonical representative Qc: among x, beta x, beta^2 x take the smallest Montgomery representative, then
negate y if its representative is odd. Then Qc=[eps lambda^i]Q. The step is Q -> [g_j]Qc, with j=bits 90..100 of
x(Qc). Since [g]([eps lambda^i]Q)=[eps lambda^i]([g]Q), the successor class depends only on the class, so
this is a walk on Z/M (logarithms in E modulo mu6 form a cyclic group, and the walk adds log g_j). The scalar of
the walked point is tracked exactly: e' = e * (eps lambda^i g_j), one F_r multiplication using a table of
6*2048 precomputed constants. The point [e']base is computed from scratch by fixed-base signed-window
multiplication (below), where base is P0 for tame and B for wild trails.

Trails start at e0=h^(2^20) (h random, 20 squarings, uniform in E), tame at [e0]P0 and wild at [e0]B. A
point is distinguished (DP) when the low 90 bits of x(Qc) are zero (theta=2^-90). At a DP the trail stores
(x(Qc), type, ec) with ec the canonical scalar. Same-type merges are discarded. When a tame DP (ec_T) and a
wild DP (ec_W) have the same x(Qc), the canonical points are identical, so [ec_T]P0=[ec_W]B=[ec_W v]P0 and
v = ec_T / ec_W mod r exactly (no residual ambiguity). Only tame-wild merges are useful.

Fixed-base multiplication. Signed digits of w=32 bits cover a 255-bit scalar in 8 windows. For each base and
window i, store the points j*2^(32i)*base for j=1..2^31 in affine form (96 bytes each): 2^35 points in all,
2^41.585 bytes (3.3 TB). A multiplication is 8 table fetches (negating y for negative digits) and 7 batched
incomplete affine additions. Exceptional cases (equal or opposite operands) have probability below 2^-250 per
addition, below 2^-127 over the whole run.

## 3. Success probability

Let m be the total number of steps, split equally between tame and wild trails. In the random-mapping model
the number of tame-wild merge events is Poisson with mean (m/2)^2/M, so P(success)=1-exp(-m^2/(4M)). For
P=0.39 the mean is ln(1/0.61)=0.49430, giving m0=sqrt(4*0.49430*M)=1.40612*sqrt(M)=2^116.628. The run uses
m=1.1*m0=2^116.765 steps to cover walk non-randomness (Teske 2001 measured r-adding walks with r>=16 within a few
percent of random; we take 10 percent in step count, and r=2048 is far above 16). With the margin the model
gives P=0.45, so the claim of at least 0.39 holds if the walk loses at most a factor 0.826 in collision rate.
Trails that have merged but not yet reached a DP are run on to their DP; this adds at most 2^10 trails * 2^90 steps
(mean) = 2^100 steps, a 2^-16.8 fraction. The DP-set is detected by sorting or hashing and no collision is missed.
Fruitless cycles of the Bos et al. 2010 / Bernstein-Lange-Schwabe 2011 kind come from the negation
map acting as a reflection; here mu6 acts by translations of the logarithm, so none arise. A cycle needs a
sum of at most L=20*2^90 random jumps to vanish modulo M, probability about L^2/M<2^-43 per trail, below 2^-16
over the 2^27 trails; trails are capped at 20/theta steps and restarted, which costs at most 2^94 steps each time.

## 4. Phase 2: from v to s

Compute s0 with s0^(2^20)=v by 20 successive square roots in F_r (Tonelli-Shanks, r-1=2^32*u). Since v is
a d-th power, v=g^l with 2^20 dividing l; both square roots have exponent divisible by 2^19 at the first stage (because
(r-1)/2 is divisible by 2^31), and so on, so every choice of root succeeds.
Cost at most 2^11 F_r multiplications per root. All solutions are s0*omega, omega in mu_d. Take omega0 of exact order d
(c^N for random c, test omega0^(d/2)=-1). Baby-step giant-step over i=2^10*b+a: store X_a=[s0 omega0^a]P0 (a<2^10, fixed-base
table of P0, 7 additions each) and test Y_b=[omega0^(-2^10 b)]T (b<2^10, each a generic scalar multiplication,
at most 2^10 units counting a doubling as 2 units). A match gives s with [s]P0=T, which is the attack-side check.
Total at most 2^20.1 units.

## 5. Cost ledger (1 unit = 2000 instructions = one batched incomplete affine addition)

Per step: 7 additions = 7 units. Other instructions: two F_p multiplications for beta x and beta^2 x
(2*383=766, using the measured blst fp_mul price for a 6-limb multiply; the 4-limb F_r multiply is cheaper but
is also charged 383), three-way 6-limb compare and y-sign normalisation 70, index and constant fetch 20, one F_r
Montgomery multiplication 383 (constants kept in Montgomery form, so no conversion), DP test and counter 10,
signed-digit recoding 48, eight table fetches (address, 6 load-pairs, conditional negation, store) 240, operand
handling around the seven additions 210. Sum 1747 instructions, charged as 2000 = 1 unit. A step costs 8 units.

| Item | Units |
| --- | ---: |
| Steps: 8 * 2^116.765 | 2^119.7653 |
| Trails running on to their DP | at most 2^103 |
| Table build: 2 bases * 8 windows * 2^31 additions, with doublings and exceptional cases, charged 2^36 | 2^36 |
| Jump constants (2048 * 20 squarings and 6*2048 products) | below 2^20 |
| Trail starts: 2^27 trails * (7 additions + 20 squarings) | below 2^31 |
| DP storage, sorting and matching: 2^26.8 DPs at 1000 instructions plus sort | below 2^28 |
| Phase 2, validation of B, square roots | below 2^21 |

W(2^20) <= 2^119.7653 * (1+2^-16) + 2^36 + ... < 2^119.7654. The submitted bound is log2 W = 119.8.
Table-build memory transfers are loads and stores counted above. The same table serves all trails, so a
single processor suffices; if P processors each hold a copy, add P*2^36 units (P<=2^40 changes nothing).

## 6. Heuristics and limits

H1 (random walk): the class walk collides like a random mapping on M classes, within the 10 percent step
margin (Teske 2001; van Oorschot-Wiener 1999 for distinguished points). H2: low bits of x(Qc) are uniform, so
DPs occur at rate 2^-90. H3: tame and wild trails are independent apart from merges. The instruction counts
are analytic upper bounds, not a built implementation; loads count one instruction and latency is not
charged, as the model fixes. Memory: tables 2^41.585 bytes, DP store 2^33.4 bytes, other state below 2^20;
total below 2^41.6 bytes. The memory is large but needs no new hardware class. Nothing larger than the
toy experiment was run. This is an exploratory analytic bound, not a practical attack.
