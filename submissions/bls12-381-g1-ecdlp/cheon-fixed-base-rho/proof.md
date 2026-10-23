# bls12-381-g1-ecdlp: Cheon's attack with a tame-wild class walk

Status: analytic cost bound for exploratory review. Heuristics are named in section 6.
No quantum step, restricted distribution or uncharged advice is used.

## 1. Setting and what is used

Challenge: s uniform in 1..r-1, input T=A_1=[s]P0 and A_i=[s^i]P0, i<=n, n=2^20. The attack reads only
B=A_d with d=2^20 (one point, 96 bytes) and T. Here r=z^4-z^2+1 with z=0xd201000000010000, and
v2(r-1)=32, so d=2^20 divides r-1. Put N=(r-1)/d (log2 N=234.857089, N even and divisible by 3)
and let E=F_r^*{}^d be the subgroup of d-th powers; |E|=N. Cheon (Eurocrypt 2006; J. Cryptology 23, 2010)
reduces recovery of s to finding v=s^d in E from B=[v]P0, then extracting s from v (section 4).
Valid range: every n with 2^20<=n<r; W does not depend on n beyond needing A_{2^20}. For n<2^20 the bound is not claimed.

The class walk uses the j=0 endomorphism of G1: the map phi(x,y)=(beta x,y), with beta a cube root of unity
in F_p, equals [lambda] for a primitive cube root lambda in F_r (the matching pair is fixed by one scalar
multiplication of P0).

## 2. Phase 1: collision search for v

Let mu6 be the sixth roots of unity in F_r. They lie in E (-1 because N is even; the cube roots
because 3 divides N). The class of a point Q is {+-phi^i(Q)} (six distinct points). Let S={[e]P0 : e in E}. Both
the tame points [e]P0 and the wild points [e]B=[e v]P0 (e in E) lie in S, and S has M=N/6 classes (log2 M=232.272127).

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

Fixed-base multiplication. The scalar e' is below r<2^254.857. Signed digits of w=37 bits, each in
[-2^36, 2^36), cover every integer of absolute value below about 2^259 in 7 windows (7*37=259; the carry out of the
top window is absorbed by the top digit, since 2^254.857 is far below 2^259). For each base and
window i<7, store the points j*2^(37i)*base for j=1..2^36 in affine form (96 bytes each): 2 bases * 7 windows *
2^36 = 14*2^36 = 962,072,674,304 = 2^39.807355 points in all, 2^39.807355 * 96 = 92,358,976,733,184 bytes =
2^46.392317 bytes (about 92 TB, 84 TiB). A multiplication is 7 table
fetches (negating y for negative digits, no fetch for a zero digit) and 6 batched incomplete affine additions
(a window w=36 would need 8 windows and 7 additions, so w=37 is the smallest window giving 7 windows; w=43 is the
smallest window giving 6 windows, 5 additions and 2^52.170 bytes of tables, and w=37 trades 2^5.8 of storage for
one more addition per step). Exceptional cases (equal or opposite operands) have probability below 2^-250
per addition, below 2^-127 over the whole run. The additions of all trails advance in lockstep rounds: a batch of
2^20 trails (of the 2^26.8 trails the run needs) is in flight at any time, with 64 bytes of state each,
2^26 bytes, and each of the 6 addition stages is done as one batch across
those trails, so the batch-inversion amortisation that the 6M convention (3M for the formula, 3M for batching)
assumes is available; the state of a trail between rounds is its scalar, step counter and type, charged below.

Table storage and transfers. The 2^46.392317 bytes of tables are one array sharded across a storage cluster. A
trail fetches table points and nothing is copied per processor, so the build is paid once, and one table serves
all trails. Each fetch reads one padded 128-byte record holding a 96-byte affine point, as 16 eight-byte loads.
The fetches are batched: the 2^20 trails in flight each issue one fetch per window, so one window stage is a
single batched remote read of 2^20 * 128 bytes = 128 MiB, and a full round of seven stages moves 896 MiB, spread
evenly over the shards. The model prices a load at one instruction and charges no latency and no block
granularity, so the transfer is 16 instructions per fetch, 112 per step for seven fetches; that charge is itemised
and priced in section 5.

## 3. Success probability

Let m be the total number of steps, split equally between tame and wild trails. In the random-mapping model
the number of tame-wild merge events is Poisson with mean (m/2)^2/M, so P(success)=1-exp(-m^2/(4M)). For
P=0.39 the mean is ln(1/0.61)=0.49430, giving m0=sqrt(4*0.49430*M)=1.40612*sqrt(M)=2^116.627787. The run uses
m=1.1*m0=2^116.765291 steps to cover walk non-randomness (Teske 2001 measured r-adding walks with r>=16 within a
few percent of random; we take 10 percent in step count, and r=2048 is far above 16). With the margin the model
gives P=0.4501, so the claim of at least 0.39 holds if the walk loses at most a factor 0.826 in collision rate.
A merge is seen at a distinguished point, so the run is not stopped the instant a collision exists: the 2^20
trails in flight each continue to their next distinguished point, on average 1/theta = 2^90 steps (the geometric
waiting time of the 2^-90 test), so the run-on adds up to 2^20 * 2^90 = 2^110 steps. That is 2^-6.765291 =
0.9193 percent of the step total and is charged at the full step price in section 5. If those in-flight samples
were written off instead of run on, the success probability would be 0.4501 * (1 - 0.009193) = 0.4460, still
above 0.39, so either treatment supports the claim; the run-on charge is the conservative one and is the one
made. The DP-set is detected by sorting or hashing and no collision is missed.
Fruitless cycles of the Bos et al. 2010 / Bernstein-Lange-Schwabe 2011 kind come from the negation
map acting as a reflection; here mu6 acts by translations of the logarithm, so none arise. A cycle needs a
sum of at most L=20*2^90 random jumps to vanish modulo M, probability about L^2/M<2^-43 per trail, below 2^-16
over the 2^26.8 trails; trails are capped at 20/theta steps and restarted, which costs at most 2^94 steps each time.

## 4. Phase 2: from v to s

Compute s0 with s0^(2^20)=v by 20 successive square roots in F_r (Tonelli-Shanks, r-1=2^32*u). Since v is
a d-th power, v=g^l with 2^20 dividing l; both square roots have exponent divisible by 2^19 at the first stage (because
(r-1)/2 is divisible by 2^31), and so on, so every choice of root succeeds.
Cost at most 2^11 F_r multiplications per root. All solutions are s0*omega, omega in mu_d. Take omega0 of exact order d
(c^N for random c, test omega0^(d/2)=-1). Baby-step giant-step over i=2^10*b+a: store X_a=[s0 omega0^a]P0 (a<2^10, fixed-base
table of P0, 6 additions each) and test Y_b=[omega0^(-2^10 b)]T (b<2^10, each a generic scalar multiplication,
at most 2^10 units counting a doubling as 2 units). A match gives s with [s]P0=T, which is the attack-side check.
Total at most 2^20.1 units.

## 5. Cost ledger (1 unit = 2000 instructions = one batched incomplete affine addition)

Per step: 6 additions = 6 units (the additions are the reference operation and are not also counted in
instructions). Other instructions, each line priced by its implementation:

| Line | Instructions |
| --- | ---: |
| Two F_p multiplications for beta*x and beta^2*x, 2*383 measured for a six-limb multiply | 766 |
| Three-way six-limb compare and y-sign normalisation of the class representative | 70 |
| Class index extraction and the 2048-entry jump-constant fetch | 20 |
| One F_r Montgomery multiplication (constants kept in Montgomery form) | 383 |
| Distinguished-point test on the low 90 bits and the step counter | 10 |
| Trail state load and store (scalar, counter, type) with loop control | 40 |
| Signed-digit recoding of seven digits at 8 each | 56 |
| Seven table fetches at 64 each | 448 |
| Operand handling around the six additions at 30 each | 180 |
| Itemised total | 1973 |

The table fetch is priced component by component. A fetch at 64: request and address 2, the transfer of the
padded 128-byte record holding a 96-byte affine point as 16 eight-byte loads 16, conditional negation of y by the
digit's sign 12, store into the addition operand slot 6, addressing and loop 8, which is 44 itemised and 64
charged; the seven fetches of a step therefore carry 140 instructions of cushion above their itemisation. The
transfer is batched over the trails in flight, as section 2 states, so the 16 loads per record are the whole
per-record transfer charge. The addition formula and its batched inversion are the reference unit and are not
charged twice.

Charged work: the itemised sum is 1973 instructions, charged as 1973, which is 1973/2000 = 0.9865 unit. A step
costs 6 + 0.9865 = 6.9865 units.

| Item | Units |
| --- | ---: |
| Steps: 6.9865 * 2^116.765291 | 2^119.569861 |
| Trails in flight running on to their next distinguished point: 2^110 steps * 6.9865 | 2^112.804570 |
| Table build (2^39.807355 batched additions, 962,072,674,304 points) with ladder heads and stores, and jump constants, charged 2^40 | 2^40 |
| Trail starts: 2^26.8 trails * (6 additions + 20 squarings) | below 2^32 |
| DP storage, sorting and matching: 2^26.8 DPs at 1000 instructions plus the sort | below 2^28 |
| Phase 2, validation of B, square roots | below 2^21 |

W(2^20) <= 2^119.569861 * (1 + 2^-6.765291) + 2^40 + 2^32 < 2^119.5831.

Arithmetic check: log2(6.9865) = 2.804570 and log2 m = 116.765291, so the step total is 2^119.569861. The run-on
row is 2^-6.765291 = 0.9193 percent of it, and log2(1 + 0.009193) = 0.013202, so the walk total is 2^119.583063.
The table build with jump constants is 2^40 units, 2^-79.58 of the total; trail starts, distinguished-point work
and phase 2 together stay below 2^32 units, 2^-87.58 of the total; the largest of them moves the logarithm by less
than 10^-22. W < 2^119.5831, and the submitted bound is log2 W = 119.59 = 11959/100, rounded up to two decimals.

Preprocessing. Preprocessing is reported as 40.0. The 2^39.807355 table points are built as 14 ladders (2 bases
times 7 windows), each advanced in lockstep across 2^18 independent sub-ladders with Montgomery batch inversion,
so one batched incomplete affine addition per point covers both the formula and the inversion; the ladder heads
add 2^21.8 points, the stores are inside the round-up, and the 2^15.5 units of jump constants (2048 * 20
squarings and 6*2048 scalar products) sit far below it, so the build is charged 2^40 units and that charge is
inside W above. A sequential build of the same point count, one dependent affine addition at a time with no batch
inversion, is at least 2^40.66 units; even charging 2^42 (4.5 units per point, above that sequential figure)
moves the logarithm of W by 2^-77.58 and leaves W < 2^119.5831, so the bound does not depend on how the build is
priced. Table-build memory transfers are loads and stores counted in the build. One table, sharded across the
machine's storage, serves all trails (trails fetch table points, nothing is copied per processor), so the build is
paid once.

Sensitivity of the two priced lines. Doubling both the record fetch (448 to 896 instructions) and the operand
handling (180 to 360) adds 628 instructions per step, for 7.3005 units and W < 2^119.6465, 0.056 bits above the
submitted bound of 119.59; the bound recorded here is the undoubled one.

## 6. Heuristics and limits

H1 (random walk): the class walk collides like a random mapping on M classes, within the 10 percent step
margin (Teske 2001; van Oorschot-Wiener 1999 for distinguished points). H2: low bits of x(Qc) are uniform, so
DPs occur at rate 2^-90. H3: tame and wild trails are independent apart from merges. The instruction counts
are analytic upper bounds, not a built implementation; loads count one instruction and latency is not
charged, as the model fixes. Memory: tables 2^46.392317 bytes (92,358,976,733,184 bytes, about 92 TB or 84 TiB),
DP store 2^33.4 bytes, other state below 2^30;
total below 2^46.4 bytes (about 92 TB). This is a large storage cluster but far short of raw-storage scale, and
2^5.8 smaller than the 2^52.170 bytes of tables a 6-window w=43 multiplication would need; the higher window
would give 5 additions per step but would move the table into the petabyte range, and w=37 was chosen to keep the
table build and the storage footprint modest. The cost model gives memory no scalar weight and counts a load as one
instruction with no latency charge, which is what the ledger uses.
Nothing larger than the toy experiment was run. This is an exploratory analytic bound, not a practical attack.
