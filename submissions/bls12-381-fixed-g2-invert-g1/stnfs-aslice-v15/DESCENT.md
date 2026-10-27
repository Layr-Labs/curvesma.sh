# The window-sieved descent (v14)

This replaces the v13 descent. The v13 recurrence charged a 2^26-instruction exact
norm screen on every coefficient tuple of every special-q box; the advisory
diagnostics priced that honestly at two-log work 2^126.37 units (CF128/CG24), which
exceeded the whole claim. The repair is not a constant: it is a different node
algorithm, specified here, priced by measured kernels, and validated by exhaustive
toy-scale differentials. The v13 parameters (CF64/CG12 in code, CF128/CG24 in prose,
the +160/+28 ECM operand offsets and the -32-bit survivor slack) are withdrawn.

## Node structure (every node, including the boot)

A node is (side, q, root pair (s0,t0), smoothness target n < log2 q, box k). Let L_q
= { (a,b) in Z^12 : a(s0) - t0 b(s0) = 0 mod q }, determinant q. The production basis
is LLL-reduced (LLL12 = 2^36 instructions per node, certification as in v13: exact
Gram-Schmidt size reduction and Lovasz condition, Gram determinant q^2). The box is
c in [-k,k]^12 with cardinality 2^S = (2k+1)^12, S on the grid 12 log2(2k+1) for
k in {2,3,4,5,6,8,10,12,16,24,32,48,64,96,128}; it is processed in bounded segments
(2^16-cell sum arrays), never as a whole-box allocation.

1. **Certified envelope screen.** Along c0-lines (the shortest reduced basis
   direction), a and b are affine in c0 with integer coefficients, so the v13
   embedding-bound machinery applies verbatim with interval a and b: per place the
   degree-12 homogenised f-factor and the g-factor are evaluated in outward-rounded
   binary64 intervals over 64-cell blocks, adaptively subdivided to width 128
   (8 bits). Blocks whose certified upper bound (log2 N_f + log2 N_g, 16-scaled)
   exceeds the envelope are dropped before sieving; retained cells have true norms
   at or below the envelope, so crediting proba(envelope) per retained cell is
   conservative. Envelopes: log2 N_f <= 72c + 76 and log2 N_g <= 12c + 6 log2|u| + 9,
   c = (log2 q + S)/12, before dividing out q. The measured distribution over actual
   reduced special-q lattices (16 ideals, q in [2^70,2^260], 8 box sizes, 6144
   exact norm pairs, plus the v13 study's 19,200 pairs, maxima 94.0/13.2 and
   102.4/17.0) supports a retained fraction of at least 0.99 at these envelopes
   (stats14.json; cf = log2 N_f - 72c: mean 2.3, p99 64.4, p99.9 78.6); the
   charged ENV_FACTOR is the conservative 0.97.
   ENV_FACTOR = 0.999 multiplies the credited yield.
2. **Bucket sieve.** For each sieving prime p <= 2^LX (LX = 69), p != q, and each
   root pair (s',t') of f or g mod p (about NBAR = 2 per prime over both sides),
   the divisibility condition is the linear congruence w.c = 0 mod p with
   w_i = a_{B_i}(s') - t' b_{B_i}(s') mod p, the basis images evaluated at s'.
   The hit set {c in box : w.c = 0 mod p} is the full lattice L_p = ker(w mod p)
   (determinant p) intersected with the box, enumerated by a Schnorr-Euchner walk
   over the reduced basis of L_p. Per (p, root pair): the systematic basis
   (p e_1, e_i - (w_i/w_1) e_1), 12-dimensional reduction and functional
   composition, PAIR_SETUP = 2^24.8 instructions (the measured sound 12-dimensional reduction: 11.6M, 21.0M, 30.1M Ir at p = 2^20, 2^40, 2^62); per hit: the walk emits
   (segment, offset, 16 log2 p) into the segment's bounded record buffer;
   per cell the whole sieve (walk, records, apply, scan, bound screen) is charged
   SIEVE_CELL = 2^21.64 instructions per cell: the measured 2946 Ir per leaf
   (sound build, callgrind) times the smallest measured leaves-per-hit (254, at
   p = 2^20 where the hit mass lives; the measured range over p = 2^20..2^62 is
   254..1663), plus the record append, times 4.32 hits per cell, plus 10.9
   zero/apply/scan and 1200 for the certified bound screen (a widened allowance;
   KERNELS.md). The production leaves-per-hit is the one extrapolated quantity
   (see the recurrence section for both-end sensitivities).

3. **Threshold and survivors.** A cell is a survivor if its unsieved residual per
   side (certified norm upper bound minus the sieve-found mass) is at most
   k n + 96 bits. The 96-bit slack covers the sub-2^8 valuations, unsieved prime
   powers (the sieve credits one log2 p per prime), and the 16-scaled rounding.
   Survivors' exact norms (2^26 instructions) are divided by the sieved primes
   and the residual (at most k n + 96 bits) is factored by ECM at level
   ceil(n / 3.32) digits: a relation requires every residual prime <= 2^n and at
   most k per side; the relation's children are exactly its residual primes in
   (2^LX, 2^n]. k in {2,3,4,6} is chosen per node by the recurrence. Survivors are
   processed until STOP_FACTOR x MIN_RELS relations are found and verified.
4. **The corrected window yield.** A relation at smoothness target n with at
   most k window primes above the sieve bound requires the sub-2^LX part of the
   norm to be 2^LX-smooth; the yield is the semismooth probability
   pr = sum_{j=0}^{k} (L^j / j!) * rho((log2 N - j n + alpha) / LX) per side,
   L = ln(n/LX), each window prime charged at its maximum mass n bits (the
   conservative form of the standard large-prime-variant estimate; the reviewer's
   formula of the v14 rejection). At n = LX the window is empty and this reduces
   to the plain Dickman probability used at the bottom bins. The v14 model
   (plain n-smoothness times a Poisson CDF on the window count) ignored the
   sub-2^LX smoothness cost and is withdrawn. The false-positive rate is
   rho((nf - k n - 96)/LX) rho((ng - k n - 96)/LX) per cell against the true rate
   proba(nf, alpha_f, n) proba(ng, alpha_g, n) P(<= k window primes). With loose
   windows (k n approaching the norm) the threshold admits almost-smooth junk at
   a rate that drownds the ECM. The junk model rho((nf-kn-96)/LX) rho((ng-kn-96)/LX)
   is an overestimate of the survivor rate in simulation (so the survivor ECM is
   over-charged, never under-charged), and the recurrence charges the residual ECM
   on min(all survivors, the stopping bound). The toy differentials (evidence14
   [8] and [8b]) confirm the machinery: exact congruence functionals, exact
   scores, no false negatives over the exhaustive thin-window box, no false
   credits and exactly extracted children on true relations through the full
   survivor path, with unfactored survivors conservatively rejected exactly as
   the production ECM does.
5. **Children and recursion.** Children of a passing relation land in bins
   [lB, n); their distribution is the kappa-weighted Dickman model conditioned on
   at most k window primes (scaling pois_cdf(lam, k-1)/pois_cdf(lam, k)). The conditional child means are capped at k window primes (scale
   min(1, k/E[#unconditional])). The subtree cost uses PAD_CHILD = 2 and
   PAD_NODE = 4 as in v13. Failed branches
   backtrack to the next relation of the parent node; nodes hold MIN_RELS = 8
   expected relations before the box is accepted by the recurrence.

6. **No small-unit duplicates were found.** The lattice L_q is closed under the
   unit group of O_h, and a unit multiple of a relation is the same relation
   (identical norms and ideals), so the box could in principle contain duplicate
   credits. Measured on actual reduced special-q lattices (q = 2^120 g, 2^200 f,
   2^70 g, boxes k = 96): the in-box rate of the mates by every one of the 502
   units of O_h with coefficient vectors in [-2,2]^6 (found by exhaustive
   enumeration of norm-+-1 sextics; y^+-1,+-2,+-3 included) is zero over 120,480
   sampled mate tests (stats14.json). The unit action's anisotropy in the reduced
   basis coordinates pushes every mate out of the box; units with larger
   coefficient vectors scatter strictly further. O_h has signature (2,2), so the
   unit group has rank 3 and larger units exist and are untested; the residual
   duplicate risk from them is inside H5's allowance.

## Boot

The boot is one p^6-determinant lattice (a B + b A = 0 in F_p[y]/h for the batch
target T' = T^e, e uniform; the log of T' divided by e mod r gives log T; refresh
after at most 2^32 trials). The boot element a - b theta has a single norm (the g-side form a^2 + |u| y b^2 over the places of h: there is no f-side for theta^2 = -|u| y), so its children are g-side nodes and the recurrence descends the g side only. Its norm envelope is 12 c + 6 log2|u| + 9 with c = (6 log2 p + S)/12, which at the operating box (S = 91.1, scale 96) is 2786.4 bits, covered by the 2800-bit envelope; the experiment measures boot norms at scales 32 and 96. A single
2^91-cell box at scale 96 replaces the v13 design of 2^53.6 ECM candidates: the
boot is sieved like any node at the corrected yield (lB1 = 160, k = 12; the
ECM-per-trial boot is the evaluated control at every lB1). The 2^36-instruction batch
setup (lattice construction, certification, T^e) is unchanged from v13. Boot norms
at scale 96 are bounded by the envelope (the supplied study measured 2048 exact
norms at scale 8, max 2725.9; the experiment measures further scales).

## ECM and the honest no-sieve control

ECM prices are the v13 model (GMP-ECM table, measured ecm12 Montgomery
multiplication, 99% per-prime success). The full schedule at level ceil(n/3.32)
digits on every trial is honest: after it, every remaining cofactor prime exceeds
2^n, so the trial is decided. The recurrence retains this control at every bin
(e.g. 2^139.82 two-log at CF128/CG24 without sieve).

## Recurrence and expectation

descent14.solve: 2-bit bins, the k-grid, S-grid and target grid as above; per node
the cheaper of the ECM control and the window-sieve (the root-pair count for the setups is 2 per rational prime over both sides, NBAR = 2, so pairs = 2 li(2^LX)); children bins as in 5. At the
operating point (CF 76, CG 9, ENV_FACTOR 0.97, LX 69, MIN_RELS 8, SIEVE_CELL
2^21.64, PAIR_SETUP 2^24.8) the recurrence gives **two-log work 2^105.91 units**
(boot own 2^112.7 instructions, children subtree 2^123.7; at the unverified optimistic leaves/hit = 8 the recurrence gives 2^105.42, and at the worst measured end leaves/hit = 1663 it gives about 2^111.5).

**H5 (expectation allowance).** The recurrence is a model. The explicit allowance
of 128 for backtracking, repeated or rejected relations, correlated children, ECM
variety, rank/unknown-log dependencies and implementation overheads gives expected
two-log work below 2^112.91 units. The descent meter caps both logs at **2^116
units**; Markov gives P(cap exceeded) <= 2^-3.09, and the modelled overall success is
at least 0.9 (1 - 2^-3.09) = 0.794, above the required 0.39. The differential
evidence ([7], [8]) supports the node mechanics; the finite-tree viability beyond
MIN_RELS and the pads remains the disclosed heuristic limitation, as in v13.
