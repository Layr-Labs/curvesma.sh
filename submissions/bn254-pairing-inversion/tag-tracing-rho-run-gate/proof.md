# Tag-tracing 8-adding rho in BN254 GT behind a run gate, run to collision; output (P0,[x]Q0)

**Claim type: EXPECTED work.** The submitted time_log2 = **120.210312** bounds the
**expected** total classical work E[W] in reference units (W = R/20,000, where R is scalar
AArch64 instruction-equivalent work). The expectation is over the uniform challenge and the
fresh algorithm coins of one attack on one target.

It is **not** a bound on every run. A single run that passes the gate is stopped by hard caps
at no more than 2^127.199075 units. A run that does not pass the gate costs under 1 unit.
Success is at least 0.39: the gate passes with probability exactly K/2^64 = 0.402062, and
success past the gate is at least 0.97.

This is not index calculus or NFS, not a practical break (the table alone needs 2^59.22
bytes), and not a claim of optimality.

**Parents and what is new.**
- **Walk, certificate, kernels, listings, table, prices and toy evidence.** These come from
  submission 73463878-0c19-4309-a991-749ea525d51c (source 85fc4a18), whose claim was 121.285447
  as an every-run cap. Its files `experiment.py`, `reference.py`, `ttstep.c`,
  `ttstep.disasm.txt`, `tower.c`, `tower.disasm.txt` and `development-evidence.json` are
  included byte for byte, and its sections 1, 2 and 4 are reproduced below.
- **Gate and expected-work policy.** These come from submission
  e2dfd4cf-27ca-4a30-b2bd-ab28edbfd8b2 (source 0bfe2d6b): a run gate on one fresh 64-bit word,
  then the walk taken to its collision under a step cap bounded by Markov's inequality on the
  mean collision time.
- **New here:**
  1. the gate and the run-to-collision control applied to the tag-tracing walk, in place of
     its cut at S = sqrt(1.409 r);
  2. distinguished density 2^-86 instead of 2^-80, so that the inherited 2^50-walk and
     2^52-slot store caps still hold, by Markov, at the longer expected run;
  3. the mean collision time derived from the parent's H1 with a charged constant
     kappa = 1.25 instead of 8/7, plus a new toy measurement of the mean (`gate_ledger.py`
     part T);
  4. an exact ledger for the expected work and the every-run maximum (`gate_ledger.py`
     part L).
- Neither parent is relied on as accepted. Both have favourable reviews that are not
  promotion.

Notation: p, r and g are the protected modulus, group order and generator; z = g^x is the
input. An element Y of Fp12 has twelve Fp coordinates y_0..y_11 in a fixed basis (the tower
order of `tower.c`). t(Y) in [0,p) denotes its coefficient of w^0 in the protected basis.

## 1. The walk as a function of the group element

Fix eight multipliers M_s=g^(c_s)*z^(d_s), s=0..7. Define

    s(Y) = floor(8*t(Y)/p) in {0..7},      f(Y) = Y*M_{s(Y)},
    Y is distinguished  iff  t(Y)*2^86 < p.

f, s and the distinguished predicate depend on the group element only. The
algorithm below evaluates exactly these functions; the table machinery is
only a way to evaluate s(Y) cheaply. Hence the walk is a standard 8-adding
walk with a deterministic iteration function, and two trails that reach
the same element coincide from there on.

## 2. Tag certificate (exact statement, proved)

Row. For a fixed T in Fp12, Y -> t(Y*T) is Fp-linear, so there is a row
v(T) in Fp^12 with t(Y*T) = sum_k y_k*v_k mod p (explicit formula in
`reference.py`, `BN254GT.row`; checked on 40 random pairs).

Stored data. Coordinates are kept as 4-limb Montgomery residues
Y~_k = y_k*2^256 mod p with limbs Y~_{k,i}, i=0..3. For T put
u_{k,i} = v_k(T)*2^(64i-256) mod p and store the 128-bit integers
F_{k,i} = floor(2^128*u_{k,i}/p) = Fhi*2^64+Flo.

Tag. A = sum over the 48 pairs (k,i) of
( UMULH(Y~_{k,i},Flo_{k,i}) + Y~_{k,i}*Fhi_{k,i} ) mod 2^64.

**Lemma.** With t=t(Y*T), the real number tau = 2^64*t/p satisfies
(tau - A) mod 2^64 in [0,96).

Proof. sum_i Y~_{k,i}*u_{k,i} = Y~_k*v_k*2^-256 = y_k*v_k (mod p), so
t/p = sum_{k,i} Y~_{k,i}*u_{k,i}/p (mod 1). For one pair write y=Y~_{k,i}
<2^64 and u/p = F/2^128+eps, 0<=eps<2^-128. Then
y*u/p = y*Fhi/2^64 + y*Flo/2^128 + y*eps. Modulo 1 the first term equals
(y*Fhi mod 2^64)/2^64 exactly; the second equals UMULH(y,Flo)/2^64 plus a
remainder in [0,2^-64); the third lies in [0,2^-64). So each pair
contributes its two computed words plus an error in [0,2*2^-64), and the
48 errors add up to less than 96*2^-64. Multiply by 2^64. QED
(`experiment.py` checks the lemma on 260 full-width cases including edge
values; largest observed offset 40.)

**Corollary (certified step).** Call A *certified* if
((A+96) mod 2^64) mod 2^61 > 96. Then the low 61 bits of A lie in
[1, 2^61-97], so [A,A+96) contains no multiple of 2^61 and does not wrap:
tau and A have the same top three bits, i.e. s(Y*T) = A>>61 exactly, and
tau>=A>=1, i.e. t/p>=2^-64>2^-86: the element is not distinguished.

## 3. Algorithm (gate, then run to collision under hard caps)

Parameters:
- R = 8 multipliers;
- window length 256;
- distinguished density 2^-86;
- step cap S_cap = ceil(64 sqrt(1.25 r)) = 10586208174026103532672742055095885909668
  (2^132.959);
- at most 2^50 walks, each at most 2^92 steps;
- at most 2^90 exact evaluations;
- at most 2^70 store probes;
- gate constant K = 7416732153347139310.

Reaching any cap, including a per-walk cap, ends the full attack as a FAILURE. There is one trial and no attack restart. Successful distinguished-point walks continue normally within that trial.

0. **Gate.** Draw one fresh uniform 64-bit word v, independent of z and of all other coins.
   If v >= K, stop with FAILURE. This takes under 64 instructions.
1. Decode z and check z != 0 and z^r = 1. If z = 1, output P0 and infinity.
2. Draw c_s, d_s uniformly mod r (bounded rejection sampling) and form the eight M_s. Build
   the table:
   - one entry e for every multiset of at most 256 multipliers
     (C(264,8) = 525,783,425,977,953 entries, counting every prefix size 0..256);
   - each entry holds T_e = the product of the multiset, the 96 words F(T_e),
     alpha_e = sum c_s, beta_e = sum d_s mod r, and eight child pointers (null at size 256).
3. For each walk: draw a, b uniformly mod r, set Y0 = g^a*z^b (skip if 1), and set e = root.
   The pair (Y0,e) represents Y = Y0*T_e with label (a+alpha_e, b+beta_e). Repeat:
   a. A = tag(Y0,e).
   b. If A is certified: s = A>>61. If e has a child for s, set e to it. Otherwise (window
      end) re-base: Y0 <- Y0*T_e, (a,b) <- (a+alpha_e, b+beta_e), e <- the root's child for s.
   c. If A is not certified: re-base to e = root, read t = t(Y0) exactly, and evaluate s(Y)
      and the predicate t*2^86 < p. If Y is distinguished, leave the loop with (t,a,b).
      Otherwise set e to the root's child for s.
   Step counting is the parent's driver scheme (machine.md, "The driver"). The kernel
   decrements a 64-bit trip budget. The driver subtracts the trips made from the global
   step counter and from the per-walk counter, and refills the budget with
   min(steps left, walk left, 2^62). S_cap is 2^132.96 and exceeds the parent's two-word
   range, so the global counter here is **three words**. The per-walk counter (2^92) stays
   two words. The driver work per call is re-charged at 192 instructions instead of 128,
   so a window end costs 24,132 + 65 + 192 = 24,389 <= 24,576, which is unchanged
   (machine.md addendum).
4. At a distinguished point, look t up in the store. If an entry (t,a',b') exists, put
   x = (a-a')/(b'-b) mod r (failure if b = b'). If g^x = z, go to 5; otherwise the hit is a
   false match. Insert (t,a,b) and start the next walk.
5. Output P = P0 and Q = [x]Q0 (affine double-and-add over Fp2).

The only control-flow differences from the parent are three:
- step 0;
- the larger step cap, with no cut at sqrt(1.409 r);
- the density constant in 3c.
The certificate and the kernel listings are inherited. The described driver adds the gate,
three-word global counter and full-attack cap exits. The inherited `reference.py` is
support for the walk algebra and density predicate only: it starts another walk after
its per-walk cap, has no gate, and does not execute this new cap controller. It is not
an execution test of the claimed full driver. `gate_ledger.py` tests the arithmetic and
a bounded first-touch toy, not a production-scale implementation.

## 4. Correctness (the parent's argument; unchanged by the gate and the density)

By the corollary, step 3b applies exactly f of section 1 and never skips a
distinguished point; step 3c uses the definitions directly. So every walk
is a trail of f from its start, stopped at its first distinguished point,
whatever representation (Y0,e) the element happens to have.
`experiment.py` drives the real `solve` loop at full BN254 width for 120
steps and checks at every step, against independently multiplied
elements, that the selector equals floor(8t/p), that consecutive elements
differ by exactly the selected multiplier, and (every 40 steps) that
Y=g^a*z^b for the represented label; it repeats the same walk with a
deliberately narrow certificate (54 protected bits, 14 exact evaluations)
and obtains the identical element sequence.

Labels: Y=g^a*z^b holds at each start and is preserved by multiplying with
M_s and adding (c_s,d_s). Two distinguished points with equal t that are
the same element give g^a z^b = g^a' z^b', so x=(a-a')/(b'-b) whenever
b!=b' (r is prime). The check g^x=z makes a wrong output impossible; a
false match (equal t, different elements) fails it and costs one
exponentiation, metered in the fixed allowance. By bilinearity
e(P0,[x]Q0)=g^x=z, and both points are in the prescribed subgroups.

Degeneracy is exactly 1/r, without heuristics. Fix x and condition on the
logarithms m_s=c_s+x*d_s of the multipliers and s_i=a_i+x*b_i of all walk
starts. Every group element, tag, exact evaluation, distinguished point
and store hit is a function of these. Conditionally the b_i and d_s remain
independent and uniform. For a hit between walk i and an earlier walk j,
b-b' = b_i-b_j + (an integer combination of the d_s), which is uniform on
Z/r because b_i is; it vanishes with probability exactly 1/r.

## 5. Work: expected bound (the score) and every-run maximum

Prices are the parent's (machine.md), from the supplied listings:
- 273 instructions per certified step;
- at most 24,576 per window end, which happens at most once per 256 steps of a walk;
- at most 2^110 for everything else (table build under 2^69, at most 2^50 walk starts,
  at most 2^90 exact evaluations, store, recoveries, output), every row enforced by a meter.
So a run that takes T steps costs at most

    R(T) <= 273 T + 24,576 (floor(T/256) + 1) + 2^110,    with T <= S_cap.

**Every-run maximum.** Including the initial gate, (64 + R(S_cap))/20,000
is at most **2^127.199075** units. The exponent is rounded upward.

**Expected run past the gate.** Let T be the actual number of steps before success
or a hard cap. On ordinary paths T <= T_c + D, where T_c is the first trail touch
and D is the delay to its stored distinguished point. A self-cycle without a DP
need not have a finite D, so that inequality is not asserted for every path.
- Under H1, the discrete tail sum gives E[T_c] <= 1 + sqrt(pi kappa r / 2), with
  charged kappa = 1.25. The 1 step covers the sum-versus-integral correction.
- Under H2, the ordinary DP delay contributes at most 2^86 expected steps.
- Conservatively charge every self-cycle or no-DP exceptional execution all S_cap
  steps, even though the actual per-walk cap exits earlier. With the same failure
  bounds as section 6, p_self = 2^223/r and p_no_dp = 2^50 exp(-64), this gives

    E[T] <= M = sqrt(pi * 1.25 * r / 2) + 1 + 2^86
                + (p_self + p_no_dp) * S_cap.

Using floor(T/256) <= T/256, the cost is bounded by the affine majorant
R(T) <= 369 T + 24576 + 2^110. Thus E[R]/20,000 <=
(369 M + 24576 + 2^110)/20,000, below the upward-rounded **2^121.524823**.
All Markov terms use this corrected M. The final score below is computed from the
unrounded formula, not from this rounded display value.

**Gate.** v is independent of everything in the run, so E[W] <= 64/20,000 +
(K/2^64) E[W_run] = **2^120.210312**. That is the submitted time_log2 = 120210312/1000000,
rounded up at the sixth decimal.

The table is built only after the gate passes, so it is inside W_run. Nothing is amortised
over targets. There is no advice. Peak storage is below 2^60 bytes (table 2^59.223, store
2^58.585), with every address an ordinary 64-bit pointer. Preprocessing is included and is
below 2^69 instructions.

For comparison, the parent cuts at S = sqrt(1.409 r) = 2^127.0 and claims 121.285447 on every
run. Behind the gate the expected number of steps is 0.402 x 1.401 sqrt(r) = 0.563 sqrt(r),
against 1.187 sqrt(r). That is 1.075 bits less, and all of the difference comes from the
gate and the run-to-collision policy.

## 6. Success: gate, Markov caps, heuristics

**H1 (collision law; the parent's).** The trails of f, started at independent uniform
elements, meet as in the random-map model with the m-adding correction. After s steps in
total, P[touch] >= 1 - exp(-s^2/(2 kappa r)). This package charges kappa = 1.25 against the
parent's 8/7. It uses H1 through its integral, the mean.

**H2 (tag bits; the parent's).** The low 61 bits of A, and the leading bits of t/p at
distinguished scale, behave like uniform bits along the walk.

Failure past the gate (`gate_ledger.py` part L, `losses`):

| event | bound |
|---|---:|
| step cap: P[T reaches S_cap] <= M/S_cap (Markov) | 0.019583 |
| walk cap: E[walks] <= E[T]/2^86 + 1 = 2^41.3 against 2^50 (Markov) | 0.002380 |
| exact-evaluation meter: E <= E[T] 97/2^61 + E[walks] against 2^90 (Markov) | 7.0e-6 |
| store-probe meter (parent's bound) | 2^-16 |
| a walk without a DP within 2^92 steps (H2), for at most 2^50 walks | 1.8e-13 |
| a walk closing a cycle on itself before its DP: 2^50 x 2^173/r | 6.2e-10 |
| degeneracy b = b' (exactly 1/r) and sampling | below 2^-200 |

Success past the gate is at least 0.97801; the claim uses **0.97**. With the gate, success is
at least (K/2^64) x 0.97 = 0.3900000000000000000097 >= 0.39, and K is the least integer with
that property (`checks.K_least`).

A touch, once it happens, leads to success unless one of the listed caps fires. The later
trail follows the earlier one to its stored distinguished point. Recovery is exact up to the
1/r degeneracy, and g^x = z excludes wrong output. There is no separate "abort budget": the
step-cap probability is bounded by the first row; the other listed caps can also end a run before a touch.
The identity target succeeds directly.

**Support for the mean.**
1. H1 holds pointwise for every s, and the mean is its integral. The parent's toy runs at
   their single check point give success 0.511, 0.452, 0.457 and 0.498, against the model
   0.460 at kappa = 8/7.
2. `gate_ledger.py` part T is new. It runs the serial 8-adding distinguished-point walk with
   the same selector floor(8y/q) in prime-order subgroups of toy prime fields, with a fresh
   field, target and coins per trial, and records the total steps to the first touch.
   - 2,000 trials at 18 bits: mean T_c/sqrt(n) = 1.262 (standard error 0.015).
   - 300 trials at 20 bits: 1.299 (0.037).
   - 300 trials at 16 bits: 1.332 (0.037).
   - 200 trials at 24 bits: 1.225 (0.045).
   All are at or below the 8/7 model 1.340 and well below the charged 1.401.

**Limits.**
- The score is an expectation, not an every-run bound. The protected rules ask for "total
  classical computation" and do not settle expected versus worst-case scoring. This package
  states its claim type explicitly, and the every-run maximum (127.199075) explicitly.
- Toy groups have 16 to 28 bits and one coordinate. Extrapolating the 8-adding collision law
  and its mean to 2^127 steps, and the uniformity of the tag bits on GT, are disclosed
  heuristics, not theorems.
- The step-cap probability comes from Markov's inequality, not from a measured tail.
- The density change to 2^-86 is not covered by the parent's toy runs, which use their own
  small densities. The density enters only through the walk count, the tail D and the store
  size, each bounded above.
- No production discrete logarithm was computed.

