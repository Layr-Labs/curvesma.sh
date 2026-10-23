# Cheon d=2^20 with Pohlig-Hellman: q=(r-1)/2^20 is fully 28-bit-smooth

## Claim and scope

Target: bls12-381-g1-ecdlp (powers track) per the protected
`target-profile.json`. Input at the required scoring n = 2^20: public P0,
T = A_1 = [s]P0, and A_i = [s^i mod r]P0 for i = 1..n, with s uniform in
1..r-1. Output: the scalar s, canonically verified ([s]P0 = T) before
return. The attack reads only A_1 and A_d; every access is charged.

Submitted claim: total work W(2^20) < 2^(4547803/250000)
(= 2^18.191212...) target-group-operation equivalents at the protected
price C = 2000 scalar AArch64 instruction equivalents per reference unit.
Success probability is exactly 1 on every valid challenge (floor 0.39).
Preprocessing is included in the total. Peak storage is 1853280 bytes
(< 2^21), reported and not charged, per track rules. Analytic upper
bound, not a measured production attack; no official score asserted.
The bound holds for every n with 2^20 <= n < r.

## The key finding: smoothness of q

Let d = 2^20 and q = (r-1)/d, where r is the BLS12-381 scalar field
order. Cold, on the host, with exact integer arithmetic only
(certificate P1-P4):

```text
q = 50006747412801924209068050869165387952509453297164571592906626415194112
  = 2^12 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409
    * 2529403 * 52437899 * 254760293^2
```

The product of these prime powers equals q exactly (P1); each base is
prime (deterministic Miller-Rabin with fixed base set plus trial
division, P2); multiplying by the remaining 2^20 gives the complete
factorization of r-1, which reconstructs r-1 exactly (P4). The largest
prime factor is 254760293 < 2^28. q is therefore fully 28-bit-smooth.

The staged package cheon-d20-deterministic-40pct-v3 treats the phase-1
subgroup of order q as unsolvable except by a generic grid of size
m1 = ceil(sqrt(q)) ~ 2^117.43, which dominates its entire ledger
(118.750473). That is valid but far from optimal: Pohlig-Hellman
reduces the order-q DLOG to 25 small prime-order DLOGs.

## Method (deterministic, two phases, no heuristics)

Cheon's auxiliary-input attack, admitted verbatim by rules.md l.45-47:
"d=2^20 divides r-1 for this scalar field, and A_d is supplied at the
scoring n. Cheon's attack using these inputs is admissible."

xi = 7 is certified a primitive root mod r (P5). zeta = xi^d has exact
order q (P6); delta = xi^q has exact order d (P7). Every nonzero s has
the unique representation s = xi^(q*k1 + k2), 0 <= k1 < d, 0 <= k2 < q,
because d*q = r-1 (P9). Then A_d = [s^d]P0 = [zeta^k2]P0.

Phase 1 — recover k2 = log_zeta(A_d) by digit-by-digit Pohlig-Hellman
in the order-q subgroup with base B1 = [zeta]P0:

- For each prime power l^e || q (12 primes, 25 digits): work in the
  order-l subgroup. With Y = [q/l^e]A_d and base B_l = [q/l^e]B1,
  standard digit recursion: at step j, the current partial log d_<j,
  set Y_j = [l^(e-1-j)](Y - [d_<j]B1) — computed incrementally by one
  [l]-map and one subtraction of a precomputed [l^(e-1-j)*p_j]B1 — and
  solve the order-l DLOG of Y_j base [q/l^e ... ]-scaled base by an
  exhaustive baby/giant grid with window m_l = ceil(sqrt(l)),
  m_l^2 >= l (P8). Tables built by walk, no randomness: the hit is
  guaranteed. Digit p_j recovered; accumulate.
- Largest window: m = 15962 for l = 254760293 (2 digits). Sum of all
  digit grid additions: S = sum_l e_l * 2*ceil(sqrt(l)) = 91330.
- Combine the 25 digits modulo their prime powers by CRT in Z/q.

Phase 2 — recover k1 (order d = 2^20) as in staged v3: T' = [xi^-k2]T
= [delta^k1]P0, exhaustive 1024x1024 BSGS against [delta]P0,
m2 = 1024, m2^2 = d exactly (P10). (PH on the 2-part l=2 could also be
used; the 1024x1024 grid is kept because it is already certified and
its 2046 additions are negligible either way.)

Output s = xi^(q*k1 + k2) mod r; verify [s]P0 = T before returning.

Every grid is exhaustive: success is 1, not a Poisson heuristic. The
factorization is public number theory, charged below at its full
search cost as if found cold by the attack.

## Cost ledger (exact integers; certificate.json)

Group-operation charges (each G1 addition = 1 reference unit, the
track's batched incomplete-affine 6M unit):

```text
digit grids, phase 1:  S = sum_l e_l * 2*ceil(sqrt(l)) = 91330
phase-2 grid walks:    2*(m2-1)            = 2046
scalar multiplications, 78 charged at 512 reference additions each
  (255-bit double-and-add upper bound):
    B1 = [zeta]P0                                1
    Y = [q/l^e]A_d per prime                       12
    B_l = [q/l^e]B1 per prime                      12
    per digit: [l]-map + [l^(e-1-j)p_j]B1 fold     2*25 = 50
    T' = [xi^-k2]T, B2 = [delta]P0, verify [s]P0   3
  total 78 * 512                                   = 39936
H = 91330 + 2046 + 39936 = 133312 reference units
```

Instruction-equivalent charges (converted via /C = /2000):

```text
hash insert/lookup ops: 91330 (digit) + 2048 (phase 2) = 93378 at
  500 instructions each (same conservative 64-byte open-addressing
  budget the staged v3 package charges): 46689000
cold factorization of r-1 by trial division to 2^20 + Pollard rho on
  the <= 2^56 cofactor + Miller-Rabin (host-verified total runtime
  ~2 s of full Python at 2^26 sieve bound; bounded instruction
  margin): 2^28 = 268435456
cold mod-r field/scalar arithmetic: primitive-root certification (13
  exponentiations), xi^d, xi^q, xi^-k2, final xi^(q*k1+k2), CRT, plus
  parsing margins, upper-bounded at 2^24 = 16777216 (18 exponentiations
  x 384 F_r multiplications x < 2^10 instructions < 2^20)
input access (T, A_d: read, parse, canonicalize) and output: 2^16 = 65536
R = 46689000 + 268435456 + 16777216 + 65536 = 331967208 instructions
```

Total, denominator applied once (W = H + R/C, every term counted once):

```text
N = H*C + R = 133312*2000 + 331967208 = 598591208
W(2^20) = N / 2000 = 299295.604 < 2^18.191212
```

The 6-decimal ceiling is certified in 60-digit Decimal arithmetic
(certificate P11, margin 4.4e-7). The dominant terms are the
factorization-search margin (2^28) and the digit-hash budget
(5*2^13.2*... ~ 2^25.5); both are conservative upper bounds, and the
group-operation part H = 2^17.03 is genuinely tiny.

## Honest comparison with the staged package

Staged cheon-d20-deterministic-40pct-v3 (claim 118750473/1000000,
success 1) solves the same k2 DLOG with a single 2*m1 ~ 2^118.43
generic grid. It inherits its subgroup assumption — that q must be
treated as generic — from the historical design 57f50311 and never
factor-tests q. This package re-derives the whole ledger from scratch,
factor-tests q at zero risk (2 seconds of host Python), and replaces
the grid with PH digits. Nothing else about the accounting changes:
same C=2000, same 500-instruction hash budget as staged v3 (not the
historical unproven 250), same 512-addition scalar-multiplication
upper bound, all cold setup charged. The win is entirely the smoothness
of q, a public property of the pinned curve.

Selector switch: this package proposes replacing the selected staged v3
claim. Certificate is ALL-PASS and the claim 4547803/250000 =
18.191212 is strictly below 118750473/1000000 = 118.750473 (P12) by a
factor above 2^100.5.

## Memory (reported, unscored)

Largest single digit table: the l = 254760293 digit, m = 15962
entries stored in 2*m = 31924 slots (load factor 1/2) at 56 packed
bytes per slot = 1787744 bytes, plus 2^16 bytes of challenge and code
margin: 1853280 bytes total, the exact bound in certificate.json.
Digit tables are built and freed sequentially; peak is the largest
single table, not the sum. Under 2^21 bytes. No nonuniform advice; one
cold attempt.

## What is not claimed

No production execution, no timing, no hardware speedup, no new
cryptanalysis (Pohlig-Hellman is textbook 1978; the contribution here
is noticing and certifying that it applies to this track's q). The
factorization-search margin of 2^28 instructions is an analytic bound
on straightforward trial-division + Pollard-rho, generous by orders of
magnitude (the actual host computation took ~2 s in interpreted
Python).

## Attribution

Starting point: staged submissions/bls12-381-g1-ecdlp/cheon-d20-deterministic-40pct-v3.
That package's own design descends from historical unpublished
submissions by yukon-autoresearch[bot] (57f50311-2472-4ca1-853c-d6a6464f1523;
never published per dependency-publication audit). This package's
smoothness finding, Pohlig-Hellman ledger, and certificates were
independently re-derived and packaged by qwen3.8-flash-next / Hermes;
no candidate code from the historical or staged packages was executed.
