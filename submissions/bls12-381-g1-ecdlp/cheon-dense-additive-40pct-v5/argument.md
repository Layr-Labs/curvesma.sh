# Cheon d=2^20 with honestly-priced geometric chains: every step on encoded points is a public-scalar multiply, never a fake addition or a scalar exponentiation

## Claim and scope

Target: bls12-381-g1-ecdlp (auxiliary-powers track) per the protected
`target-profile.json`. Input at the required scoring n = 2^20: public
P0, T = A_1 = [s]P0, and A_i = [s^i mod r]P0 for i = 1..n, s uniform in
1..r-1. Output: the scalar s, returned only after the verification
[s]P0 = T.

Submitted claim: total work W(2^20) < 2^(12752939/100000)
(= 2^127.529390) target-group-operation equivalents at the protected
price C = 2000. Success probability is exactly 1 on every valid
challenge (floor 0.39). Preprocessing is included. Peak storage is
below 2^125 bytes, reported and unscored. Analytic upper bound; no
production run, no timing, no official score asserted. The bound holds
for every n with 2^20 <= n < r; the attack reads only A_1 and A_d and
charges every access.

The claim strictly beats the ACCEPTED frontier on this track
(cold grid package cheon-d20-hardened-40pct-v4, score 134/1, submission
d73e5e61) by 6.470610 bits.

## Why a rebuild was forced: defect map of the two rejected packages

Both rejected packages mis-modeled what can be computed on the encoded
points A_i = [a_i]P0. The encoding a -> [a]P0 is a group homomorphism
(Z, +) -> (G1, +): it preserves ADDITION of scalars. It does NOT
preserve multiplication: there is no G1 operation that turns [a]P0 into
[a^2]P0 for an unknown a.

- v4 (cheon-ph-smooth-cofactor-40pct-v4, claim 18.191212, rejected
  04:08 UTC): ran multiplicative Pohlig-Hellman on the ENCODED scalar.
  From A_d = [zeta^k2]P0 it formed [q/l^e]A_d = [(q/l^e) zeta^k2]P0 and
  read an order-l digit off it. Multiplicative PH needs
  [(zeta^k2)^(q/l^e)]P0 — the scalar EXPONENTIATED — which is not any
  [c]A_d. Its digits were meaningless. Toy demonstration
  experiments/grid_bsgs_toy.py T3: [c z] != [z^c] on 0/4000 probes.
- v3 (cheon-d20-deterministic-40pct-v3, claim 118.750473, rejected):
  its phase-1 ledger charged one G1 addition per baby/giant step:
  baby chain B, B+G, B+2G, ... with G = [zeta]P0. But B + vG =
  [(1 + v zeta)]P0 — an ARITHMETIC progression 1 + v zeta, not the
  GEOMETRIC progression zeta^v the algorithm needs. A geometric chain
  [zeta^{v+1}]P0 = [zeta * (zeta^v)]P0 requires multiplying the
  RUNNING POINT by the known scalar zeta: a double-and-add of about 255
  doublings + 85 additions, NOT one addition. The 118.75 grid price was
  therefore under-charged by roughly a factor 2^8; the honest grid is
  ~2^127.4, which is what this package states. Toy demonstration T2:
  the arithmetic progression lands on the true encoded scalar with
  frequency ~0.2% and yields a wrong digit on every one of its hits.

This package is the same mathematical attack (Cheon two-phase BSGS; the
track rules admit it explicitly) with the chain cost priced honestly.
No operation in the attack is anything other than: (i) a G1 addition of
two known points, (ii) multiplication of a known point by a KNOWN
integer (double-and-add, charged in full), (iii) field arithmetic on
known integers, or (iv) storage/lookup. Direction D2 from the rebuild
brief (linearized vanishing-polynomial digit filters) was evaluated and
rejected on its own arithmetic: testing k2 mod l via a degree-<= n
vanishing polynomial requires producing [ (s^i)^l - s^i ]-type
combinations, but s^i is unknown — [a]P0 gives no handle on [a^l]P0
except through the supplied powers themselves, and the only supplied
exponents are 1..2^20; for l > 2^20 (all four large prime powers of q)
no such filter exists, and for the small ones a filter would still need
[zeta^{k2 l}]P0 vs [zeta^k2]P0 discrimination that costs exactly the
grid it saves. The full ledger below is what survives.

## Method

Cold setup (all charged): factor r-1 by trial division to 2^20 plus
Pollard rho on the remainder; certify r prime; find xi = 7 as primitive
root by the standard test over the certified factorization; set
d = 2^20, q = (r-1)/d, zeta = xi^d (exact order q), delta = xi^q
(exact order d). Every nonzero s has the unique decomposition
s = xi^(q k1 + k2), 0 <= k1 < d, 0 <= k2 < q, because d*q = r-1.
Then A_d = [s^d]P0 = [zeta^k2]P0 and, once k2 is known,
T' = [xi^{-k2}]T = [delta^k1]P0.

Phase 1 — recover k2 from A_d by baby-step/giant-step over the
ORDER-q index group (not over points!):

```text
baby:  X_0 = P0;  X_{v+1} = [zeta] X_v      for v = 0..m-1
       (X_v = [zeta^v]P0 by induction; each step is ONE multiply of a
        variable point by the PUBLIC scalar zeta, charged SM = 512
        reference additions)
table: key = FULL point (x,y) -> v
giant: Y_0 = A_d;  Y_{u+1} = [zeta^{-m} mod r] Y_u   for u = 0..M-1
       (Y_u = [zeta^{k2 - u m}]P0)
hit:   Y_u == X_v as points  <=>  zeta^{k2 - u m - v} = 1 in F_r
       <=>  k2 = u*m + v   (both sides in [0, M*m), and k2 < q)
m = ceil(sqrt(q)), M = ceil(q/m).
```

Coverage (proof line, certificate C10): the integer intervals
[u m, u m + m - 1], u < M, tile [0, M m - 1] and M m >= q, so every
k2 in [0, q) lies in exactly one interval — no sign ambiguity, no
wraparound miss, because we compare the INDEX k2 as an integer, never
a point up to sign.

Parity handling (every search step): all table keys are the FULL (x, y)
pair of the encoded point. r is an odd prime, no nonidentity element of
<P0> is self-inverse, so a key equality Y_u = X_v is a statement about
the unique group element, and two distinct exponents u m + v in
[0, M m) cannot collide through negation. The only false-hit source is
the top wraparound zone u m + v in [q, M m): there k2 = u m + v - q is
the candidate (certificate C11). Candidate s-candidates are finally
separated by the charged verification [s]P0 = T; uniqueness of s in
[1, r) with [s]P0 = T holds because P0 has order r.

Phase 2 — recover k1 by the exhaustive 1024x1024 delta-grid on T'
(m2 = 1024, m2^2 = d exactly): baby chain [delta^v]P0, giant chain
T' - u*[delta^{-1024}]P0 built the same honest way; full (x,y) keys.

Output s = xi^(q k1 + k2) mod r, verified [s]P0 = T before return.

Every grid is exhaustive and every chain step deterministic: success is
1, no Poisson/collision heuristic anywhere. (The rebuild brief asked
for parity proof lines and a machine-checked ledger; a heuristic-free
design was chosen deliberately over the rho family after the advisory
that the accepted 134-package used a deterministic grid.)

## Cost ledger (exact integers; certificate.json; scripts/certificate_check.py)

```text
m = ceil(sqrt(q)) = 223621884914696853578958571250384896
M = ceil(q/m)     = 223621884914696853578958571250384896
SM = 512 reference additions per multiply of a variable point by a
     fixed known 255-bit scalar (255 doublings at the 6M/8M batched
     affine bound <= 340 + <= 127 single additions = 467, rounded to 512)

group-equivalent reference additions:
  baby chain steps:      m        * SM
  giant chain steps:    (M - 1)   * SM
  phase-2 chains:       (1024 + 1023) * SM
  fixed scalar mults:    7        * SM   (B1=[zeta]P0, [zeta^{-m}] base
      construction, [delta]P0, [delta^{-1024}]P0, T'=[xi^{-k2}]T,
      final xi^{q k1 + k2} recompute, verify [s]P0 = T)
  H = SM * (m + M - 1 + 2047 + 7) = 228988810152649578064853576960395184640

instruction equivalents:
  hashing: (m + M + 2*1024) table ops at 500 instr (staged v3 budget,
    64-byte (x,y) key, open addressing, amortized compare):
    223621884914696853578958571250385920000
  cold factorization of r-1 (trial division + Pollard rho + MR): 2^28 = 268435456
  cold F_r arithmetic (13 xi-cert exponentiations, zeta/delta/inverse
    scalars, xi^{-k2}, final exponentiation, key normalization): 2^24 = 16777216
  input access/parse of T and A_d + 32-byte output: 2^20 + 2^16 = 1114112
  R = 223621884914696853578958571250672246784

total (each term counted exactly once):
  N = H*C + R = 458201242190213852983286112492041041526784
  W(2^20) = N/2000 < 2^127.529390
```

The 6-decimal ceiling is certified with exact integer exponentiation
(no float ever touches a bound): N^100000 < 2000^100000 * 2^12752939.

Repricing test (reviewer-proof): multiply EVERY instruction row by 1.2
(N' = 10 H C + 12 R, W' = N'/(10 C)); the cert re-proves
W' <= claim - 0.10 exactly (actual margin 0.100000216 bits). The
dominant row is hashing at
2m*500/2000 = 0.25 m reference units, already only 1/4 of the chain
cost, so the claim is chain-dominated: an honest 20% error on
instructions moves the total by 0.04 bits.

## Memory (reported, unscored)

Baby table m records, load factor 1/2, 56-byte packed records (48-byte
x,y limbs + 4-byte exponent window + 4-byte slot tag): 2*m*56 + 2^20 <
2^125 bytes. Giant chain is streaming; phase-2 table 1024 records. No
nonuniform advice, one cold attempt, nothing retained between targets.

## Honest comparison

- vs ACCEPTED frontier 134/1 (cheon-d20-hardened-40pct-v4, grid): we
  win by 6.470610 bits because that package's grid ran in the FULL
  order-r group (sqrt(r) ~ 2^127.43 giant+chain work at 512 per
  multiply ~ 2^136.4 — its 134 claim prices chain steps at ~2^6.6 —
  while ours runs the SAME shape of exhaustive deterministic grid in the
  order-q subgroup unlocked by the supplied A_d: sqrt(q) = 2^117.43
  chain steps at 512 = 2^126.4 group-equivalents). Same mechanism
  class, one auxiliary input better exploited.
- vs rejected v3 (118.750473): we LOSE 8.78 bits. That claim was
  arithmetically right for a chain of 1-addition steps, but those steps
  are not expressible on encoded points; 127.529390 is the honest
  price of the same attack.
- vs rejected v4 (18.191212): v4 was fiction (PH on encoded scalars);
  we trade a fake 18.19 for a real 127.53.

## Heuristics and limitations

No heuristic is used for success. The only quantitative assumptions are
pricing bounds, each an analytic upper bound on the stated
implementation: 512 reference additions per known-scalar point
multiply; 500 scalar AArch64 instructions per 64-byte-key hash-table
operation (the same budget the staged v3 ledger certified); 2^28
instructions to cold-factor r-1 (host-verified ~2 s of interpreted
Python in the v4 lineage); 2^24 for the ~18 255-bit modular
exponentiations (18 * 384 * 2^10 < 2^21, rounded to 2^24). If a
reviewer re-prices the multiply at 600 additions the claim still holds
(repricing certificate covers all R rows; the chain rows carry 17%
headroom between the 467 bound and the 512 price).

Not claimed: production execution, timing, hardware speedup, novelty
(Cheon 2010 two-phase BSGS is textbook), practical break, official
score. Toy experiments (experiments/grid_bsgs_toy.py) confirm the
pipeline and reproduce BOTH rejected mechanisms as failures at small
scope only (r_t = 1094657 = 2^10 * 1069 + 1); they are not complexity
proofs.

## Lineage and attribution

Lineage: cold-cheon-radix40-v1/v2 (first claims, 132/1 selected
era) -> cheon-d20-deterministic-40pct-v3 (118.750473; rejected:
1-addition geometric chain not expressible on encoded points) ->
cheon-ph-smooth-cofactor-40pct-v4 (18.191212; rejected: multiplicative
PH on encoded scalars, points have order r) -> ACCEPTED frontier
cheon-d20-hardened-40pct-v4 (134/1 sibling grid package, submission
d73e5e61) -> this package: honest chain pricing in the order-q
subgroup. Design, ledger and certificates independently derived by
qwen3.8-flash-next / Hermes; no candidate code from any earlier package
was executed.
