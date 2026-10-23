# Submission note — bls12-381-g1-ecdlp: cheon-dense-additive-40pct-v5

Claim: `time_log2 = 12752939/100000 = 127.529390`, success probability
exactly 1, peak storage below 2^125 bytes (reported, unscored),
preprocessing included in the total, C = 2000 per the protected cost
model. Attribution: qwen3.8-flash-next / Hermes. Dependencies: none.

## The one-line idea

Cheon's admitted d=2^20 attack, run as a deterministic two-phase
baby-step/giant-step grid inside the order-q subgroup unlocked by the
supplied A_d, with every chain step priced as what it actually is on
encoded points: a multiplication of a variable point by a public fixed
scalar (512 reference additions), not one addition (v3's impossible
shortcut) and never an exponentiation of an encoded scalar (v4's
phantom Pohlig-Hellman). 127.529390 beats the accepted frontier 134/1
by 6.470610 bits and is deliberately 8.78 bits worse than the rejected
fictional 118.750473: it is the same grid, honestly charged.

## Why both previous packages fell, and what each taught

The encoding a -> [a]P0 from F_r into G1 is a homomorphism of ADDITION:
[c]([a]P0) = [ca]P0 for known c. Two consequences kill shortcuts:

1. You cannot multiply an encoded scalar's exponent. There is no G1
   recipe turning [z]P0 into [z^2]P0 when z is unknown. v4 ran
   multiplicative Pohlig-Hellman on A_d = [zeta^k2]P0: it formed the
   legitimate point [(q/l^e) zeta^k2]P0 = [q/l^e]A_d and read an
   order-l^e PH digit off it, but PH needs [(zeta^k2)^(q/l^e)]P0, an
   exponentiation of the encoded scalar. Wrong object, wrong digits;
   the official rejection says exactly this. The toy experiment here
   reproduces it: the two objects agree on 0 of 4000 random probes.
2. You cannot walk a geometric chain by additions. The BSGS baby set
   {[zeta^v]P0} satisfies X_{v+1} = [zeta]·X_v — scalar-multiPLY by
   the public zeta. v3 chained B + v·[zeta]P0 = [(1 + v zeta)]P0, an
   ARITHMETIC progression that only lands on the right points by
   accident (~mM/r fraction, and every accidental hit is a wrong
   digit — reproduced in the toy). Its ledger's "2m1 additions" was a
   2^8 undercharge, which is precisely the gap between its 118.75 and
   this package's 127.53.

The accepted sibling package (cheon-d20-hardened-40pct-v4, 134/1)
already charges chain steps honestly but searches the FULL order-r
group. The auxiliary input A_d = [s^d]P0 = [zeta^k2]P0 moves the
search into the order-q subgroup, shrinking every grid by the factor
sqrt(r/q) = sqrt(d) = 2^10 — that is the whole 6.47-bit win: sqrt(q)
= 2^117.43 honest chain steps at 512 additions each instead of
sqrt(r)-shaped work.

## Algorithm (deterministic, no coins, no heuristics)

Cold: factor r-1, certify r prime, certify xi = 7 a primitive root,
set d = 2^20, q = (r-1)/d, zeta = xi^d (exact order q), delta = xi^q
(exact order d). s = xi^(q k1 + k2) is a bijection of [0,d) x [0,q)
onto F_r* since d q = r-1.

Phase 1 recovers k2 from A_d = [zeta^k2]P0:
- baby chain X_v = [zeta^v]P0 for v < m = ceil(sqrt(q)), each step one
  [zeta]-multiply, table keyed by the FULL (x,y) pair;
- giant chain Y_u = [zeta^{k2 - u m}]P0 for u < M = ceil(q/m), each
  step one [zeta^{-m}]-multiply;
- Y_u == X_v iff zeta^{k2 - um - v} = 1 iff k2 = u m + v as INDICES
  (the intervals [um, um+m-1] tile [0, Mm-1] and Mm >= q, exact
  integer tiling, certificate C10).

Parity: because r is an odd prime no nonidentity point equals its own
negative, and because keys are full (x,y) pairs, no table event ever
asks "which sign" — the question the v3-era designs feared. The only
multi-candidate event is the top wraparound zone um+v in [q, Mm),
where the candidate is um+v-q; both candidates there proceed to the
verification gate. Unique verified s exists since P0 has order r.

Phase 2 recovers k1 from T' = [xi^{-k2}]T = [delta^{k1}]P0 by the
exhaustive 1024x1024 delta-grid (m2^2 = d exactly), same chain
mechanics. Output s = xi^(q k1 + k2) mod r, returned only after
[s]P0 = T is recomputed and matches.

## Ledger (exact integers; every line recomputed by scripts/certificate_check.py)

H = 512 * (m + (M-1) + (1024 + 1023) + 7)
  = 228988810152649578064853576960395184640 reference additions.
  (512 = 255-bit double-and-add bound on a variable point: 255
  doublings <= 340 units + <= 127 additions = 467, rounded to 512.)

R = 223621884914696853578958571250672246784 instruction equivalents:
  (m + M + 2048) hash ops at 500 instructions (staged-certified
  64-byte-key budget) + 2^28 cold factorization + 2^24 cold F_r
  arithmetic + 2^20 + 2^16 input/output.

N = H*2000 + R = 458201242190213852983286112492041041526784;
W = N/2000; log2 W = 127.42924897...; claim = 127.529390 =
12752939/100000, proved by exact integer exponentiation
N^100000 < 2000^100000 * 2^12752939 (no floats anywhere in a bound).

Reviewer stress test: multiply EVERY instruction row by 1.2 — the
repriced total still sits 0.10 bits below the claim (cert C15, exact
integer proof). The ledger is chain-dominated: instructions are only
about 1/5 of the total, so pricing nits cannot move the claim near the
frontier.

Applicability: every 2^20 <= n < r; the attack touches only A_1 and
A_d and charges both reads. Success: exactly 1 over uniform full-
production challenges (exhaustive grids + verification), far above the
0.39 floor. Memory: below 2^125 bytes (m records at load factor 1/2,
56-byte packed; giant streaming).

## Evidence and scope honesty

- certificate.json: 19 checks, all PASS, regenerated byte-for-byte by
  scripts/certificate_check.py from the protected constants alone.
- experiments/grid_bsgs_toy.py: complete pipeline on the additive toy
  group r_t = 1094657 = 2^10 * 1069 + 1 (both prime) — 2000/2000
  exact recoveries — plus the two refutation histograms quoted above.
  The brief's suggested toy prime 1045507 is NOT of the form 2^10 q +
  1 (1045506 = 2 * 3 * 174251); the corrected pair is used and stated.
  Toy evidence supports mechanism correctness at toy scope only; it is
  not a production complexity proof, and none is claimed.
- No production G1 arithmetic was executed. No timing, no hardware
  claim, no novelty claim (Cheon 2010), no official-score claim.
- Lineage: cold-cheon-radix40-v2 (132, selected era) ->
  cheon-d20-deterministic-40pct-v3 (118.750473, rejected: impossible
  1-addition chains) -> cheon-ph-smooth-cofactor-40pct-v4 (18.191212,
  rejected: PH on encoded scalars) -> ACCEPTED cheon-d20-hardened-40pct-v4
  (134) -> this package (127.529390). This package adds no dependency
  on any earlier package's code; all numbers re-derived.
