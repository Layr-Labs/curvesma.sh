# Argument: disjoint DP/routing fields, fixed-cap refusal, bounded detection (v9c)

Claim: total classical work <= 12939/100 = 2^129.39 target-group-operation
equivalents at C=1000, success >= 0.39, advice 0, peak memory < 2^73 B.
Live head at dispatch: 129.46 (`5fa78d62`). Supersedes v9b (`5ed5529`, rejected).

## 1. Target, normalization, algorithm

NIST P-256, p = 2^256-2^224+2^192+2^96-1, a = -3, cofactor 1, order r prime. C=1000
is the borrowed secp256k1 unit; P-256 arithmetic is charged by its own schedule
(section 2). H=0, so score = log2(T/1000) with T in instructions.

Construction inherited unchanged from v9b: negation classes {P,-P}, N = (r-1)/2
classes, key = x, canonical member (x,|y|), transition canon(2X_i+T_j) with
j = H(X_i) mod 2^16 and T_j from a +/- twisted table, an escape that doubles the
minimum-fingerprint class of the detected cycle (so merged walks stay merged), and
recovery x = (a'-a)/(b-b') mod r at a class collision, verified by [x]P = Q.

## 2. Complete step price (3023 -> 3087 instructions)

canon(2X+T_j) = X + (X+T_j): two batched affine additions, each 2M+1S+6A, its
denominator inverted in a batch slot at 3M. Rows are instructions, with loads,
stores, addressing and control counted:

  two affine additions, core     4M+2S+12A                     1484
  batch-inversion products       2 slots x 3M                  1320
  batched Fermat inversion share 2 x ceil(51856/2^16)             2
  coefficient update             (2a+k, 2b) mod r                63
  canonicalisation of y          p-y, cmp, csel, store           14
  degeneracy tests               x2==x1 for both adds            12
  bounded cycle detection        2 fp compares + ring + 1/64
                                 checkpoint min-select            16
  coefficient negation           2 x (negate 16 + select/store 8) 48
  control (itemised, 6 rows)     hash/DP/route/cap/load          128
  total                                                         3087

M=220, S=176, A=21 are itemised in certificate.json. Detection and coefficient
negation are new rated rows: v9b charged only the escape doubling and its
14-instruction canonicalisation row covered y alone. Inversions are batched across
>= 2^15 concurrent trails (2^16 slots), not across steps; the Fermat row is
256S x 176 + 30M x 220 + 200 = 51856 per inversion.

## 3. Ledger (P4: one constant, imported)

r-1 factorises completely (certificate.json); ord_r(2) = N.

KAPPA = 11/10 lives in `experiments/protocol.py`, imported by the toy and this
ledger: S_ops = ceil(KAPPA sqrt(2 N' ln(1/0.60))) = 2^127.652966, budget ratio
1.111844. v9b's toy used 1.06128 (the 1.05 allowance) while its formula used 1.10.

  phase 1      = 3.087 (S_ops + 2^60)           = 2^129.279171
  reports      = 2^63.652966 x 128 instructions = 2^60.687
  setup        = 2^16 table entries + node init = 2^25.045
  total                                        = 2^129.279171
  stress (every row x 6/5)                     = 2^129.542206
  claim                                        = 129.39

Escapes are arrivals: a doubling of at most 1620 instructions consumes a step slot
priced at 3087, so every escape is charged inside S_ops, with no escape rate
assumed.

## 4. Reports, capacity, routing (P1)

DP width 64, so expected reports = 2^63.652966 <= 2^64 = cap/4. Capacity is 2^16
nodes x 2^50 records = 2^66 records = 2^72 B < 2^73 B; cap/mean per node is
5.087772. Chernoff: log2 P[reports >= cap] <= -8.8e19, and node-wise with a 2^16
union bound the same bound is below 2^-100.

Field spec (replaces v9b's). The key is the full 256-bit class key
K = x. The DP predicate at width b tests mix64((K mod 2^64) ^ DPSALT) mod 2^b ==
0: it reads ONLY bits 0..63 of K. The router computes
rw = mix64(((K >> 128) mod 2^64) ^ RSALT), node = rw >> 48: it reads ONLY bits
128..191 of K, giving 16 bits -> 2^16 nodes. The fields are disjoint, so the DP
event never constrains the route; mix64 is a bijection, so at b = 64 the predicate
pins exactly one value of bits 0..63 and leaves bits 128..191 free. v9b fed the
SAME low 64-bit field to both functions, so at b = 64 its route was constant -
the collapse the reviewer identified.

Evidence at production width (phase G): 150000 DP-conditional 256-bit keys through
the real predicate and router reach 58884 of 65536 nodes, max load 11 against
expected 2.2888, union bound 2^-30.74 at tolerance 24; 1.5M keys attain all 65536
route values; the v9b spec reaches 1 node with load 150000 on the same keys.

## 5. Fixed-cap refusal (P2)

Each node cap is fixed BEFORE the run from the DECLARED expected report count,
never the observed one: cap = max(1, ceil(5.087772 x declared reports / nodes)).
A report finding its node full is REFUSED, counted, and flags the trial cap-loss;
it is not stored, matched or overwritten, and the trail walks on. Every pipeline
cell runs that path with 0 cap violations; the same seeds replayed without refusal
show 74 cap violations in the q26 dp6 cell, exactly what the protocol prevents.
v9b derived its cap after execution from the observed total, so its "0 overflows"
tested nothing.

## 6. Success arithmetic with named debits

success >= 0.39 - P_halt - P_cap - P_refuse, with
  P_cap    <= 2^-8.8e19                    (Chernoff, section 4)
  P_refuse <= P_cap                        (same event: a full node)
  P_halt   <= extras/headroom = 2^-64.95   (Markov; extras = report row + setup,
            headroom = 2^129.39 - 2^129.279171)
No successful trial is truncated: trial cost is a hard function of the step
counter, 3.087 x S_ops + extras = 2^129.279171 <= 2^129.39. Modelled success at
S_ops is 1-0.60^1.21 = 0.4610; the toy measures walk efficiency 0.864-0.919
(kappa 1.043-1.076) in the refusal-on production cells, giving 0.42-0.44; claim
0.39.

## 7. Bounded detector and memory (P3)

Detection is O(1) per trail: a 2-slot ring compares fp(K) with fp(K-2) and one
checkpoint register carries the running minimum fingerprint of the previous 64-step
block. No visited set survives: no history scan or minimum-select is charged.
Memory: 2^16 trails x 4 entries x (8 B fingerprint + 64 B point + 64 B coefficient pair) = 2^25.09 B, against the 2^15 x 2^64 = 2^79
live records (2^85 B) carried history would require; peak is 2^72 B of records +
2^25.09 B of detector, below 2^73 B. Detection costs 16 instructions per step.

## 8. Toy (production protocol end to end)

q22/q26, 400 attacks per cell (3200 in the volume cell), fresh challenge and table
per attack; one DP endpoint per trail, endpoint match on the class key, exponent
reconstruction, exact [x]G = Q, the shipped detector, router and fixed-cap refusal.
Refusal ON / same seeds OFF: q26 dp6 0.4275 / 0.4350 against model 0.4626, 144 of
8209 reports refused in 84 of 400 trials; dp4 0.4750 / 0.4700; wide table R=1024
0.4150 / 0.4350; q22 dp4 0.4275 / 0.4275; double budget 0.8900 (model 0.9166). At
8x volume the fixed cap rises to 46 (5.09x the declared mean) and refusals fall to
189 of 65928 (0.29%) with recovery identical to the no-refusal replay (0.4153).
No verification mismatch occurred in any cell. Adverse results are kept: dp=9
recovers 0.0000 with the fixed cap (floored at 1 at toy volume, 519 of 1587
refused) against 0.2325 without refusal, and the declared cap-stress cells (cap =
half the mean; a 16-node saturation control) recover 0.0000/0.1725 and 0.2650 with
0 violations - stress, not the production regime. Run: 116 s, 2.02e7 steps.

## 9. Corrections and limits

v9c repairs: (P1) disjoint DP/route fields, collapse demonstrated and removed;
(P2) fixed-cap refusal executed, refusals and post-refusal recovery reported;
(P3) bounded 2-cycle + checkpoint detector, charged detection and negation rows,
memory no longer 2^79 records; (P4) one 1.10 constant imported by toy and ledger.
Cost: step 3023 -> 3087, claim 129.36 -> 129.39, gate gap 0.07 bits. ord_r(2) = N
is even, so 2^(N/2) = -1 and the shortest sign-ambiguous gap is N/2 = 2^126.414 x S_ops. H1 carries the success bound and is
supported by the toy census, not proved. Instruction counts are analytic. Not a
production run, not a measured kernel, not a practical attack.
