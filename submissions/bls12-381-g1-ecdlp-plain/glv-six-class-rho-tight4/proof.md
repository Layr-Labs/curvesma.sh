# Six-class parallel Pollard rho on BLS12-381 G1 (negation + GLV)

**Claim.** time_log2 = **126.302** = 126302/1000, the exact three-decimal ceiling
of log2 W (section 5). x is recovered from T = [x]P0 with probability
>= 0.3901659 per uniform instance and fresh coins; total work is
<= W_max = 2^126.302 units on every run (1 unit = one batched affine G1 addition,
C = 2,000 scalar AArch64 instructions), enforced by a global meter charging
an upper-bound price per action. Peak memory < 2^72 bytes (2^71.8074);
preprocessing <= 2^59 units, inside W; no advice. Success rests on
H1 and H2 (section 4). Restated from public 1c600df4 (six-class rho,
126.36); ledger, rules, store bounds, floor and failure budget re-derived here.

## 1. Algebra

p = 0x1a0111ea...ffffaaab, r = 0x73eda753...00000001, E: y^2 = x^3 + 4, and
N = (r-1)/6 classes (log2 N = 252.2721269). beta is a primitive cube root of 1 in
Fp, phi(x,y) = (beta x, y) = [lambda] with lambda^2 + lambda + 1 = 0 mod r. The
units u_t = (-1)^s lambda^j, t = 3s + j, act by [u_t](x,y) = (beta^j x, (-1)^s y);
no non-identity G1 point is fixed by a unit, so G1\{O} is exactly N orbits of 6.

**Key.** K = |y| = min(y, p-y) is a complete class invariant. On six 64-bit limbs,
y in [1, p-1]:

```
d  = p - y                        # 6 subs, borrow propagated
br = borrow(y - ceil(p/2))        # 6 sbcs; br = 1 iff y < ceil(p/2)
K  = br ? y : d                   # 6 csel
```

p is odd, so y >= ceil(p/2) iff p-y < y: exact for every y, no tie fallback, 18
instructions. Bits 0..15 of K give h = K mod R; bits 64..127 are the distinguished test.

**Canonical member.** c = (beta^j x, |y|) with j = argmin(x, beta x, beta^2 x) is
unique and equals [u_t]P. Only t is needed, from bits 384..447 of
x*floor(2^448 beta/p): 12 limb products give L', true bits in [L', L'+12]. With
PH = floor(p/2^256) = PH1*2^64 + PH0,

```
V1 = ((PH0*L' >> 64) + PH1*L') >> 64
```

is exact in five instructions:

```
umulh t0, PH0, L'
mul   t1, PH1, L'
umulh t2, PH1, L'
adds  t3, t0, t1       # carry c out
adc   V1, t2, xzr      # V1 = t2 + c
```

PH1*L' = t2*2^64 + t1, so T = t0 + PH1*L' = (t2 + c)*2^64 + t3 with c = 1 iff
t0 + t1 >= 2^64, and T >> 64 = t2 + c: the carry out of `adds` enters through
`adc`, which also lands the final shift, so no `lsr` is needed. Both products are
full (mul and umulh each), and the high half of PH0*L' meets the low half of PH1*L'
with its carry: the multiword term an earlier 4-instruction row omitted. The
truncation is exact, as T <= (2^64-1) + (2^64-1)^2 < 2^128. V1 is within 3 of
(beta x mod p)/2^320 and V2 within 5 of (beta^2 x mod p)/2^320; if one of X (top
limb of x), V1, V2 is below both others by more than 16 it is the exact argmin, else
the exact rule runs as a metered 0.25-unit action.

## 2. Walk

R = 2^16; theta = 2^-64; M = 2^20 walks; L_max = 20 theta^-1; checkpoint 64.
Table: c_h, d_h uniform from fresh coins, M_h = [c_h]P0 + [d_h]T (redraw if O),
Tab[h][t] = [u_t^{-1}]M_h with coefficients u_t^{-1}(c_h,d_h) mod r. A step is
F(C) = class(c(C) + M_h), h = K mod R, done by P <- P + Tab[h][t], which is F
since c + M_h = [u_t](P + Tab[h][t]). Trails start at [a]P0 + [b]T with fresh a, b
stored in full; starts are hard-capped at 2^63.25. At each arrival: report and end
if K is distinguished; abandon at length >= L_max; ESC on pending_escape, on
fp(K) = K_{-2} while fp(K) < K_{-1}, or on fp(K) = K_cp after a <= 64-step re-walk
from the cycle's minimum-fp member; at length mod 64 = 0 set K_cp <- fp(K); then
step, shifting the history K_{-2} <- K_{-1} <- fp(K) (fp = limbs 0, 2, 3). ESC is
P <- [2]P, a class map because doubling commutes with every unit.

**Merging lemma.** Two trails meeting one class reunite: an escape needs a repeat,
i.e. a cycle Z of F of length <= 64 with no distinguished class, which is left for
E(Z) = class(2 c(min Z)), a function of Z alone, and the other trail entering Z is
caught within 3 arrivals (2-cycle) or 128+64 (checkpoint) and escapes to the same
E(Z), within <= 192 operations. Escapes close no structural cycle: 2 has order
2^250.8571 mod r, so 1 - u 2^m != 0 mod r for m < 2^248. A key match re-walks both
trails from their stored (a,b), and x = (e_2-e_1)/(f_1-f_2) mod r is verified by
[x]P0 = T.

## 3. Meter and report store

Every action is charged before it runs at an upper-bound price, and a charge past
W_max halts with failure, so W <= 2^126.302 always. The meter is a 192-bit counter;
2^20 processors reserve 2^40-unit blocks atomically, remainders counting as spent
(< 2^61 units). **Phase 1** ends when a global operation counter (steps, doublings,
traversals) reaches S_ops, so S_ops <= operations < S_ops + 2^60; **phase 2**
finishes open trails.

**Report store (hard bound).** Records (K,a,b) are 112 bytes, 14 limbs; 2^65 slots,
double hashing h = H1(K), g = H2(K)|1: an insertion probes (h + i g) mod 2^65 for
i = 0..95, writes into the first free slot, and the attack
halts if all 96 are occupied. So an insertion costs at most 96 probes at 0.15 units each
(300 instructions: 128-bit address add, slot load, compare, covering a TLB and DRAM
access) and the run at most D_cap * 96 * 0.15 = 2^67.09 units over
D_cap = 2^63.25 insertions, a hard extra. The table cannot fill (D_cap < 2^65). At
load alpha = 2^-1.75 an insertion sees 96 occupied slots with probability
<= alpha^96 = 2^-168 under H2, union over D_cap insertions < 2^-104.7; section 4
claims only 2^-40. The 2^71.81-byte table needs 128-bit byte addresses, a two-word
add per probe; page faults plus zeroing over 2^59.8074 pages cost 1 unit each as a
hard extra, below the phase-2 term, and the write pass is the 2^57.84
preprocessing row under the 2^59 sub-cap.

## 4. Success probability

N' = N R/(R-1) = 2^252.2721489. T_c counts class arrivals of all trails up to the
first arrival at a class already visited by another trail.

**H1 (heuristic, at the budget).** P(T_c <= S_ops) >=
1 - exp(-(S_ops/1.05)^2/(2N')), the 5% allowance against the random-mapping law
(toy walks -0.2% to +2.1%). Also: |y|'s bits are pseudo-random in the class;
apart from detected cycles of length <= 64 (2-cycles at rate 1/(6R) per step,
longer ones < 2^-30) a trail's non-fruitless return to an earlier class has
probability <= L_max^2/(2N); and a new class's twist index is uniform on the 6
units, independent of its partition index.

**H2 (random oracles).** H1 and H2 of the store are independent random oracles into
2^65. DP indicators of distinct classes are independent because the DP test reads
bits 64..127 of K and h reads bits 0..15, disjoint ranges.

S0 = sqrt(2N' ln(1/0.6094)) = 2^126.1292332, S_ops = 1.05 S0 = 2^126.1996225,
U1 = 1.071 (S_ops + 2^60) = 2^126.2985810; the pure model gives exactly 0.3906, and
phase 1 performs >= S_ops class arrivals (hard counter), so T_c <= S_ops is reached
inside phase 1. Surcharges: 0.429 per doubling at rate 1/(6R) + 4R^-2 = 2^-18.585 +
2^-30 per step (a 2-cycle needs the next class to share the partition index with
canonical twist -1), and 0.012 per traversal step at rate < 2^-30 with <= 64 steps
and <= 16 escapes per trail: by linearity, expected 2^106.394238 and 2^99.818801.
Random extras expect < 2^73.9: twist fallbacks < 2^72.3, and words and sampling
attempts < 2^63 each.

**Failure events** (union bound; if T_c <= S_ops the lemma gives both trails one
report, so the attack succeeds unless one occurs): two-trail cap tail
2*21e^-20 = 8.66e-8; trail-start cap exceeded < 2^-40 (mean 2^62.2 starts, cap
2^63.25); store overflow or table full < 2^-40 (analytic < 2^-104.7);
DP-indicator dependence < 2^-40; a non-fruitless return 2 L_max^2/(2N) = 2^-115.6;
f_1 = f_2, 1/r; exceptional additions < 2^-114; 8 matches failing < 2^-60; a
spurious fingerprint match < 2^-64; preprocessing cap < 2^-38. All but the early
halt total < 1.0e-7 (the three leading terms sum to 8.7e-8).

**Success.** Headroom above the deterministic charges at W_max = 2^126.302 is
2^117.579312 against expected random charges 2^106.409287, so by Markov the meter
halts early with probability <= 2^106.409287/2^117.579312 = 0.000433998. Hence
success >= 0.3906 - 0.000433998 - 0.0000001 = 0.3901659 > 0.39, margin 1.6590e-4.

## 5. Cost

Per ordinary step, besides the reference addition, the ledger row by row: loads
12 (point 6 ldp, addend 6 ldp); key 18 (6 subs, 6 sbcs vs ceil(p/2), 6 csel);
L' 24 (12 products: 6 umulh, 6 mul; 12 accumulates); L' range 4 (2 cmp, 2 b);
V1 5 (umulh, mul, umulh, adds, adc); X, S and the wrap test 6 (add, sub, cmp, abs,
b); V2 4 (sub, sub, cmp, csel); argmin, two gap tests, fallback 13 (unrolled 3x);
h, t, table address 5; distinguished test on limb 1 of K 3 (cmp, b);
exceptional-addition test 7; history 26 (2 ccmp-chain compares, flag, checkpoint,
6 movs); point store and loop 12 (6 stp, length++, cap, meter add, branch).
Subtotal 139, slack 3 per step, so a step is 142 instructions.

The reference 6M covers the field arithmetic and the batch inversion; the charged
loads and store are slack. Every row is at or above its own instruction
sequence, so the 3 of slack is spare margin and a single-row audit cannot cross it.

**Prices** (units): 1.071 step, 1.083 traversal, 1.5 doubling, 0.25 twist
fallback, 0.005 word, 0.025 sampling attempt, 2^12 trail start, 0.5 report, 0.15
hash probe (<= 96 per insertion), 2^28 table, 142 phase-2 or re-walk operation,
132 small-batch inversion, 2^9 solving, 2^6 output, 2^13 verification, and the
2^57.84 report-table init of the memory paragraph.

M = C/6 = 333.33, so a doubling is 7M + 128 = 2,461.3 instructions, 2,464.3 with
the same 3 of slack, below 3,000 = 1.5 units (asserted in `experiment.py`). The 128
is the doubling's share of the Montgomery batch: its 2y denominator is queued in
the same inversion batch as the additions, so a doubling is never a small batch and
doubling-only batches number 0 (< 2^8); the only partial flushes are at most one
per processor at the run's end, <= 2^20 of them, charged 2^28 units.

**Budget.** Phase 1: 1.071 (S_ops + 2^60) = 2^126.2985810, plus surcharges
(expected 2^106.409287). Hard extras: meter remainders < 2^61; phase 2
2^20 (L_max+65) 142 = 2^95.4717; re-walks 2^79.4717; reports D_cap 0.5 = 2^62.25;
probes D_cap 96 * 0.15 = 2^67.09; starts 2^63.25 * 2^12 = 2^75.25; partial flushes
2^28; solving 2^9; output 2^6; match-lookup probes 8*96*0.15 = 116; page faults
2^59.8074; preprocessing 2^59; verification 2^13: together < 2^95.471698. So
W = 2^126.2985810 + 2^95.471698 + 2^106.409287/0.000434 = 2^126.3016727
<= W_up = 104855684481168075414614794644037112765, log2 W_up = 126.3016726813, and
126.302 is its exact three-decimal ceiling: 2^126302 >= W_up^1000 and
2^126301 < W_up^1000, checked with integers in `experiment.py`.

**Memory.** 2^65 report slots of 112 bytes = 2^71.807355, table 2^25.91, walk states
2^29: total 4,132,070,672,511,539,347,456 = 2^71.8074 < 2^72. **Preprocessing** (inside W):
report-table init 2^57.84 units, table 2^28, first starts 2^32, sampling 2^21,
expected < 2^57.85 under the 2^59 sub-cap; the random part gives a halt
probability < 2^-38 (Markov).

## 6. Experiments

`experiment.py` (stdlib) checks Part A's identities of sections 1-2 and all three
production rules at the six-limb width on 100,000 real production trail positions:
key rule and twist index 0 mismatches (0 fallbacks), the distinguished rule
on limb 1 of K 0 mismatches and 0 hits against its analytic 2^-64, and 256
constructed keys over that bit-range test. Part B walks a toy
curve over F_q with the scaled pipeline (sequential trails, so phase 2 is vacuous),
scales the rules to its one-limb field (the production bit range is the Part A
check) and compares them with the exact rules over 1,000,008 toy arrivals, and runs a two-hash store at load 1/2 with the 96-probe cap.

The declared run (q26, R = 64, 4 distinguished bits, 12,000 trials): pipeline
success **0.4198** (stderr 0.0045) and fraction with T_c <= S_ops 0.4146, both
above the 0.3906 floor; 2-cycle rate 0.002434/step against 1/(6R) = 0.002604; 387
checkpoint escapes, 5 abandoned trails, 3 exceptional additions; rules 1 fallback,
0 mismatches; store 8,192 insertions, at most 9 probes of 96, no overflow. Seeds 2
and 3 give 0.4125 and 0.4203: 0.4175 pooled over 36,000 trials, 10.5 sigma above
the floor. 
