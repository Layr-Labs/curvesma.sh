# Six-class secp256k1 rho, R = 2^32, κ = 1.000041, a pathwise waste bound on every catch, metered fe_sub corrections, a catch budget and a Markov allowance (125.834346)

**Claim.** For a uniform full-production secp256k1 input `Q=[x]P`, the classical cold six-class (negation + GLV) parallel Pollard-rho attack of `upstream-proof.md`, with the counted kernels and
arrival rows of this document, recovers `x` with probability greater than `0.39` under the explicit H1 heuristic, the disclosed borrow-rate premise C (section 5.4) and the disclosed catch premise P1 (section 5.5). A global meter charges every action before it runs, so on every run its total work is
at most

```text
time_log2 = 125.834346 = 62917173/500000
```

target-group-operation equivalents (C = 1000 scalar AArch64 instructions per unit). Peak memory is below `2^72` bytes; preprocessing is capped at `2^59` units, charged and included; nonuniform advice
is zero. This is an analytic heuristic bound, not a practical break.

## 1. What changes

This candidate builds on the accepted package c8a6819e-f955-403f-9543-45978f3a5dd9 (branch 9f06f2bf, official score 125.834417) and keeps its kernels, deposit ρ′ = 813/2^23, catch budget, Markov allowance, census support and dependencies. Three changes: (i) the partition table has R = 2^32 rows (2^46 hard row, 2^45.1 preprocessing, 2^42 bytes), shrinking the catch waste 1/(4R); (ii) with that pathwise waste the production collision constant becomes κ = 1.000041 with w = 2^-20.00, so p0 = 0.39001172 and the Markov headroom h = 0.0000116 keep the success floor at 0.39000001 > 0.39; (iii) P1 now bounds catches of every cycle length ≤ 32, not only 2-cycles: 5.5's census measures doublings of any cycle length, so 5.5's Ville lemma bounds the total catch count pathwise, and charging every catch the worst-case 16,158 wasted arrivals gives w < 2^-20.0193 for the realised run, not only in expectation (5.1, 5.2). The closest walk run is 1.000 s.e. (5.2). Kernel counts and CORR are unchanged; only the BFI slot width (33 bits) changes, and the ledger and dependent bounds are recomputed in `certificate.py`.

## 2. The counted kernels

The field schedules of the loops (`kernels/secp_loops.s`, section 2.5) come from the emitters of 2d25756c's `gen_loops.py` (not shipped), with the x3 and y3' forms of section 2.2; the experiment
executes and counts the shipped assembly. p = 2^256 − C0, C0 = 2^32 + 977, so 2^256 ≡ C0 (mod p). A field element is four little-endian 64-bit limbs; **weak** means [0, 2^256), **canonical** [0, p). Walk and table
points are canonical, except after the rare event of section 2.4. x4 holds C0 in the loops.

Sections 2.1 (fe_mul: 32 products and six carry chains, 62), 2.2 (the reduction, 23, its fused forms and their rare events) and 2.3 (fe_sqr, 64)
are in `kernels.md`, numbered as here.

### 2.4 fe_sub (9)

**fe_sub(a, b)**, a weak, b canonical: SUBS, 3 SBCS; CSEL of C0 on borrow; SUBS, 2 SBCS, SBC. On borrow the wrapped value is at least C0 + 1, so subtracting C0 cannot borrow again; the result is weak,
≡ a − b. Without a borrow the result is a − b, which is at least p only if a ≥ p and b < C0. The step uses it for dx = x_T − x_W (loop A, or the unfused pass 1), dy = y_T − y_W (prod) and u = x3 − x_W
(pass 2); b is a canonical walk coordinate and a is canonical (x3 outside the events above), so the results are canonical. In the loops its last three instructions (SBCS, SBCS, SBC) are a rare
exit (section 2.5).

**Branch-only form** (dx in loop A, dy in prod, u in pass 2; pass 1 and the doubling keep the form above): SUBS, 3 SBCS, B.HS past the correction;
on a borrow (a < b as integers), ADD of 1 to the loop's count register, SUBS of C0 from limb 0 and B.LO to the same rare exit. The CSEL form subtracted 0
without a borrow, so both compute the same value: 5 instructions on the main path, 3 more per borrow, paid from the account A (section 4).

### 2.5 The loops (37.5078125, 169.50390625 and 337.50390625 per walk)

**The batch.** Slot i (64 bytes, i = 0 … n − 1, n = 2^23) holds dx_i at +0 and e_i at +32. Each round runs, in order: loop A (`kernels/secp_fuse.s`: each walk's arrival, then dx_i into its
slot), `prod` (forward: dy_i from the table point and the walk's record, its point at x1, the prefix products and e_i), the zero test and inversion of section 4, and `pass2` (backward: λ_i and the new point).
**Invariant:** no walk's point record is written between its arrival and prod. Inside the arrival the rare path may replace the walk's point (a report's next start, also in an anchor round's arrive.s),
writing the record before the walk's dx; otherwise only pass 2 (x3, y3') and, after it, a doubling's move-back (section 4) and any restart after an abandonment write records; pass1, prod and the
zero-test re-run write none. So prod's y(W) belongs to the point of loop A's dx. With acc_(−1) = 1,

```text
prod:   acc_i = acc_(i-1)*dx_i,   e_i = acc_(i-1)*dy_i                       (i = 0 .. n-1)
pass2:  R = 1/acc_(n-1);  lambda_i = R*e_i,  R <- R*dx_i                     (i = n-1 .. 0)
```

Before λ_i, R = 1/acc_i, so λ_i = dy_i·acc_(i−1)/acc_i = dy_i/dx_i, and afterwards R = 1/acc_(i−1). Walk 0 is an ordinary slot: prod starts from acc = 1 (MOV, three MOV of xzr), so e_0 = dy_0
and λ_0 = dy_0/dx_0; pass 2's last R ← R·dx_0 (= 1) is computed and not used. That is Montgomery's 3 products per element plus λ, and every product shares one operand with another one in its pass:
acc_(i−1) in prod, R in pass 2.

**Shared rows.** For a weak b, B_0 = b and B_(k+1) = [lo, b0 + hi, b1, b2] with (lo, hi) = b3·C0 (MUL, UMULH, ADDS, B.CS: 4), so B_(k+1) ≡ B_k·2^64 because 2^256 ≡ C0. For a weak a, Z = Σ a_k·B_k
≡ a·b: row 0 is 8 products and 4 additions into five limbs, rows 1–3 each 8 products, 4 additions and 6 more (ADDS, 4 ADCS, ADC) into Z, interleaved for the two products on the same rows; no carry
is lost (a high word is at most 2^64 − 2) and Z < 2^322, so z5 ≤ 3. The fold (6) is MUL, UMULH of z4·C0, MADD hi' = z5·C0 + hi ≤ 4C0 − 1, ADDS z0 + lo, ADCS z1 + hi', B.CS: t = (Z mod 2^256) +
(z4 + 2^64·z5)·C0 ≡ Z. Each product is 32 + 34 + 6 = 72, and a pair on one b's rows is 2·72 + 3·4 = 156. The macro W2P does this for b in x25–x28 and the slot at x3: it loads the first operand
(+32) and the second (+0) two limbs at a time (4 LDP), leaves the first product in x10–x13 and the second in x25–x28 as the next b, and uses all 25 pool registers, never x0 or x18. prod's
W2Q is W2P with its first operand dy already in x5, x6, x29, x30 (its two LDP dropped), dx read forward and four registers renamed.

`kernels/secp_loops.s` holds each walk's body once as macros (W1, W2Q, W2P, W2 and their schedules, with register parameters); nested macros list the 128 copies of the fused pass and the 256 copies of each product loop.
They expand to straight-line code apart from the correction skips, counted on the main path: no borrow and no rare exit (an exit not taken costs its branch):

| loop | walk | 128 walks | per walk | tail walk | fixed |
| --- | --- | ---: | ---: | ---: | ---: |
| `pass1` (unfused; anchor rounds and the zero-test re-run) | 2 LDP x(T), 2 LDP x(W), dx 7, 2 STP: 13 | 128·13 + 64 + 1 = 1,729 | 13.5078125 | 16 | 8 |
| `prod` | 2 LDP y(T), 2 LDP y(W), dy 5; W2Q on acc (2 LDP of dx_i; 3 rows; e_i and acc_i: 2 + 12 + 144), STP of e_i (2): 169 | 256·169 + 128 + 1 = 43,393 | **169.50390625** | 172 | 28 |
| `pass2` | below: 337 | 256·337 + 128 + 1 = 86,401 | **337.50390625** | 340 | 28 |

pass1, prod and pass2 make one pointer LDP per pair, and one branch closes a group. The pointer array tp ends where bit 26 is first set (&tp[n]; 8n = 2^26 bytes), with zero words
tp[−2..−1] and tp[n..n+1], set once with tp; pass1 closes with TBZ x2, #26. prod's and pass 2's groups exceed TBZ's ±32 KiB, so they close with CBNZ on a preloaded table pointer, falling
through at a zero word (no table pointer is 0). prod's pair LDP comes last and loads the next pair's pointers (x5, x0), CBNZ x5. Pass 2 (backward) runs W2P, W2 on x0, W2P, then one LDP of
the second walk's pointer into x5 and the next pair's first into x0, then W2 on x5; CBNZ x0 closes the group, so x18 is free. An empty group loop is skipped by TBNZ x2, #26 or by the first
LDP (pass 2: LDR) and CBZ. prod ends by storing acc_(n−1) at slot n + 32 for the
inversion; its walk pointer x1 starts at the first y (ADD #32, in its fixed work). **Pass 2 walk** (backward). R stays in x25–x28, loaded once from Rinv, never stored. Per walk:

| work | instructions |
| --- | ---: |
| W2P: 2 LDP each of e_i, dx_i; R's rows B_1..B_3; λ = R·e_i, R ← R·dx_i | 4 + 3·4 + 2·72 |
| 2 LDP y(W_i) (kept to the end) | 2 |
| x3 = λ^2 − x_W − x̃_T + 2^257, fused; 2 LDP x(W_i) (kept), x̃_T's 2 LDP inside it | 73 + 2 + 2 |
| u = fe_sub(x3, x_W); 2 STP x3 | 5 + 2 |
| y3' = λ·u + y_W, fused; 2 STP y3' | 89 + 2 |
| **total**: 2·72 + 3·4 + 73 + 89 + 5 + 10 loads + 4 stores | **337** |

The tail loop's LDR x5 follows W2P. With λ = (y_T − y_W)/(x_T −
x_W) and x3 = λ^2 − x_W − x_T, the stored point (x3, λ(x3 − x_W) + y_W) is −(W + T), canonical outside the events of section 2.2.

**W2, renamed.** The x3, u and y3' code (macro W2) is 2d25756c's schedule with λ in (x22, x21, x9, x16); W2P leaves λ in (x10, x11, x12, x13), so from the off-diagonal products on, W2 and its three stubs exchange x22 ↔ x10, x21 ↔ x11, x9 ↔ x12, x16 ↔ x13. λ is the only renamed input, so the renamed code computes the same values; `evidence.md` items 2 and 5 run it.

**Loop A** (`kernels/secp_fuse.s`, two copies: `pass1f` in compare rounds, t ≡ 0 (mod 256), `pass1p` in the others). Each walk runs the arrival of section 3.2 on its current point, then: 2 LDP of x(T),
dx = fe_sub(x(T), x(W)) on the arrival's x(W) in x20–x23 (5), 2 STP of dx: 9.
- **pass1f.** LDR of the K_cp word, 3 LDP of x(W) and y0, y1, the F1 row and index (21), TBNZ, the key (2: CMP y1, #0 and ORR, or CMN y1, #2 and SUB from p0 + 1, section 3.2), the chain (5), B.NE, BFI and STR of the address (32), then the 9: **41**. TBZ x2, #26 closes 128 walks: 128·41 + 1 = **5,249 per 128 walks, 41.0078125 per walk**; a tail walk is 45 (its count kept on the stack).
- **pass1p.** No K_cp load, filter word or filter CCMP; the chain ends with B.HS: after CCMP d, C is set exactly when there is no distinguished candidate or window. Even y tests K1 = 0 by CMP y1, #0; odd y by CMN y1, #2: y1 = ~K1 − b0 (b0 the borrow of p0 − y0), so K1 = 0 gives y1 ≥ 2^64 − 2 and sets C, which sends the arrival to the rare path; it rejects the false candidates (K1 ≠ 0, below 2^-63 per arrival). The arrival is 28 (3 LDP, F1 11, select 6, TBNZ, one key test, 4 windows, B.HS, BFI); a pair's two addresses go out in one STP (the first in x12, p0 + 1 being unused), and x0 counts the tail loop (SUBS, B.NE); TBZ x2, #26 closes a group: walks 37 and 38 per pair, **4,801 per 128 walks, 37.5078125 per walk**, tail walk 40. Its rare path has no filter part.
- Both copies: fixed work 33 (prologue 22 with MOV x17, xzr; the tail-count test 2; TBNZ x2, #26; MOV x0, x17, epilogue and RET 8). x17 is the dx borrow count (returned in x0); x5, x6 and x19 are unused. Its exit (20 milli-units) charges the caller's x19 at [sp] (RCHG: LDR x13, SUB, STR, with x13 dead there).
- **Rare path.** It is that of section 3.2, placed after RET (keeping a group in TBZ's reach) and entered by one B after the walk's B.NE (B.HS): an action of 25 milli-units taken first from the same word, then at most 22 instructions, 23 with that B (the rare-path row of section 5.3), after the walk's tests and before its LDP of x(T). It may use x13, x25–x27, x29, x30 and the stack, keeps x0–x12 (pass1p: x12 holds the pair's first address), x14–x19 and the pair's other address, and leaves x(W) in x20–x23, the walk's b_j in its address register (x24; x12 for pass1p's first walk of a pair) and the BFI source in x28, for the walk's current point (after a fallback, the exact rule's j; after a report, the next trail's start, arrived inside the start's 2^13 units, written to the walk's record). prod and pass 2 read the walk's point from its record (the invariant above), so the rare path hands on no y in registers. Reports, starts, fallbacks and traversal surcharges are charged as in section 4, never from the reserve. The shipped `secp_fuse.s` has a test stub `.Ls…` there instead: after the 25-milli-unit charge it appends (x13, x1, x28, x29) to a list at [sp, #120] and branches to the BFI, so the address is the kernel's own j; it uses only x13 and the stack.

Every walk makes exactly one arrival and one operation per round; round 0's loop A arrives the starts. When phase 1 stops, `arrive.s` runs once over the last round's new points before phase 2,
so every phase-1 point is arrived. No round reserve covers it: it is one action priced **2^18 units plus 25 milli-units per walk < 2^18.85 units**, charged before it runs, in the 2^70 "other"
row (5.3), covering its 2^23·31.015625 instructions, setup < 2^10 (< 2^18·1000) and at most one rare entry per walk at 25. Its rare entries log their kind in the rare list and draw nothing from
the meter: its rare path skips the 25-milli-unit charge, prepaid in the price (`arrive.s` has no meter decrement). `certificate.py` asserts both inequalities.

**The stored member** (4595bad7). The walk keeps any member of a class (`upstream-proof.md` section 2). For S with t = 3s + j, −S has the same x and K and y of the other parity (y ≠ 0),
so s' = 1 − s, Tab[h][t'] = −Tab[h][t] and −S + Tab[h][t'] = −(S + Tab[h][t]), in the same class. Storing −(W + T) changes no class, key, index, fingerprint
or report, only the member held; doublings, traversals, re-walks and solving act on classes or recompute their members.

**Tails.** Each loop first runs its n mod 128 left-over walks in a tail loop (the table's tail walks; an empty one is skipped by one CBZ), then the groups. n = 2^23 leaves none;
section 4 still charges each round for 127 per loop.

**Rare exits** (kaankolcu's a9297fb0). y3's fused form adds y_W, x3's borrow chains end with ADC and its second subtrahend is read at [p, #192] (section 2.2). The tails below are replaced by a branch
to a stub `.Lr…` after the function's RET (one copy per walk copy; B.LO and B.CS reach ±1 MB, more than the 361,240 bytes of `secp_loops.s`'s assembled code), which ends with B to `.Lb…` just after
the exit. An exit is taken exactly when its tail would change a register:

- **fe_sub** (dx, dy, u). B.LO after the SUBS of C0 (pass1: of the CSEL'd C0) from limb 0; the stub is SBCS, SBCS, SBC. It is taken iff a − b borrowed and the wrapped low limb is below C0. Otherwise (C set, or no borrow and B.HS skipped it) the tail would subtract 0.
- **x3's and y3''s folds.** B.CS after the ADDS, ADCS of (w0, w1) into limbs 0 and 1; the stub is ADCS, ADC of 0 into limbs 2 and 3. It is taken iff s1 + w1 + c ≥ 2^64, so only if s1 ≥ 2^64 − 2 (w1 ≤ 1: s4' ∈ [0, C0 + 2] for x3 with the 2p bias, [0, C0 + 1] for y3', section 2.2); there is no sign exit.
- **Rows and their folds** (W2P, in prod and pass 2). A row's B.CS exits when b0 + hi carries (b0 ≥ 2^64 − C0 + 1, at most C0/2^64), a fold's when z1 + hi' carries (z1 ≥ 2^64 − 4C0, at most 4C0/2^64). The stub ST adds the carry into limbs 2 and 3 (ADCS, ADCS, B.CC back); if that also carries, the value wrapped 2^256 and is below 2^99, so it adds C0 into limb 0 (ADDS, ADCS, ADCS, ADC, B back). e, acc, λ and R are therefore always exact weak residues, and no failure row is added.

In `secp_loops.s` each stub first runs the macro RCHG (STR x0 pre-indexed, LDR, SUB #20, STR, LDR x0 post-indexed), taking 20 milli-units from the caller's x19 (the round's reserve, section 4),
saved at [sp] by the prologue (pass1 saves only x19) and restored by the epilogue. No stub sets flags before its tail and x0 is restored, so the tail sees the flags and registers it would have seen
inline. A stub is at most 13 instructions (ST: 5 + 3 + 4 + B). There are **3,727** exit stubs in `secp_loops.s`, 129 fused-pass copies plus 256 copies in each of prod and pass 2, with 1 in pass1 (dx), 6 in
prod (dy, three rows, two folds) and 8 in pass2 (three rows, two folds, x3, u, y3'), and **258** in `secp_fuse.s` (1 per walk copy, 129 copies in each of the two functions). A walk enters at most 15 stubs per
round (1 in loop A or pass1, 6 in prod, 8 in pass2), 300 milli-units. The doubling's fe_mul, fe_sqr and x − x3 take the same forms (section 4); its fused x3 and y3 keep their inline folds.

**Registers.** prod and pass2 use the pool x5–x17, x19–x30 (25; x19–x30 saved by the prologue and restored before RET); pass1 uses x5 and x9–x17. x0 holds the tail count
(in prod's and pass 2's groups then a table pointer), x1–x3 the arrays, x4 C0, x18 prod's dy or pass 2's u borrow count, zeroed at entry and returned in x0; no stub touches x18. The loops
call nothing and have no spills; only the stubs touch the stack. Each loop is a single basic block apart from its rare exits and correction skips (B.HS over 3 instructions; `a64.py` counts the
main path with every skip taken): every other conditional is a carry chain or CSEL, and the only other branch is the back edge. One LDP or STP is one instruction.

### 2.6 Evidence

The declared run's checks (items 1 to 9), the off-runner kernel runs and the interpreter are in `evidence.md`.

## 3. Arrival rows: 31.015625 plus the round share, executed

The arrival rows are the arrival pass's work per class arrival (section 5.1), in loop A (section 2.5) or, in anchor rounds, `kernels/arrive.s` (run in the declared run): 30 instructions on either
parity path, plus per pair one LDP of the K_cp words and one STP of the addresses (a9297fb0's pairing) and per 32 pairs one loop branch, plus the round share (section 4).

| row | this | sequence |
| --- | ---: | --- |
| F1 truncated product | 13 | 2 LDP of x; 4 UMULH, 3 ADD, 4 MADD (G_1..G_5 in registers); F0 = top limb of x |
| canonical index, slot | **10** | ADD, SUBS, CSEL, CSEL, ADDS, CSEL; CCMN ×3, CCMP (section 3.1) |
| sign, key, DP test | **3** | TBNZ on bit 0 of y0; odd y: SUBS, SBCS; even y: CMP, ORR |
| h and address | **1** | BFI of bits 0..32 of y0 into the slot base (odd half permuted, 3.2) |
| traffic | **1.5** | LDP of y0, y1; half an STP of the pair's two addresses |
| filter, rare branch | **2.5** | half an LDP of the pair's K_cp limb-0 words; CCMP; B.NE past the rare path |
| arrival loop | **0.015625** | one TBZ per 32 pairs on the pointer array's end bit |
| **total** | **31.015625** | |

### 3.1 The canonical index

As in `upstream-proof.md` section 1.1, v0 = x, v1 = βx mod p and v2 = β²x mod p are distinct in [1, p) with sum p or 2p. Let A, B, C = 2^64·v0/p, 2^64·v1/p, 2^64·v2/p, reals in (0, 2^64) with A + B +
C ∈ {2^64, 2^65}. F0 is the top limb of x and F1 the truncated product. In 64-bit words the row forms

```text
S = F0 + F1 mod 2^64        d = F0 - F1 mod 2^64 (SUBS: carry iff F0 >= F1)        m = min(F0, F1) (CSEL) e = m + S (ADDS: carry iff m + S >= 2^64) j
= 2 if e carries, else 1 if F0 >= F1, else 0
```

With F2* = 2^64 − S ∈ [1, 2^64] (2^64 when S = 0), e carries exactly when m ≥ F2*. The row falls back to the exact rule (`upstream-proof.md` section 1.1) when any of four windows holds:

```text
(W1) F1 >= 2^64 - 16     (W2) S >= 2^64 - 17     (W3) e mod 2^64 >= 2^64 - 33     (W4) d < 16
```

**One-sided errors.**
- F0 = A − a with 0 ≤ a < (r + C0)/2^192, where r = x mod 2^192. Indeed 2^256/p = 1 + C0/p, so A = (x/2^192)(1 + C0/p) = F0 + r/2^192 + x·C0/(2^192·p), and x < p. So a < 2, and a < 1 unless r > 2^192 − C0.
- F1 ≡ B − b (mod 2^64) with 0 ≤ b < 16. The true bits 320..383 of x·G are T ≡ F1 + c (mod 2^64) with 0 ≤ c ≤ 14 (they lie in [F1, F1 + 14], `upstream-proof.md` section 1.1), and T = ⌊B − ε⌋ mod 2^64 with 0 ≤ ε < 2^-64 (x·G/2^384 = xβ/p − ε/2^64). So b = (B − ⌊B − ε⌋) + c < 2 + 14. If B − b < 0, then F1 = B − b + 2^64 ≥ 2^64 − 15, which (W1) excludes. Outside (W1), F1 = B − b as real numbers.
- Then F0 + F1 = A + B − a − b. Since A + B + C ≡ 0 (mod 2^64), F2* ≡ C + a + b (mod 2^64) with 0 ≤ a + b < 18. If C + a + b > 2^64, then F2* = C + a + b − 2^64 ≤ 17, so S ≥ 2^64 − 17, which (W2) excludes. Outside (W1) and (W2), F2* = C + a + b.

So F0 ≤ A, F1 ≤ B and F2* ≥ C.

**Exactness.** Assume none of (W1) to (W4).
- j = 2 (m ≥ F2*): A = F0 + a ≥ F2* + a = C + 2a + b ≥ C and B = F1 + b ≥ F2* + b = C + a + 2b ≥ C. The v_i are distinct, so C is the smallest. No gap test is needed.
- j = 1 (F0 ≥ F1 and F1 < F2*): e = F1 + S < 2^64, so (W3) gives F2* − F1 = 2^64 − e ≥ 34, and (W4) gives F0 − F1 ≥ 16. Then A − B = (F0 − F1) + a − b > 16 − 16 = 0 and C − B = (F2* − F1) − a − 2b > 34 − 2 − 32 = 0.
- j = 0 (F0 < F1 and F0 < F2*): (W3) gives F2* − F0 ≥ 34, so C − A = (F2* − F0) − 2a − b > 34 − 4 − 16 > 0. And B − A = (F1 − F0) + b − a ≥ 1 − a, which is positive whenever a < 1.

So the rule returns argmin(v0, v1, v2) unless j = 0 and r > 2^192 − C0. Under the kernel rows' random-point premise, a phase-1 arrival has r > 2^192 − C0 with probability at most (C0 − 1)·2^64/(p − 1)
< 2^-159.99 (C0 − 1 low parts for each of the 2^64 top limbs). Section 5.3 charges the failure row (S_ops + 2^61)·C0·2^64/(p − 1) < 2^-33.28. Phase-2 and re-walk operations evaluate the exact rule
inside their 600-unit price. Part A checks `fast_j`, the bit-exact model of this rule, against the exact rule on 30,000 random, 10,000 near-tie (k = 150..249, including values near 0 and v0 + v1 near
p) and 1,000 low-r inputs, with 0 mismatches.

**Fallback rate.** (W1) covers 16 values of F1, (W2) 17 of S, (W3) 33 of e and (W4) 16 of d. If βx mod p is equidistributed (`upstream-proof.md` 1.1), each has probability at most its width over 2^64,
twice that for (W3) since e is 2F0 + F1 or F0 + 2F1: below 2^7/2^64 per arrival, so the price 0.5 and the 2^70 "other" row of section 5.3 cover fallbacks.

### 3.2 The sequence

In the arrival pass x16 holds the constant 33 (the chain's CCMN); constants are set once per batch.

- **Order.** F1 (13); ADD S, SUBS d, CSEL m, CSEL of b1 or b0 on the SUBS carry, ADDS e, CSEL of b2 or the previous choice on the ADDS carry; TBNZ on bit 0 of y0 splits the arrival into two paths that rejoin at the BFI.
- **Key.** K is the even one of y and p − y (13c3295f). Even y: K0 = y0, K1 = y1; CMP y1, #0 sets Z for the distinguished test and ORR forms the filter word y0 + 1. Odd y: SUBS from the register holding p0 + 1 gives K0 + 1 into a free register (y0 stays for the BFI), because y0 is odd and p0 is odd, so y0 ≠ p0 + 1 and the borrow is that of p0 − y0; this is the filter word. Then SBCS from 2^64 − 1 gives K1, with Z set exactly when K1 = 0. Both paths set flags only after the CSELs have used the decision flags.
- **Chain and branch.** CCMN F1, #16, #2, NE; CCMN S, #17, #2, LO; CCMN e, x16, #2, LO; CCMP d, #16, #0, LO; CCMP of the filter word with K_cp's limb-0 word, #4, HS; B.NE. A CCMN test passes with C clear, the CCMP d test with C set; a skipped test carries the failing C forward (#2, #0) and #4 sets Z. So B.NE jumps to the BFI exactly when there is no distinguished key, window or filter hit; otherwise the rare path after it runs.
- **Address.** One BFI, shared by both paths and both copies of loop A, writes bits 0..32 of y0 into bits 9..41 of b_j. Even y: y0 = K0, so bit 9 gets s = 0 and bits 10..41 get h, bits 1..32 of K0 (R = 2^32). Odd y: K = p − y gives y0 ≡ p0 − K0 (mod 2^33) whatever the limb-0 borrow, so y0 mod 2^33 = 1 + 2·((P′ − h) mod 2^32) with P′ = ((p0 mod 2^33) − 1)/2 = 0x7ffffe17: bit 9 gets s = 1 and bits 10..41 get ρ = (P′ − h) mod 2^32. The s = 1 half of row ρ therefore holds −Tab[(P′ − ρ) mod 2^32][j]; ρ ↦ P′ − ρ is an involution, so every walk reads the entry of its class's (h, s, j) and the class walk is that of section 5.1. b_j has zeros there, so the result is base + 1024ρ + 512s + 64j, rows of 1024 bytes; an entry's x + 2C0 sits at +192 (pass 2 reads [p, #192]).
- **Rare path.** Its own action at 0.025 units, charged first: LDR, SUB, STR of 25 milli-units on the saved meter (3); K1 recomputed from y0, y1 and p0, with CBZ to the report (at most 5); the four window tests again and B.LO to the fallback (5). For a filter hit it forms K2 and K3 (at most 5: LDP of y2, y3, then for odd y SUBS, SBCS, SBCS, SBC of p − y from y0, which is still held), loads K_cp's other two words (LDP) and runs CMP, CCMP, B.NE to the BFI (4): at most 22. In a pair it keeps x0 and x21, delivers its walk's address in that walk's register (x21 or x30) and finds its K_cp entry at a fixed offset from x3 (in the LDP immediate). Equal fingerprints are a true repeat and go on to the traversal; the report and the fallback are charged their prices (40.4, 0.5) first. The distinguished test precedes the filter, as in section 5.1.
- **Filter word.** K_cp's limb-0 word holds K0 + 1, odd; clearing stores 0, even, so a cleared checkpoint never hits, and for a different class the word is pseudo-random and odd: false hits at most 2^-63 per arrival.
- **Loop end.** x0 holds n: an odd n runs walk 0 alone first (own LDR of the K_cp word, STR of its address); AND x18, x0, #62 gives twice the left-over pair count mod 32, CBZ skips to the main loop when it is 0, else those pairs run one per iteration (SUBS x18, x18, #2, B.NE) and TBNZ x2, #26 ends it. The main loop is 32 straight-line pairs (LDP, two walks, STP) closed by one TBZ x2, #26: walk 0's STR advances x2 by 8 and each pair's STP by 16, so the pairs end at &tp[n], tp's first address with bit 26 set (section 2.5). With n = 2^23 (2^22 pairs, a multiple of 32) a round runs only TBZ, AND and CBZ before the main loop.
- **Test stub.** `kernels/arrive.s` (inline) and `kernels/secp_fuse.s` (after RET, section 2.5) log the rare case in place of this sequence; the plain path keeps its one branch.

Other rows (2d25756c's):
- **F1 by MADD.** G = ⌊2^384·β/p⌋ (`upstream-proof.md` 1.1); F1 = Σ low halves of x_i·G_(5−i) + high halves of x_i·G_(4−i) mod 2^64: 4 UMULH, 3 ADD, 4 MADD (xd = xa + xn·xm). G_1..G_5 are loaded once per batch.
- **Anchor rounds and the cap.** Every walk makes exactly one operation (a step, a doubling or a traversal step) per batch round. Rounds with t ≡ 0 (mod 256) run `pass1f`, the others `pass1p`. In rounds with t ≡ 0 (mod 8192) the processor runs the anchor variant of the arrival pass. For every walk it forms K2, K3 and the word K0 + 1 (7), writes K_cp ← fp(K) (STP + STR, the checkpoint-write action) and adds one to the trail's count of anchor periods. It abandons the trail when the count reaches 20·2^51 (LDR, ADD, STR, CMP, B.HS: 5). A trail's first anchor comes within 8,192 operations of its start, so it is abandoned before L_max = 20·2^64 operations. A walk starting a trail records φ = t mod 8192 (13 bits) in the trail record; its anchors fall where (length + φ) ≡ 0 (mod 8192) and its compares where (length + φ) ≡ 0 (mod 256). So a re-walk or phase-2 continuation anchors and compares where the live trail did. Anchor rounds run this arrival pass and the unfused pass 1, 31.015625 + 13.5078125 = 44.5234375 instead of `pass1f`'s 41.0078125 per walk (pass 1 keeps the 7-instruction fe_sub). Section 5.3 charges the variant's 12 instructions plus that 3.515625, 15.515625 = 993/64 per walk per anchor, as a hard row, and the checkpoint writes (2 per anchor) as another.

## 4. Prices and per-batch overhead

- **Ordinary step:** 37.5078125 + 169.50390625 + 337.50390625 + CORR + 24,552/2^23 = **549.0215812921524 + 24,552/2^23 = 549.024508**, with CORR = 37798701/2^23 the correction deposit (below), priced as that exact rational; compare rounds add pass1f − pass1p = 3.5 per walk, a hard row (section 5.3).
- **Traversal step:** step + 24 + 2 + 7 + 6 = **588.0215812921524 + 24,552/2^23**: 24 to track and resume the minimum-fp member, the charge on its own path (SUBS of 39 milli-units from x5 after the passes, B.CC to the refill, before the traversal's first step in the next round), 7 forming K2 and K3 for the fingerprint the minimum compares, and 6 for its own end test against K_cp in every round (LDR, LDP, CMP, CCMP, CCMP, B.EQ).
- **Doubling:** the batched doubling shares the product chain and inversion. Its forward part, after prod's n walks and the zero test, forms dx = 2y and dy = 3x^2 and both products by acc (acc·dx and e = acc·dy); its backward part, with pass 2, forms λ = R·e and R ← R·dx. Its arithmetic is 4 fe_mul with the loops' rare exit (4·84 = 336; the four products) + fe_sqr of x (63) + dx = 2y by one weak addition (9) + 3x^2 by two weak additions with a second carry fold (2·14) + the fused x3 = λ^2 − x − x (75) + x − x3 (7) + the fused y3 = λ·(x − x3) − y (91) = 609, with the loops' rare-exit kernels of section 2.5 (fe_mul 84, fe_sqr 63, x − x3 7; fkiene's 0a1c5772); each exit is its own 0.02-unit action, charged before its tail runs. A doubling takes at most 5 carry exits and 1 fe_sub exit, fewer than an addition's 12 carry and row exits and 3 fe_sub exits. Its slot in loop A, prod and pass 2 still runs that round's addition (below), so the rare-exit row of section 5.3 charges one more operation's exits for each expected doubling, and the doubling's stubs test their charge (SUBS of 20 from x5, the word its surcharge charges, B.CC to the refill; with a tail of at most 3 and B, at most 6 ≤ 20): they never draw on the round's reserve (meter, below).
  - **Traffic: 30.** Forward it loads x and y (4 LDP) and stores dx and e (4 STP); backward it loads e, dx, x and y (8 LDP) and stores x3 and y3 (4 STP): 20, plus at most 10 for pointer updates and loop control. R and acc stay in registers.
  - **List moves: 34.** Upstream's 16 loads and 16 stores of the walk's x and y to and from the doubling list, paired into 8 LDP and 8 STP, plus 18 other instructions.
  - **Escape logic: 3.** The clearing of the history (one STR of 0 into K_cp's limb-0 word) and the charge of the surcharge on its own path (SUBS of 677 milli-units from x5, B.CC); upstream's 30 also paid for the distinguished-bit scan and trail teardown, which checkpoint-only detection does not run.

  The arrival that ends a traversal runs inside loop A at the walk's slot in a straight-line group, so that slot still runs its addition (result discarded) while the walk doubles on the doubling list, and the doubled point is the walk's member for the next round, whose loop A arrives it. The doubling is the displaced step plus 609 + 30 + 34 + 3 + 24,552/2^23 = step + 676.0029, priced at **step + 677**, a surcharge of 677 milli-units (the run's `arithmetic_bound` and the certificate's 639 are 609 + 30). The weak additions, the 18 and the list moves are those of `upstream-proof.md` section 5.3.
- **Zero-denominator test, once per batch.** prod multiplies the n canonical dx_i into acc_(n−1); p is prime, so outside the rare event of section 2.2 acc_(n−1) ≡ 0 (mod p) exactly when some x(T_i) = x(W_i). Before the inversion acc_(n−1) is canonicalised and tested (ADDS, 3 ADCS, 4 CSEL, 3 ORR, CBZ: 12 per batch). Only when it fires, with probability < 2^-114 over the whole run, does the batch take the cold path: the unfused pass 1 is re-run on the addresses loop A (or an anchor round's arrive.s) stored, with the per-walk test (each walk with dx_i = 0 is abandoned and its dx_i replaced by 1), and then prod, which forms every dy_i and e_i again from the points. The re-run, at most 2^23·(16 + 4 + 172 + 3) = 2^23·195 < 2^31 instructions (its prod's dy borrows included, so no count goes to A) and at most 7 stubs per walk (140 milli-units), is its own action priced 2^22 units, charged from M before it runs, in the "other" random row; its stubs charge x19, which the processor saves, sets to 140·n from that price and restores after the re-run, never the round's reserve.
- **Per-batch fixed work:** inversion I ≤ 24,322 (`secp_fe.s`, `upstream-proof.md` section 5.1–5.2); loop A's prologue, tail-count and group tests, count and epilogue (33 = 22 + 3 + 8, executed, in both copies); prod's 28 and pass 2's 28 (the function minus its loop bodies, by `a64.py`, each with MOV x18, xzr and MOV x0, x18); 64 for the passes' call and return; at most 8 for prod's call (MOV of x0–x3, BL); 12 for the zero test, 3 for the context words p0 + 1 and 2^64 − 1, 7 for the meter, 3 for the round counter and the anchor-round test (ADD, TST, B.EQ), 2 for the compare-round test (TST, B.NE) and 4 for the arguments of the fused call; 38 for the account A (below): 24,552 < 2^15. All 2^23 walks form one batch: every walk makes one operation per round. With n = 2^23 exactly divisible by both group sizes, no loop runs a tail walk, so a round's fixed work is 24,552 instructions (the tail rows are not charged). The batch's scratch is 72 bytes per walk (slots and pointers), under 2^30 bytes in all, inside the memory bound of `upstream-proof.md` section 5.3.
- **Meter per round** (a9297fb0's reserve, kept in its own word). The round charge is c = 885.0215812921524·2^23 + 24,552 milli-units: 549.0215812921524 (CORR included) and a reserve of 336 per walk, plus the fixed work, plus 3.5·2^23 in compare rounds; the two values stay in the callee-saved x21 and x23, and each branch of the compare-round test charges its own. Before a round the processor charges c from its block remainder x5 (SUBS, B.CC to the refill), adds 2^23 to its operation count x6 (ADD), stores the remainder in the meter word M (STR x5, [x22]; the callee-saved x22 holds &M all phase) and sets the reserve x19 = 336·n (MOV); after the passes it reloads M and adds back the unspent reserve (LDR, ADD): 7 instructions, inside the 24,552. A plain round is charged the step price per walk plus what the reserve paid. The reserve x19 (saved at [sp] by each pass's prologue, restored by its epilogue) pays only the unchecked charges: the stubs' 20 milli-units each and a rare-path entry's 25, at most 15 exits and one entry per walk, so at most 325·n < 336·n per round, and x19 stays non-negative without a test. Between the STR and the reload every other action (report, start, fallback, the zero test's re-run) charges M with a test (inside a pass via &M at [sp, #24], the caller's x22: LDR, LDR, SUBS, B.CC to the refill, STR), and after the reload the surcharges (677 or 39) charge x5; none charges the reserve, so no checked action can spend what a later unchecked charge needs. An anchor round (section 3.2) charges c + 18·2^23 (MOVZ #0x900, LSL #16, ADD: the row is charged the exact 15.515625), paying its variant's 16 and checkpoint write's 2 per walk before they run. A round's charge is below 2^33 milli-units (an anchor round's conservative upper bound (552.530188 + 336 + 18)·2^23 + 24,552 < 2^32.83; round 0's adds 9n, below), far below a 2^40-unit block, so one refill covers it and every walk is covered at the step price before it runs. Rare-path entries and rare exits are paid on their own paths and the final arrival pass is its own action (section 2.5), so every action is charged before it runs. Blocks are reserved, refilled and counted as spent exactly as in `upstream-proof.md` section 3, and phase 1 still stops at a refill once the operation count reaches S_ops. Even 2^20 processors' last blocks would cover fewer than 2^20·2^40/0.5490 < 2^60.87 further operations, so phase 1 performs fewer than S_ops + 2^61 operations (sections 3 and 5.3).
- **Correction account A** (2.4, 5.4). CORR = 37798701/2^23 = PREM + 813/2^23 per walk (σ ∈ {dx, dy, u}; PREM = 3·(3/2 + 2^-10 + 2·2^-11), and the 813/2^23 margin is checked against the Bennett/Ville variance bound) is in the step price. A is a 192-bit triple at &M + 64, +72 and +80, a line no other word uses, read and written only by the processor between passes, never by a stub; it never exceeds the work charged (< 2^136 milli-units), so it cannot wrap. Before the passes the processor adds CORR·n = 37798701 to A and stops the attack with failure unless A ≥ 9n (LDP, LDR, MOVZ, MOVK, ADDS, ADCS, ADC, MOVZ #0x480, LSL #16, CMP, CCMP x1, #0, #0, LO, CCMP x2, #0, #0, EQ, B.EQ, STP, STR: 14); after loop A, prod and pass 2 it subtracts 3k, k the count in x0 (ADD x0, x0, x0, LSL #1, LDP, LDR, SUBS, SBCS, SBC, STP, STR: 8 each; anchor rounds run pass 1 and skip loop A's): 38. A slot runs one dx, dy and u per round (a doubling round's discarded addition included), so a round has at most 3n borrows, 9n instructions, which the check prepays: A never goes negative, and every correction is paid from deposits charged before it, on every run. A starts at 9n, in round 0's charge; after round 0 the processor charges Λ_p = 2^94 milli-units on the 192-bit meter as one action, a hard row (5.3), and adds it to A. Unused A counts as spent; the stop is a failure row (5.4), never extra work.

Every other price (reports 40.4, probes 0.05, starts 2^13, fallbacks 0.5, phase-2 and re-walk operations 600, refills 0.03, words, attempts, store initialisation) is that of `upstream-proof.md`; the checkpoint write
(0.002), the rare-path entry (0.025) and the rare exit (0.02) are their own actions; the anchor-round work is a hard row (section 3), charged with its round.

## 5. The walk at R = 2^32, H1 and the exact budget

### 5.1 The walk change

The walk is that of `upstream-proof.md` section 2 with six changes: the table has R = 2^32 rows instead of 2048, M = 2^23 walks instead of 2^20 (one batch, n = M, on one processor; rows sized for up to 2^20 processors are upper bounds), cycles are detected by the checkpoint alone, the
class key is the even one of y and p − y, anchors fall where (length + φ) ≡ 0 (mod 8192), and the checkpoint is compared where (length + φ) ≡ 0 (mod 256). Everything else is unchanged: the canonical
index, the twisted table Tab[h][t] = [u_t^{-1}]M_h, the distinguished classes (limb 1 of K zero, θ = 2^-64) and L_max = 20·2^64. At each class arrival, in order:

```
if K distinguished: report (K, a, b); end trail if length >= L_max: abandon trail if (length + phi) mod 256 == 0 and fp(K) == K_cp: traverse the cycle
once more (<= 32 steps, else abandon), resume at the member of its minimum-fp class remembered during the traversal, ESC; clear history; continue if
(length + phi) mod 8192 == 0: K_cp <- fp(K)   (phi: the trail's start round mod 8192, section 3) P <- P + Tab[h][t]
```

"Clear history" stores 0 in K_cp's limb-0 word (section 3.2). ESC is the doubling escape P ← [2]P, a class map.

- **Merging lemma.** Let trails A and B arrive at the same class C, with any histories. Until one escapes, both apply F and visit the same classes. A trail escapes only after K = K_cp: K_cp is set at an arrival and the arrivals up to the next refresh are F-steps, so a repeat of K_cp is a cycle Z of F with no distinguished class, and the traversal escapes only if |Z| ≤ 32. Every escape from Z leaves from its minimum-fp class and goes to E(Z) = class(2c(min Z)), which depends only on Z. If A escapes from Z, B follows the same F path and enters Z, sets K_cp inside Z within 8,192 arrivals, compares there again lcm(|Z|, 256) ≤ 7,936 arrivals later, before its next anchor, and escapes to the same E(Z) after a traversal of |Z| and one doubling. So A and B meet again after at most 8,191 + 7,936 + 32 + 1 = 16,160 extra operations of one of them, and by induction reach the same distinguished class unless one first hits L_max or an exceptional addition (failure events). History clears at every doubling, so K = K_cp witnesses a real cycle of F.
- **Cycle rates.**
  - A 2-cycle needs the next class to use the same partition index with the inverse twist: rate 1/(6R) = 2^-34.585 per step.
  - Other cycles of length 3 to 32 need every partition index used at least twice: rate at most Σ_L Σ_b S2(L, b)·R^-(L−b) < 2^-61.99999 at R = 2^32, over L = 3..32 and b = 1..⌊L/2⌋, where S2(L, b) counts the partitions of L steps into b blocks of size ≥ 2 (1, 1, 3, 15 for (L, b) = (3, 1), (4, 1), (4, 2), (6, 3)). `certificate.py` recomputes it exactly.
  - Cycles longer than 32 end their trail (the traversal cap if compared, else the anchor-count cap) and are rarer still: below 2^-376 per step for L = 33 to 16,384 at R = 2^32 (exact to L = 256; above it S2(L, b) ≤ C(L, b)·b^(L−b) and b ≤ 8,192, so each L adds at most 2^L·(8,192/R)^⌈L/2⌉). Such trails number below 2^-232 per run and cost below 2^-100 units in expectation; `certificate.py` asserts both.
  - A cycle with some partition index used once closes by chance: a fresh step lands on one of the trail's t earlier classes with probability at most t/N under H1's pseudo-random key fields, at most L_max^2/(2N) = 2^-117.8 per trail (H1's non-fruitless-return clause), under 2^-54 per run; their doublings and traversals cost under 2^-47 units, in the 2^70 "other" row (5.3).
  - So the expected number of doublings is below (1/(6R) + 2^-51.99999)·(S_ops + 2^60) + 2^-54 < 2^-34.5849·(S_ops + 2^60). A 2-cycle is traversed in 2 steps, other short cycles in at most 32, so the expected number of traversal steps is below (2/(6R) + 32·2^-51.99999)·(S_ops + 2^60) + 32·2^-54 < 2^-33.5849·(S_ops + 2^60). These use only linearity of expectation over the pseudo-random partition and twist of each new class; 5.5's budget also
  uses premise P1.
- **Detection delay and pathwise waste.** A cycle Z is caught within 8,191 + lcm(|Z|, 256) ordinary steps of entering it (8,449 with the 2-step traversal for |Z| = 2; 16,158 for the largest |Z| ≤ 32), ordinary work counted in S_ops, with its arrivals in T_c. The expected waste rate is 8,449/(6R) + 16,158·2^-61.99999 < 2^-21.54, but 5.2 needs a pathwise bound. 5.5's Ville lemma is one, and P1 now covers every catch of length ≤ 32 (5.5's census measures doublings of any cycle length): outside a 2^-140 event, phase 1 catches fewer than (1 + 2^-11)·q′·(S + 2^61) cycles of length ≤ 32, q′ = 1/(4R) (premise P1, 1.5 times the uniform 1/(6R)). Charging every such catch the worst-case 16,158 wasted arrivals (no expected count is used) gives the realised waste w ≤ (1 + 2^-11)·16,158/(4R) < 2^-20.0193 < 2^-20.00; the 2^-140 event is in 5.3's 1.1e-7 budget and `certificate.py` checks the inequality.
- **Table.** 2^32 rows of 1024 bytes (section 3; the odd half of each row permuted as in 3.2), 2^42 bytes, aligned to 2^42. Each entry also holds x̃ = x + 2C0, so a row is redrawn (metered, as for M_h = O) when one of its three x-values is ≥ p − C0 (below 3C0/p < 2^-222 per row); the redraw depends on M_h only, so section 5.2's premises hold. Even with a copy on each of the 2^16 nodes, each node stays below 2^56 bytes and peak memory below 2^72 bytes (`upstream-proof.md` 5.3, "Memory"). Building it takes 2^32 draws of [c]P + [d]Q at ≤ 2^13 units each plus twists and x̃ fields: below 2^45.1 units, in the hard extras (< 2^46) and in preprocessing, below its 2^59 cap.
- **Key and partition.** The class key is the even one of y and p − y; h is bits 1..28 of K (section 3), disjoint from the distinguished bits 64..127, routing bits 128..143 and bucket bits 144..193. The anchor phase φ only shifts where a walk's anchors and compares fall; each walk is anchored once per 8,192 of its own operations and compared once per 256, as the lemmas use.

### 5.2 H1 at R = 2^32

H1 keeps its form: P(T_c ≤ S_ops) ≥ 1 − exp(−(S_ops/κ)^2/(2N')), with N' = N·R/(R − 1) and now κ = 1.000041, on the hazard census of 12ccd81d, which this package ships (`census_hz.py`,
`census-hz-q29.json`). The rare exits change no result; the stored member (2.5) does not either. T_c counts every class arrival, including repeated arrivals inside
fruitless cycles, during the detection delay and in traversals; its other parts are as in `upstream-proof.md` section 4 (pseudo-random key fields, the non-fruitless return bound, a uniform twist index
independent of the partition index).

**Why T_c runs slow on a toy.** T_c ends at the first arrival at a class visited by *another* trail: with D classes visited, `own` by the current trail, the collision chance is (D − own)/N, so E[T_c] ≈ √(πN/2) + ā, ā = (1 − θ)/θ: +0.45 % on q29 dp_bits 5, but ā/√N = 2^-62.7 in production (at most S_ops·L_max/N = 2^-58.4
of the hazard).

**Hazard census.** `census_hz.py` runs the walk of section 5.1 on the toy groups (2d25756c's walk code) under the rule (1, 256). F and E(Z) are the same under the production rule, which changes only
when a cycle is caught (`evidence.md` item 7), hence only the repeated arrivals, the waste (1 + w) covers.
- At every fresh-sourced arrival (a trail start, an escape, or a step from a class first visited at the previous arrival) it adds the random-mapping hazard (D − own)/N; a merge is an arrival at another trail's class.
- For a uniform random mapping the hazards are the merges' compensator, so E[merges] = E[Σ hazards], and ρ = merges/Σ hazards measures the walk's collision intensity against the law with the own-trail exclusion removed exactly. Repeated arrivals add no hazard: they are the waste w.
- In multi mode a run of 64·√N arrivals continues after each merge, which ends only the merging trail; in tc mode each trial stops at its first merge (T_c). `--pure` runs uniform random functions through the same trail logic: the control.

`census-hz-q29.json` holds each multi-mode run's ρ and κ_eff = ρ^-1/2; `census.md` gives the pools, the control at T_c and the first-collision distribution, which
follows the random-mapping law at the quantile H1 uses as well as in its mean. The rule below was fixed before each run.

**κ = 1.000041.**
- In production a run of S arrivals has at least S/(1 + w) fresh ones, with w = 2^-20.00: 5.1's pathwise bound is w ≤ 2^-20.0193 outside the 2^-140 event, so this is the realised waste fraction, not an expectation; `certificate.py` asserts the bound from `census-hz-q29.json`.
- Its own-trail term is below 2^-58.4 of the hazard, and its collision intensity is ρ.
- So P(T_c ≤ S) ≈ 1 − exp(−ρ·(S/(1 + w))^2/(2N)), which is H1 with κ_eff = (1 + w)/√ρ.
- κ rests on the multi-mode runs, against the law's ρ = 1, which the compensator makes exact for a uniform mapping. At R = 2^20 and dp_bits 5 the 1,024-table run gives κ_eff = 0.99978 ± 0.00010, so (1 + w)(κ_eff + 2 s.e.) = 0.99999. Pooled over the eleven walk runs, κ_eff = 0.999887 ± 0.000066 and the same bound is 1.000019. κ = 1.000041 is the smallest six-decimal value at or above the pooled bound 1.000019 and within one standard error of every walk run's (1 + w)·κ_eff; the closest is q29 at R = 2^25, dp_bits 3, whose excess is 0.0002800, 1.000 standard errors. The three-run pool in `census.md` is scatter only (χ²/df = 2.5).
- The 512-table ideal control: ρ = 0.99976 ± 0.00028 (0.86 s.e.).
- No walk run's (1 + w)·κ_eff exceeds κ by one standard error. The closest, q29 at R = 2^25, is 0.0002800 above κ, 1.000 standard errors; its (1 + w)(κ_eff + 1 s.e.) = 1.000601 is above κ.
- N' = N·R/(R − 1) is kept. The twisted walk's 1/(6R) anti-collision is inside ρ, so N' over-covers it.
- This is the standard random-mapping heuristic; no theorem bounds the constant.

**Other checks**:
- **The hosted toy pipeline** draws a fresh table per trial at R = 2048 and adds the toy's own-trail term ā to its scaled budget: evidence of the whole procedure, not of κ. Its 2-cycle rate 2^-13.6 is 2^17 times production's, and a caught 2-cycle holds a walk past the toy cap 640, so the trail is abandoned: **640·2^-13.6 = 0.052** is the abandonment the toy rate predicts (`abandoned_visit_fraction`: 0.0502 and 0.0481 in two local runs of the declared input). In production the same product is **(1 + 2^-11)·16,158/(4R) < 2^-20.0193**, the pathwise bound on every catch of length ≤ 32 (5.1, 5.5): L_max = 20·2^64 is far above 8,447, so production buys a delay inside the meter and abandons nothing. The rule check (`evidence.md` item 7) is toy-scale only.

### 5.3 The exact budget

H1 and κ = 1.000041 remain the disclosed heuristic (5.2). Every catch of length ≤ 32 is paid from B_c (5.5), the Markov headroom is h = 0.0000116 and the budget for
every other failure event is 1.1·10^-7, so p0 = 0.39001172 keeps the success floor at 0.39000001 > 0.39 with W below 2^125.834346. Values below are rounded up.

```text
N' = ((r-1)/6) * 2^32/(2^32 - 1)
S0 = sqrt(2*N'*ln(1/(1-0.39001172))) < 2^126.6992709
S_ops = 1.000041*S0 < 2^126.6993300
phase1 = (549.0215812921524 + 24552/2^23)/1000*(S_int + 2^60) < 2^125.8342725
hard extras < 2^110.878832
compare rounds 0.0035*ceil((S_int + 2^61)/256) 2^110.540901
anchor rounds 0.015515625*ceil((S_int + 2^61)/8192) 2^107.689196
checkpoint writes 0.002*ceil((S_int + 2^61)/8192) 2^104.733546
correction seeds 2^20*(2^94 + 9*2^23)/1000 2^104.034216
overflow insert/lookup and remaining hard rows (table < 2^46) 2^107.20
catch budget B_c = 0.805*(1 + 2^-11)*(S_int + 2^61)/(4R) < 2^92.387095
residual random extras E_r < 2^93.769828
doublings 0.677*2^-61.99999*(S_ops + 2^60) 2^64.137
traversals 0.039*32*2^-61.99999*(S_ops + 2^60) 2^65.019
rare-path entries 0.025*(2*2^-61.99999 + 2^-56 + C0/2^64)*(S_ops + 2^60) 2^89.377
rare exits 0.02*(25*C0/2^64 + 2*2^-63)*(S_ops + 2^61 + (1/(6R) + 2^-61.99999)*(S_ops + 2^60)) 2^93.699
starts 2^13*(S_ops/2^64 + 2^23) 2^75.699
followers (5.5) 2^80
words, attempts, probes, small batches, fallbacks, zero-test re-runs, final arrival 2^70
allowance A_r = ceil(E_r/0.0000116) < 2^110.165344
W_ledger = phase1 + hard + B_c + A_r < 2^125.834346
```

The rare-path row counts filter hits (true ones of cycles of length 3 to 32, at most twice their rate, an over-count; false ones below 2^-63), window fallbacks (below 2^-57),
distinguished keys (2^-64) and, in loop A, a C0/2^64 term per arrival kept as an upper bound on its distinguished candidates (at most 3·2^-64 in pass1p; odd borrows do not enter). The rare-exit
row counts, per operation:
- the three fe_sub exits (dx in loop A, dy in prod, u), each taken only when the low limb of a − b wraps below C0: at most C0/2^64 each;
- the two carry exits of x3's and y3''s folds, each taken only when s1 ≥ 2^64 − 2: at most 2^-63 each;
- six row exits (three in prod, three in pass 2; C0/2^64 each) and four fold exits (e and acc in prod, λ and R in pass 2; 4C0/2^64 each), section 2.5.

So the expected number of exits is at most 25·C0/2^64 + 2·2^-63 < 2^-27.356 per operation, at 0.02 units each, over fewer than S_ops + 2^61 operations, plus as many again per expected doubling,
whose slot also runs its addition (section 4). The rates rest on the random-point premise of the kernel rows. Every exit is charged before its tail runs (section 4) whatever the inputs, so the
premise affects only this expectation. E_r's rows keep 5.1's expected rates; B_c and the realised waste use P1. `certificate.py` uses the exact rates and asserts E_r < 2^93.77; `experiment.py`'s float ledger recomputes them (log2 W_ledger ≈ 125.8343455) but is not a certificate.

In group-operation units W_ledger < 2^125.834346, the cap, which `certificate.py` bounds below; the next
lower six-decimal value, 125.834345, is below W_ledger and is not claimed. Against the headroom left after B_c (cap − phase 1 − hard − B_c > 2^110.189559), the Markov ratio is
0.0000114070 (exact, rounded up) < 0.0000116, so Markov's inequality on the nonnegative U gives P(U > A_r) ≤ E_r/A_r < h and the meter halts early with probability below h + 2^-140 (5.5). The 1.1·10^-7 covers every other failure event: the rows of `upstream-proof.md` section 4 (trail cap < 10^-7, preprocessing cap <
2^-30 and smaller), the three kernel rows of section 2 and the canonical-index row of section 3.1: the dropped carry 7·(S_ops + 2^61)·(C0^2 + 3C0)/p < 2^-62.4, at most 7 reductions per batched
operation; a doubling's 7 come on top at rate 2^-34.58); a doubling's negative fused fold (S_ops + 2^61)·C0/2^256 < 2^-97.2; a fused result in [p, 2^256) for x3 or y3, 2·(S_ops + 2^61)·C0/p < 2^-96.2; and a wrong j = 0, (S_ops + 2^61)·C0·2^64/(p
− 1) < 2^-33.28; and the account's rows of section 5.4, a table that is not good 6·e^-2048 < 2^-2951.9 and the early stop 2^20·2^-140 = 2^-120; and 5.5's 2^-140. Rare exits return inline results and add no failure row. Together these sum to 1.0103·10^-7. So

```text
success >= 0.39001172 - 0.0000116 - 1.1e-7 = 0.39000001 > 0.39.
```

`certificate.py` re-derives in exact arithmetic the kernel identities, loop counts (main paths also by `a64.py`), fixed work, tails, prices, CORR and the
account's seed, DKW and Ville rows, rates, B_c and its lemma's constants, w at q′, ledger, failure rows, κ bounds, stubs, reserve, round charges, re-run price, preprocessing and memory, and the
claim 125.834346 against 125.834345. The candidate certificate completed in the pinned isolated image and printed `CERTIFICATE: PASS`; the declared initial experiment ran separately.

### 5.4 The correction account: premise C and its census

**Premise C** (disclosed; trail averages only). Fix the table T, good as below, and an anchor phase φ. For a trail from fresh start coins let L be its
slot-steps (every operation of its slot, traversal steps and a doubling round's discarded addition included) and B_σ the slot-steps at which site σ
borrows (a < b as integers). Then E[B_σ] ≤ (½ + π₀_σ)·E[L], π₀_dx = 2^-10, π₀_dy = π₀_u = 2^-11. Nothing is claimed per step, no independence between
steps and no exact ½; the pads are π = 2π₀, a margin ρ′ = 813/2^23 per slot-step. Support:
- **dy.** ±W share the class and (h, j), and the s = 1 entry is −Tab[h][j] (3.2); with Y = y(M_h) the borrow is [Y < K] for even y_W, [Y > K] for odd: one sign borrows. So the rate is ½ whenever the held member's sign, that of −(W + T), is a fair coin given the class (the random-point premise); 
- **u.** x3 = x(W + T) is a fresh pseudo-random coordinate: P(x3 < x_W) = E[x_W/p] = ½.
- **dx.** x_T = β^−j·x(M_h) is a draw from column j of T (h pseudo-random given x_W: H1's key fields), so the rate Σ_j P(j)·E[F̂_j(x_W/p) | j] (F̂_j that column's empirical CDF) lies within max_j D_j(T) = max_j sup_v|F̂_j(v) − v| of E[x_W/p] = ½, whatever x_W's law given j (for j = 0, x_W is its orbit's minimum). T is good when max_j D_j ≤ 2^-11; by the DKW inequality on three columns of R = 2^32 independent rows, P(not good) ≤ 6·exp(−2R·2^-22) = 6·e^-2048 < 2^-2951.9, a failure row (5.3). The other 2^-11 of π₀_dx is model slack.
- **Cycles.** Steps trapped in a cycle before detection are under w < 2^-20.0193 of all (5.1); they move a site's average by under 2^-21.02.

**Concentration.** Order one processor's trails by start time; let X_k = 3·B_k − CORR·L_k, B_k the three sites' borrows over trail k. X_k depends only on
T, φ_k and trail k's own coins (merges do not couple counts; φ_k, the start round mod 8192, is fixed by the earlier trails G_(k−1)), so premise C gives
E[X_k | G_(k−1)] ≤ −ρ′·E[L_k | G_(k−1)]. A slot-step adds 3b − CORR ∈ [−a, 9 − a] ⊂ [−a, a], a = CORR, so X_k − E[X_k | G_(k−1)] ≤ c = 9·L_max and
Var(X_k | G_(k−1)) ≤ a²·L_max·E[L_k | G_(k−1)], as L_k < L_max = 20·2^64 (3.2). Bennett's supermartingale exp(λM − ψ_c(λ)V) (M the centred sum, V the summed
variances, ψ_c(λ) = (e^(λc) − 1 − λc)/c²) with Ville's inequality at λ = 2^-85: as e^t − 1 − t ≤ (t²/2)/(1 − t/3), ψ_c(λ)/λ·a²·L_max ≤ λa²L_max/(2(1 − λc/3))
= approximately 0.0000968181 < ρ′ = 813/2^23 ≈ 0.00009691715, the variance term is absorbed, and every partial sum Σ_(i≤k) X_i is at most 140·ln 2/λ < 2^91.60 except with probability 2^-140,
whatever the number of trails. After round 0, A before a round is at least 9n + Λ_p + CORR·n minus the processor's realised X; at most n trails are in progress and
a realised part exceeds its trail's whole X by at most CORR·L_max. As 2^91.60 + n·CORR·L_max < 2^91.60 + 2^93.49 < 2^94 = Λ_p, A ≥ 9n outside that event:
the early stop has probability at most 2^-140 per processor, 2^-120 over 2^20, plus the DKW row. `certificate.py` checks these inequalities exactly. The
work bound does not use premise C.

**Census** (`census_corr.py` on production secp256k1 with this walk rule, `census-corr.json`, described in `census.md`; runs and rule declared before
any ran). Rule: for R = 2^26 and for R = 2^22 separately, each site has q̂ + 4·SE ≤ ½ + π₀ and each table q̂_dx,T ≤ ½ + π₀_dx + 4·SE_T, SE the walk-cluster
standard error; a failure meant no package at these pads; R = 2^16 is reported, never judged. Production R = 2^32 is larger, so the dx pad only improves (DKW 6·e^-2048). Rates over 8 tables × 64 walks of 200,000 slot-steps (R = 2^16:
40,000; the support above does not depend on a step's trail age, apart from trapped steps), recomputed exactly by the declared run (`evidence.md` item 8):

| R | dx | dy | u |
| ---: | --- | --- | --- |
| 2^26 (decisive) | 0.500095 ± 0.000049 | 0.499974 ± 0.000050 | 0.499972 ± 0.000028 |
| 2^22 (decisive) | 0.499926 ± 0.000048 | 0.500069 ± 0.000052 | 0.499985 ± 0.000027 |
| 2^16 (diagnostic) | 0.501772 ± 0.000895 | 0.500722 ± 0.000679 | 0.500127 ± 0.000063 |

Every decisive q̂ + 4·SE is below 0.500292, below ½ + π₀ = 0.5009765625 (dx) and 0.50048828125 (dy, u), and every table passes: the rule holds.

### 5.5 The catch budget: premise P1

A caught 2-cycle costs at most 805 milli-units beyond its ordinary steps and rare exits: a doubling (677), two traversal steps (2·39) and
two filter hits (2·25; an over-count, as only the detecting compare hits). A longer caught cycle (length ≤ 32) costs at most 677 + 32·39 + 32·25 = 2,725; Section 5.3 pays every catch 805 from the
fixed row B_c and charges the excess of the longer cycles, with the other extras, to E_r and the Markov allowance (5.3).

**Premise P1** (disclosed; trail averages, as premise C). Fix the table T and an anchor phase φ. For a trail from fresh
start coins let L be its slot-steps (as in 5.4) and Z the cycles of length ≤ 32 it catches. Then E[Z] ≤ q′·E[L], q′ = 1/(4R), 3/2
times 5.1's uniform 1/(6R). Nothing is claimed per step. Z is the measured quantity: the P16 census records `doublings`, the escapes from cycles of any length, per slot-step (5.5, Census). At q′, (1 + 2^-11)·16,158·q′ = 2^-20.0193 stays within 5.2's w = 2^-20.00; since the catch bound is the lemma's Ville bound, this covers the realised waste of every catch, not just the mean of the 2-cycles, with the worst-case delay charged to each.

**Lemma.** Order phase 1's trails by start time (ties by processor, slot); G_(k−1) is the table and the coins of trails
1 to k − 1. Z_k and L_k depend only on T, φ_k and trail k's own coins, fresh and so independent of G_(k−1) (a trail's path ignores other trails;
φ_k is fixed by G_(k−1)), so P1 gives E[Z_k | G_(k−1)] ≤ q′·E[L_k | G_(k−1)] and X_k = Z_k − (1 + η)q′L_k, η = 2^-12, has drift at most
−ηq′·E[L_k | G_(k−1)]. Catches are at least 257 slot-steps apart (history clears at each doubling, and a compare follows
an anchor by 256), so Z_k ≤ L_k/257 + 1; as L_k < L_max, X_k − E[X_k | G_(k−1)] ≤ L_max/257 + 1 + (1 + η)q′L_max ≤ c = L_max/256 + 1 and
Var(X_k | G_(k−1)) ≤ v·E[L_k | G_(k−1)], v = (c + (1 + η)²q′L_max)·q′, as Z_k² ≤ c·Z_k. With 5.4's Bennett supermartingale and Ville's
inequality at λ = 2^-72, ψ_c(λ)/λ·v ≤ λv/(2(1 − λc/3)) < ηq′ absorbs the variance, so every partial sum Σ_(i≤k) X_i is
below 140·ln 2/λ < 2^79 except with probability 2^-140. The trails started in phase 1 have fewer than K + M·L_max
slot-steps, K = S_int + 2^61 (section 4; S_int is an integer ≥ S_ops), as each of the M = 2^23 slots leaves at most one
unfinished. So phase 1 catches fewer than (1 + η)q′(K + M·L_max) + 2^79 ≤ (1 + 2^-11)·Kq′ cycles of length ≤ 32 outside that event.

**Conclusion.** Phase 1's extras are at most 0.805·N_c + U units, N_c its catches of length ≤ 32 and E[U] ≤ E_r (its 2^80 follower
row over-counts; the long cycles' traversal and filter rows are in E_r). An early halt needs them above cap − phase 1 − hard ≥ B_c + A_r with B_c ≥
0.805·(1 + 2^-11)·Kq′, so N_c exceeds (1 + 2^-11)·Kq′ or U > A_r: probability below 2^-140 + h (5.3). P1 enters only
this row; the work bound holds on every run.

**Census** (`census.md`, `census-p1.json`; declared before any run). P16 measures P1 on production secp256k1:
`census_corr.py`'s doublings (of any cycle length, the quantity Z counts) per slot-step of 200,000-slot-step trail prefixes from fresh starts on
64 fixed tables at R = 2^16, against the rule's undercount F. The toys measure the per-arrival mechanism: Y = 1 when a step from C1 first reaches a
class C2 with 5.1's 2-cycle condition (same partition index, inverse twist), on census_hz.py's walk on q29. For D16 (R = 2^16, 2,048 tables), D18 (R = 2^18, 4,096) and P16
separately, with q = 1/(6R) and 16 groups: Y/(qF) + 4·SE ≤ 3/2 (SE the larger of the Poisson and group-spread errors),
χ² ≤ 16 + 5√32 and, on the toys, no Y = 1 followed by a class other than C1; a failure meant no package. Y/(qF) ± SE
(χ²): D16 1.0086 ± 0.0241 (14.39), D18 1.0752 ± 0.0451 (31.54), P16 1.0537 ± 0.0465 (14.85); the rule holds (`evidence.md`
item 9). Scope: trail prefixes at R = 2^16 against whole trails at R = 2^32 here, and toy groups (N ≤ 2^24.83, R ≤ 2^18).

## 6. Provenance, files and limitations

Design provenance: the construction, report store, ledger structure and F1 truncated product are the public line of `upstream-proof.md` (carried by 02feb284). bradcush's line gives the counted kernels, loops, rows, checkpoint-only walk, census tool and κ rule, the reserve word (445f9a1c), pass 2's shared rows (a00774d8), the permuted odd table half (886e6f54), λ = R·e and walk 0's slot (7dd2c482), dy in prod and R = 2^26 rows (02cc35d7), 2^23 walks per batch and compares every 256th (0cd0ad73), and the branch-only fe_sub corrections, account A, premise C and census (34a30cf5). kaankolcu's: rare exits, reserve and paired traffic (807e4bf5, a9297fb0), −(W + T) forms (4595bad7), fused arrival (cbea87bf), strided compares (4ac17b66), R = 2^26 rows (850e11fe). aburan28's: headroom retunes (b78fb4f6, 4b4defe7), arrival pairs (2bd168b5, 01488655), R = 2^22 (610c8e49). fkiene's: doubling recount (0a1c5772), R = 2^21 (6a6cb91c). Ulysseus's census reruns (daf2c741); b5de673a's fused x3 subtractions; P-256 provenance only. The direct predecessor 7680a8ca provides B_c, P1 and its census; 671caf26 the 256-walk grouping. This proposal changes the allowance, p0, R, κ and the waste argument. No prior acceptance, review outcome or score is a premise.

Files: `kernels.md` (2.1–2.3), `evidence.md`, `upstream-proof.md`, `experiment.py`, `certificate.py`, the `census_*` scripts and `census-*` JSON, `census.md`, and `kernels/`.

The limitations are upstream. H1 is a random-mapping heuristic supported at toy scale; κ = 1.000041 rests on the toy-group hazard census. The step price is the executed
main-path count of pass1p, prod and pass 2, CORR and the round's fixed work over 2^23 walks; compare rounds' 3.5 per walk is a hard row. Rare exits and the kernel
and canonical-index failure rows rest on the random-point premise; the corrections' rates on premise C (symmetry and census, no theorem), only through the early-stop row; the
catches of length ≤ 32 on P1 (5.5), only through the early-halt row. The attack is generic and impractical, not a claim about deployed secp256k1.
