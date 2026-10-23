# bls12-381-g1-ecdlp: Cheon's reduction with an exact baby-step giant-step

Status: analytic worst-case bound for exploratory review. Success is an exact count, not a probability model:
there is no random walk, no distinguished point, no hashing and no collision heuristic. Every number below is
recomputed by `experiment.py` (part 1), which prints it.

## 1. Inputs used and the range of n

The challenge gives T = A_1 = [s]P0 and A_i = [s^i]P0 for i <= n, with s uniform in 1..r-1. The attack reads
exactly two points: B = A_d with d = 2^20, and T. Everything else is public profile constants. So W(n) is the same
constant W for every n with 2^20 <= n < r, and the claim is made at the scoring value n = 2^20. For n < 2^20 the
bound is not claimed.

Facts checked by `experiment.py` with exact integer arithmetic:

- r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 * 2529403 * 52437899 * 254760293^2
  (every factor is checked prime by trial division and the product is checked equal to r - 1).
- 7 generates F_r^* (7^((r-1)/q) != 1 for every prime q above). Write zeta = 7.
- d = 2^20 divides r - 1. Put N = (r-1)/d (log2 N = 234.857089). 6 divides N. Put M = N/6
  (log2 M = 232.272127; M = 8334457902133654034844675144860897992084908882860761932151104402532352).
- r < 2^255.
- lambda = zeta^((r-1)/3) mod r is a primitive cube root of unity, beta is the cube root of unity in F_p with
  [lambda]P0 = (beta * x(P0), y(P0)); the program checks this equality on P0 by scalar multiplication. Since G1 is
  cyclic of prime order r and phi(x, y) = (beta x, y) is an endomorphism of the curve mapping G1 to itself,
  phi acts on G1 as [lambda].

Let E be the subgroup of F_r^* of order N (the d-th powers) and g = zeta^d, a generator of E. The sixth roots of
unity mu6 = {+-1, +-lambda, +-lambda^2} lie in E because 6 | N, and mu6 = <g^M>, so g^a is in mu6 exactly when
M divides a.

## 2. The algorithm

Notation. For a scalar 0 < e < r and a base point Q in {P0, B}, FB(Q, e) is the fixed-base product [e]Q computed
from tables (section 3). key(X) for a point X = (x, y) is min(x~, (beta x)~, (beta^2 x)~), where z~ is the
Montgomery representative z * 2^384 mod p, fully reduced to [0, p). A record is 64 bytes: the 6-word key and a
2-word index (index below 2^116) whose top bit is the list tag (0 baby, 1 giant).

Parameters (exact integers, printed by the program):
m = 57738922408142166625121974744311379 = ceil(sqrt(2M/5)) (log2 m = 115.475099),
G = ceil(2M/(5m)) = 57738922408142166625121974744311378 = m - 1,
m + G = 115477844816284333250243949488622757 (log2 = 116.475099), and m*G <= M.

```
Preprocessing (section 3): tables Tab[Q][k][c] = [c * 2^(85k)]Q for Q in {P0, B}, k in {0,1,2}, 1 <= c < 2^85.
Constants: kb = g * 2^256 mod r, kg = g^(-m) * 2^256 mod r, omega = zeta^N (order d), beta~.

Phase A (baby steps), in L = 2^16 lanes, lane t handling j in [t*ceil(m/L), (t+1)*ceil(m/L)) intersected with [0, m):
  e <- g^(first j of the lane) mod r
  for each j of the lane:  write record (key(FB(P0, e)), j, tag 0);  e <- MontMul(e, kb)        # e = g^(j+1)
Phase B (giant steps), same lane layout over i in [0, G):
  e <- g^(-m * first i) mod r
  for each i of the lane:  write record (key(FB(B, e)), i, tag 1);  e <- MontMul(e, kg)        # e = g^(-m(i+1))
Phase C (matching): LSD radix sort of all m + G records on the full 384-bit key (24 stable passes of 16 bits),
  then one scan comparing each record's key with the previous one.
  If no two adjacent keys are equal: output 1 (failure; all work above has been done and is charged).
  Else let the pair hold baby index j and giant index i.
Phase D (recovery of v = s^d): u0 = j + i*m. For k = 0..5: if FB(P0, g^(u0 + kM)) = B then u = u0 + kM mod N.
Phase E (recovery of s): s0 = zeta^u. For l = 0 .. d-1: if FB(P0, s0 * omega^l) = T, output s0 * omega^l.
```

All d candidates of phase E and all 6 tests of phase D are charged even when the answer appears early.

## 3. Fixed-base products and the absence of exceptional additions

Three windows of 85 bits cover every scalar below r < 2^255. For 0 < e < r write e = c0 + c1 2^85 + c2 2^170
with 0 <= ck < 2^85. FB(Q, e) starts with an empty accumulator; for k = 0, 1, 2, if ck = 0 the window is skipped,
if the accumulator is empty it becomes Tab[Q][k][ck] (a copy), otherwise one incomplete affine addition
acc + Tab[Q][k][ck] is made. So FB uses at most two additions, and the ledger charges two in every case.

Lemma 1 (no exceptional addition). Q has order r (P0 by the profile; B = [s^d]P0 with s^d != 0 mod r). When the
addition at window k happens, acc = [a]Q and the table point is [b]Q with a = sum_{i<k} ci 2^(85i) and
b = ck 2^(85k). Then 1 <= a < 2^(85k) <= b and a + b <= e < r. The incomplete formula fails only when the two
x-coordinates are equal, that is a = +-b mod r. a = b mod r is impossible because 0 < a < b < r; a = -b mod r
is impossible because 0 < a + b < r. Neither operand is the identity. The result is [e]Q, never the identity
because 0 < e < r. This is a proof for every scalar, not a probability; the toy run raises on any equal-x
addition and records zero such cases.

Every scalar used is nonzero mod r: baby scalars g^j and giant scalars g^(-im) lie in E, and the phase D and E
scalars are nonzero elements of F_r. So no step input is the identity, and every point a step produces is finite.

Table build. Each table row k of base Q has 2^85 - 1 entries. Qk = [2^(85k)]Q comes from 85k doublings (170
doublings per base). Entry 1 is Qk and entry 2 is a doubling. Entries 3 .. 2^85-1 are built in L = 2^16 chains, each
chain adding Qk to its previous entry; a chain's first entry [c]Qk is computed by left-to-right double-and-add.
Exceptions are excluded by the same argument: an addition [c]Qk + Qk with 2 <= c <= 2^85 - 2 would need c = +-1
mod r; in double-and-add, a doubling of [a]Qk with 0 < a < r is never exceptional (r is odd, so no point of G1 has
y = 0), and [2a]Qk + Qk with 2 <= 2a < 2^86 would need 2a = +-1 mod r. The batched chain additions use Montgomery's
batch inversion across the L chains; the chain heads, window bases and entry 2 are charged at the full price of an
unbatched affine operation (section 6).

## 4. Exactness of the matching

Lemma 2 (key). No point of G1 has x = 0 (the points (0, +-2) have order 3) or y = 0 (order 2). For X in G1 \ {O}
the six points +-phi^i(X) = [+-lambda^i]X are distinct and their x-coordinates form the 3-element set
{x, beta x, beta^2 x}. Two points X, X' of G1 \ {O} have equal keys exactly when these sets meet, that is when
x(X') = beta^i x(X), that is when X' = +-phi^i(X). Montgomery representation z -> z 2^384 mod p is a bijection and
is F_p-linear, so beta~ x~ (one Montgomery multiplication) is the representative of beta x, and
-(x~ + (beta x)~) mod p is the representative of beta^2 x (because 1 + beta + beta^2 = 0). So the key computed is
exactly a function of the class, and different classes give different keys. No hash and no truncation is used.

For e in E, [e]P0 and [e']P0 are in the same class exactly when e' / e is in mu6. Write v = s^d = g^u with
0 <= u < N. Baby j has scalar g^j; giant i has [g^(-im)]B = [v g^(-im)]P0 = [g^(u-im)]P0. They match exactly when
g^(u - im - j) is in mu6 = <g^M>, that is when u = im + j (mod M).

Lemma 3 (at most one pair, found exactly). Baby keys are pairwise distinct (j -> j mod M is injective on [0, m),
m <= M); giant keys are pairwise distinct (im with i < G are distinct mod M because Gm <= M). So equal keys occur
only as a baby-giant pair. The values im + j for i < G, j < m are the distinct integers 0 .. mG - 1, all below M, so
at most one pair (i, j) satisfies u = im + j (mod M), and it exists exactly when (u mod M) < mG. After a full sort,
equal keys are adjacent, so the scan finds the pair exactly when it exists.

Phase D. From the pair, v = mu g^(u0) with mu in mu6 = {g^(kM)}, so exactly one k in 0..5 passes
FB(P0, g^(u0+kM)) = B, and then v = g^u. Phase E. s = zeta^t for some t, and s^d = v = zeta^(du) means
d t = d u (mod r - 1), that is t = u (mod N), so s = zeta^u omega^l for exactly one l in 0 .. d-1, and phase E
finds it with the attack-side check [s]P0 = T. The output satisfies 1 <= s < r.

## 5. Success probability: exact count

s -> s^d is a homomorphism from F_r^* onto E with every fibre of size d, so for s uniform in 1..r-1, v is uniform
on E and u is uniform on 0..N-1, and (u mod M) is uniform on 0..M-1. By Lemma 3 the attack succeeds exactly when
(u mod M) < mG. The number of secrets s it recovers is exactly d * 6 * m * G out of r - 1 = d * 6 * M, so

P(success) = m G / M = 0.400000000 (exact rational printed by the program; it exceeds 2/5 by 2.46e-36).

The algorithm uses no random coins, so this holds for every choice of coins, and it is at least 0.39. On failure
the algorithm has still done every step charged below; nothing is rerun.

## 6. Cost ledger (1 unit = 2000 instructions = one batched incomplete affine G1 addition)

Pricing. A batched incomplete affine addition is one reference unit; its formula, its share of the batch inversion
and its internal arithmetic are not charged again. Other work is counted in scalar AArch64 instructions and divided
by C = 2000. An F_p multiplication is the measured blst_fp_mul price, 383. An F_r multiplication (4 limbs) is also
charged 383: an itemised 4-limb CIOS Montgomery multiplication is 174 instructions (4 rounds of 8 mul/umulh for
a_i*b, 10 add/adc, 1 mul for the quotient digit, 8 mul/umulh for q*r, 10 add/adc = 148; 13 loads; 9 for the final
conditional subtraction; 4 stores), so 383 is an over-charge. All F_p values are kept fully reduced in [0, p).

Large stores. The tables and the record arrays exceed 2^64 bytes, so their addresses and the sort's bucket
pointers and counters are two 64-bit words. Every word read or written in these stores is charged 2 instructions
(one to form the high address word, one access), and every pointer or counter update is charged as two-word
arithmetic. The machine model gives memory no scalar weight; this is the stated model for transfers.

Per step (baby or giant), worst case, no early exit: 2 additions = 2 units, plus these instructions:

| Item | Instructions |
| --- | ---: |
| Scalar update: one F_r Montgomery multiplication by the stored constant kb or kg | 383 |
| Split the 4-word scalar into three 85-bit two-word digits, 8 per digit | 24 |
| Zero test per digit (orr, compare, branch), 3 per digit | 9 |
| Table address per digit: two-word index times 96 plus window base, with carries, 10 per digit | 30 |
| Table reads: 3 entries x 12 words x 2 per large-store access | 72 |
| Operand staging into and out of the batch buffers, 60 per addition, 2 additions | 120 |
| beta x: one F_p Montgomery multiplication | 383 |
| beta^2 x = -(x + beta x) mod p: modular addition 19, modular negation 21 | 40 |
| Canonical reduction of x and beta x to [0, p), 13 each, charged even when already reduced | 26 |
| Minimum of three 6-word values: two compare-and-select of 12 | 24 |
| Record write: 8 words x 2, plus two-word output pointer increment 2 | 18 |
| Loop control: two-word counter increment 2, two-word compare 2, branch 1 | 5 |
| Step total | 1134 |

A table entry is a packed 96-byte affine point (x and y, 6 words each, no padding), stored and read as 12 words;
the storage figure in section 7 uses the same 96 bytes.

Sort and scan, per record. 24 LSD counting-sort passes of 16-bit digits cover the 384-bit key (key < p < 2^381).
Each pass is a histogram sweep and a stable scatter sweep:

| Item, per record per pass | Instructions |
| --- | ---: |
| Histogram: read the key word holding the digit 2, extract digit 2, counter address 1, two-word counter load/add/store 6, loop 5 | 16 |
| Scatter: read the key word 2, extract digit 2, bucket pointer address 1, two-word pointer load 2, copy the 8-word record (8 reads, 8 writes at 2) 32, two-word pointer increment 2, two-word pointer store 2, loop 5 | 48 |
| Per pass | 64 |

24 passes give 1536 per record. The scan: read the 6-word key 12, compare with the previous key (6 eor, 5 orr,
1 branch) 12, read and test the tag word 4, keep the previous record pointer 2, loop 5, total 35. Sort and scan:
1571 instructions per record. Per pass, the 2^16 buckets cost 8 instructions each (clear, two-word prefix sum,
store, loop): 24 * 2^16 * 8 instructions in all.

Per record (one step plus its share of the sort and scan): 2 + (1134 + 1571)/2000 = 3.3525 units.

Batch inversion. Each round of the L = 2^16 lanes has two addition stages, each with one F_p inversion by Fermat
(380 squarings and 228 multiplications at 383, plus 50) = 232,914 instructions, charged separately even though the
unit already amortises it. Rounds: ceil(m/L) + ceil(G/L) = 2 * 881026037721895853044463725957.

| Item | log2 units |
| --- | ---: |
| Steps: (m + G) * (2 + 1134/2000) | 117.835183 |
| Sort and scan: (m + G) * 1571/2000 plus bucket work | 116.126783 |
| Batch inversions: 2 * rounds * 232,914 / 2000 | 108.338753 |
| Table build (preprocessing, below) | 87.650917 |
| Lane starts: 2^17 exponentiations in F_r at 2*256*383, plus 64 | 23.615887 |
| Phases D and E: 6 + 2^20 fixed-base candidates at 2 units + (1134 + 25)/2000, the 6 with one extra F_r multiplication, plus three F_r exponentiations and 70 instructions to read the pair, form u0 and write the output | 21.367257 |
| Constants and reading B: 8 F_r exponentiations, 12 loads, 2 F_p multiplications | 9.616132 |
| Total W(2^20) | 118.221865 |

Arithmetic check: (m + G) * 3.3525 is 2^(116.475099 + 1.745237) = 2^118.220337; the batch inversions add
2^108.338753, a factor 1 + 2^-9.88, giving 2^118.221865; the remaining rows add less than 2^88 units, below
2^-30 of the total. The program sums the exact rationals and prints log2 W = 118.221865. The submitted bound is
log2 W(2^20) <= 118.23 = 11823/100, rounded up from 118.221865.

Preprocessing (included above). 6 tables of 2^85 - 1 entries: each entry 1 unit for the batched addition plus 90
instructions (60 staging, 24 for the 12-word write at 2, 6 address and loop); per table, 2^16 chain heads at
170 unbatched operations of 235,312 instructions each (inversion 232,914 plus 6 multiplications plus 100), one
inversion per round of chains, and entry 2; plus 340 unbatched doublings for the window bases of the two bases.
Total 2^87.650917 units.

## 7. Memory, advice and preprocessing

Peak storage: the tables, 6 * (2^85 - 1) * 96 bytes = 2^94.169925 bytes; the record array and its sort buffer,
2 * (m + G) * 64 bytes = 2^123.475099 bytes; lane state, counters (2^16 two-word counters) and code below 2^30
bytes. Peak at most 2^123.475099 + 2^94.17 + 2^30 bytes, below 2^123.48 bytes (claimed memory_log2_bytes 123.48).
Nonuniform advice: none (0 bytes); every table is built inside the charged work from P0 and B. Preprocessing
2^87.650917 units, inside W.

## 8. Assumptions, limitations and what is not claimed

- No heuristic is used for success or work: success is the exact fraction m G / M, and the work is the same for
  every challenge, computed for the worst case (all steps, all sort passes, all candidates in phases D and E).
- The instruction counts are analytic itemisations, not a built and measured implementation. The two
  multiplication prices come from COST_MODEL.md (383 measured for blst_fp_mul); the F_r multiplication is charged at
  that price, above its 174-instruction itemisation.
- The memory is astronomically large (about 2^123.5 bytes). The track gives memory no scalar weight and counts a
  load or store as an instruction; this submission charges 2 instructions per word access and two-word pointer
  arithmetic to address it. Latency, bandwidth and physical feasibility of such storage are not modelled.
- The batched-addition unit is used as the track defines it: 2^16 independent lanes advance together, so every
  addition stage is one batch; lanes whose digit is zero copy instead of adding and are charged the full addition.
- The bound is claimed for 2^20 <= n < r; only A_(2^20) and T are read.
- The toy run (part 2 of `experiment.py`) is evidence at its scope only: it uses a 12-bit subgroup (rt = 2689,
  d = 8, M = 56), plain integers instead of Montgomery representatives, one lane, windows of
  ceil(bitlength(rt)/3) = 4 bits instead of 85, radix passes of 4 bits over the toy key width (12 bits) instead of
  16 bits over 384, and rebuilds every table per secret. The order of operations is otherwise the pseudocode of
  section 2. It runs every secret, and the recovered count equals d * 6 * m * G exactly.

## 9. Experiment

`experiment.py` reads JSON {curve, target_permille} from standard input (initial input {"curve": 0,
"target_permille": 400}). Part 1 prints the production numbers above. Part 2 runs phases A to E on the chosen toy
curve for every secret s in 1..rt-1 with the program's own tables, sort and scan; it checks the key lemma on every
point of the toy subgroup, raises on any exceptional addition, asserts every output satisfies [s]P0 = T, and
compares the number recovered with d * 6 * m * G. At the initial input: m = G = 5, M = 56, 1200 of 2688 secrets
recovered, equal to 8 * 6 * 5 * 5 and to the fraction 25/56 of all secrets, with zero exceptional additions.
