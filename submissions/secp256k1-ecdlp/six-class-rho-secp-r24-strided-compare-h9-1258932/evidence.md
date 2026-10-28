# Evidence (proof.md section 2.6)

`experiment.py` with `"kernels": true` (the declared initial run) does all of the following:

1. **Limb ports.** Python ports of the generated fe_mul (the six chains, as data), the reduction and its fused
   form, fe_sqr and fe_sub assert every no-overflow step above (each carry-only ADC, the exact product identity,
   s4 ≤ C0, w1 ≤ 1, s4' ∈ [−2, C0 + 2]) and the predicted value of every rare event of sections 2.2 and 2.4,
   against big integers on all pairs of 44 edge values and on random values, with 4,420 triples through the three
   fused forms (a^2 − b − c, a·b − c, a·b + c) and 4,135 through x3's a^2 − b − (c + 2C0) + 2^257 (40 with
   c = p − C0 − 1), 200 constructed carries and 50 constructed negative folds. A 24-walk batch built from the
   ports is checked against −(W + T) by affine addition.
2. **The shipped loops, executed** on production points in batches of 1, 3, 6 and 129 walks (3 and 6 run the
   tails), every output against −(W + T) by affine arithmetic; the executed count from 1 to 129 walks must be
   exactly 99.515625 (pass 1) and 444.515625 (pass 2) per walk, with no exit taken.
3. **The emitters, executed.** `kernels/secp_tests.s` (the same emitters as functions) on edge and random inputs: 730 products, 444 squares, 818 fe_sub pairs, 1,670 fused a^2 − b − c
   (40 constructed negative folds) and 1,670 fused a·b − c, and 402 reductions, at least 200 of which carry.
   Every result must equal the port's value under the same rule. The loops' x3 and y3' forms have no
   wrapper; items 1, 2 and 6 check them.
4. **Counts and ledger.** It parses the counts from the shipped files (stubs after RET are in no loop), checks 12,738 = 128·(84 + 7 + 4 + 4) + 64 + 2 and
   56,898 = 128·(3·84 + 73 + 89 + 2·7 + 12 + 4) + 64 + 2, and recomputes the doubling bound, the per-round work with
   its tail walks and the ledger of section 5.
5. **The arrival rows, executed.** `kernels/arrive.s` on production inputs in batches of 1 to 251 walks: random x and y, part A's near ties and near-zero values, x with
   x mod 2^192 > 2^192 − C0, keys forced distinguished (5%) and filter words forced to hit (15%). For every walk
   the stored table address and the rare-path code (logged by the file's test stub) must equal the bit-exact
   model `fast_j`, and every non-fallback index the exact rule. Rare-free batches of 64, 65, 128, 129, 192 and
   193 walks (walk 0 and the last pair on the odd path) must differ by exactly 31.015625 instructions per walk.
6. **Rare exits, executed.** Pass 1 and pass 2 (pass 2 also alone, on chosen slots and Rinv) on 26 constructed
   batches of 12 to 129 walks (129 reaches every copy of the groups, 12 and 20 the tail copies): second borrows (a < b, (a − b) mod 2^64 < C0); products 2^128·a whose reduction carries
   (s4 = 1, s1 = 2^64 − 1), chained by walks that reset acc or R to 2^128; λ = 1 with x3 − x_W wrapping below C0
   (also with y_W near p); λ = 2^128 with x3 − x_W = a (y3's carry); and λ = 2^128 with x_W + x_T = p − m (x3's
   carry), each kind rotated through every walk position. Every output word must equal the limb ports of the
   inline kernels, x19 must drop by exactly 20 per exit, and all 1,167 stubs must be entered: 1,913 walks, 2,345
   exits, 0 mismatches.
7. **The fused pass 1, executed** (`kernels/secp_fuse.s`, random words): rare-free batches of 1, 129 and 130
   give the fixed work 80, the group rate 129.03125 and the tail walk 133. In mixed
   batches of 2, 3, 7 and 130 walks, walk i gets kind (i − n) mod 6: a distinguished key, a window, an odd-y
   borrow, a filter hit, a constructed fe_sub exit or none. So walk 0 gets a filter hit, a borrow, an fe_sub exit
   (its stub `.Lr1_w0`) and a window, the tail copy every kind (walks 1 to 6 of the 7-walk batch), and the group
   copies of the 130-walk batch every kind at rotated positions. Batches of 3, 12 and 131 force product carries
   (acc_(i−1) = 2^128) at walks 2; 2, 6, 11; and 2, 64, 130 (walk 0 forms no product). Every address must
   equal the rule (on the rare path, the test stub's: the kernel's own j), every dx and prefix product the
   big-integer value, the rare list the predicted walks, and the meter 20 per exit and 25 per rare entry. Each
   batch runs `pass1f` and `pass1p` (fixed work 78, 127.03125, tail 131; a filter hit is not rare there): 548 walks,
   79 and 54 rare entries, 40 exits each (one an fe_sub exit in walk 0), 0 mismatches.
8. **The rule** (`rule_check`): 600 toy starts at R = 8, each walked under the rule (1, 256) (compare at every
   arrival, anchors every 256) and the production rule with the same phase and a cap of 2^15 operations: every
   start that ends distinguished under both, 599 of 600, ends at the same class and point after the same escapes.

The interpreter (`kernels/a64.py`) executes every opcode and addressing form the
shipped files use, tracks the C and Z flags (all that EQ, NE, HS, LO, CS and CC read) and expands macros
before running or counting, as the assembler does.

Off the hosted runner (not shipped): `kernels/secp_fuse.s`, built with aarch64-linux-gnu-gcc and run under
qemu-aarch64 user mode (no AArch64 hardware), table 2^34 bytes aligned to 2^34: for each of `pass1f` and `pass1p`,
21,578 walks in 107 batches of 1 to 385 (775 fe_sub exits and 275 forced product exits taken; 3,222 and 2,412 rare
entries, `pass1p`'s 810 filter hits going straight) with 0 mismatches against big-integer dx, prefix products, the
25-bit address model and the meter. Each rare kind (distinguished, window, odd-y borrow, filter hit) and the
fe_sub exit is entered in all 130 walk copies of each, the rare path on both parity paths of every copy, and a
product exit in all 129 multiplying copies. `pass1f` assembles to the same words as before. Largest branch
distances (objdump, 43,456 instructions): TBZ/TBNZ 76 bytes (±32 KB); CBZ 79,892 (±1 MB); B.cond 86,988 (±1 MB);
B 87,008 (±128 MB). The assembler accepts every immediate.
