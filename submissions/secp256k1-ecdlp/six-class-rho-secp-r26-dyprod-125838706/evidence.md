# Evidence (proof.md section 2.6)

`experiment.py` with `"kernels": true` (the declared initial run) does all of the following:

1. **Limb ports.** Python ports of the generated fe_mul (the six chains, as data), the reduction and its fused
   form, fe_sqr and fe_sub assert every no-overflow step above (each carry-only ADC, the exact product identity,
   s4 ≤ C0, w1 ≤ 1, s4' ∈ [−2, C0 + 2]) and the predicted value of every rare event of sections 2.2 and 2.4,
   against big integers on all pairs of 44 edge values and on random values, with 4,420 triples through the three
   fused forms (a^2 − b − c, a·b − c, a·b + c) and 4,135 through x3's a^2 − b − (c + 2C0) + 2^257 (40 with
   c = p − C0 − 1), 200 constructed carries and 50 constructed negative folds, and 1,504 products on shared rows
   (section 2.5), every row path taken: a shift and a fold that carry, and each that wraps 2^256. A 24-walk batch
   built from the ports (prod's dy_i, e_i and acc_i, pass 2's λ = R·e_i) is checked against −(W + T) by affine addition.
2. **The shipped loops, executed**: pass1, prod and pass2 in sequence on production points in batches of 1, 3, 6
   and 129 walks (3 and 6 run only the tail loops), every output against −(W + T) by affine arithmetic; the
   executed counts from 1 to 129 walks must differ by exactly 128 times 13.515625 (pass1), 171.515625 (prod) and
   339.515625 (pass2), with no exit taken: 139 walks, 0 mismatches.
3. **The emitters, executed.** `kernels/secp_tests.s` (the same emitters as functions) on edge and random inputs: 730 products, 444 squares, 818 fe_sub pairs, 1,670 fused a^2 − b − c
   (40 constructed negative folds) and 1,670 fused a·b − c, and 402 reductions, at least 200 of which carry.
   Every result must equal the port's value under the same rule. The loops' x3 and y3' forms have no
   wrapper; items 1, 2 and 6 check them.
4. **Counts and ledger.** It parses the counts from the shipped files (stubs after RET are in no loop), checks
   1,730 = 128·(7 + 4 + 2) + 64 + 2, 21,954 = 128·(7 + 2·72 + 3·4 + 6 + 2) + 64 + 2 and
   43,458 = 128·(2·72 + 3·4 + 73 + 89 + 7 + 10 + 4) + 64 + 2, and recomputes the doubling bound, the per-round
   work (fixed work: loop A 34 from item 7, prod 26, pass 2 26; pass1 9; the 127-walk tail charge 1,135.0625, 25,647 per round), the step
   550.546875 + 25,647/2^20 and the ledger of section 5.
5. **The arrival rows, executed.** `kernels/arrive.s` on production inputs in batches of 1 to 251 walks: random x and y, part A's near ties and near-zero values, x with
   x mod 2^192 > 2^192 − C0, keys forced distinguished (5%) and filter words forced to hit (15%). For every walk
   the stored table address and the rare-path code (logged by the file's test stub) must equal the bit-exact
   model `fast_j`, and every non-fallback index the exact rule. Rare-free batches of 64, 65, 128, 129, 192 and
   193 walks (walk 0 and the last pair on the odd path) must differ by exactly 31.015625 instructions per walk.
6. **Rare exits, executed.** pass1, prod and pass2 of `kernels/secp_loops.s` on 32 constructed batches of 12 to
   129 walks (129 reaches every copy of the groups, 12 and 20 the tail copy):
   - pass1's dx and prod's dy: borrows (a < b, (a − b) mod 2^64 < C0; for about half, limbs 1 and 2 of
     a − b + 2^256 are zero, so the borrow reaches limb 3), in every walk or the even ones; x19 must drop by exactly
     20 per predicted exit;
   - prod and pass2: rows whose three shifts carry (in one batch the first wraps 2^256) and operands a0 with
     a0·b in [T, T + a0), so both folds carry and wrap, cycled through the walks (prod's dy enters as y(T) with
     y(W) = 0);
   - pass2 with dx = 1, R = 1 and e chosen as λ: λ = 1 with x3 − x_W wrapping below C0 (also with y_W near p),
     λ = 2^128 with x3 − x_W = a (y3's carry), λ = 2^128 with x_W + x_T = p − m (x3's carry), and plain walks.

   Each kind is rotated through every walk position. Every output word must equal the limb ports of the inline
   kernels, x19 must drop by a multiple of 20, the exits counted by the meter must equal the stub entries, and all
   1,935 stubs must be entered: 2,227 walks, 2,600 exits, 0 mismatches against the limb ports.
7. **Loop A, executed** (`kernels/secp_fuse.s`, both copies, random words): rare-free batches of 1, 128 and 256
   give the fixed work 34, the group rates 43.03125 (pass1f) and 39.515625 (pass1p) and the tail walks 47 and 42;
   a 129-walk batch in which every walk takes the dx exit (about half of them reaching limb 3) enters all 258
   stubs. In mixed batches of 2, 3, 7 and 130 walks, walk i gets kind (i − n) mod 6: a distinguished key, a
   window, an odd-y borrow, a filter hit, a constructed dx exit, or none. Every address must equal the rule
   (on the rare path, the test stub's: the kernel's own j), every dx the big-integer value, the rare list
   the predicted walks, and the meter 20 per exit and 25 per rare entry. Each batch runs `pass1f` and `pass1p`
   (a filter hit is not rare in pass1p), odd y0 > p0 and K0 > p0 keys injected (54 odd borrows, entering no rare
   path): 656 walks, 53 and 28 rare entries, 155 exits each, 0 mismatches against the address rule of proof.md 3.2
   (BFI of bits 0..26 of y0, P′ = 0x3fffe17).
8. **The rule** (`rule_check`): 600 toy starts at R = 8, each walked under the rule (1, 256) (compare at every
   arrival, anchors every 256) and the production rule with the same phase and a cap of 2^15 operations: every
   start that ends distinguished under both, 599 of 600, ends at the same class and point after the same escapes.

The interpreter (`kernels/a64.py`) executes every opcode and addressing form of `secp_loops.s`, `secp_fuse.s`, `arrive.s` and `secp_tests.s` (`secp_fe.s` is counted, not executed), tracks the C
and Z flags and expands macros before running or counting, as the assembler does. Off the hosted runner, llvm-mc
(aarch64) assembles every shipped `.s` file, and for `secp_loops.s` and `secp_fuse.s` the assembled .text equals
that of the interpreter's macro expansion byte for byte; branch distances are inside their reaches (largest B.cond
175,140 bytes in `secp_loops.s`, 39,504 in `secp_fuse.s`). No AArch64 hardware or emulator ran them.
