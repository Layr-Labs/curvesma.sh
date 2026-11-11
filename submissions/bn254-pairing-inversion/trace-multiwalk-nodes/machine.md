# Scalar AArch64 programs and executed instruction counts

All work is counted in scalar AArch64 instructions on 64-bit general-purpose registers. Each emitted
instruction is one architectural instruction (LDP, STP, BL, RET, LDAR, STLR, STADD, LDAXR and WFE count one
each). `progs.py` generates every program; `experiment.py` executes them on an interpreter that models
64-bit registers, the carry and zero flags and real branches, compares their effects with Python, and prints
the counts. Each count below is an executed count. BL writes x30; every RET must return to the address in
x30, and every routine must return with x30 and sp as it found them (the interpreter asserts both).

**Encodability.** The generator asserts every immediate: ADD/ADDS/SUBS immediate 0 to 4095; LDP/STP
-512 to 504 in multiples of 8; LDR/STR 0 to 32760 in multiples of 8; MOVZ 16 bits with LSL 0, 16, 32 or 48;
SP adjustments 0 to 4095 x 4096 in multiples of 4096; shift amounts and EXTR lsb 0 to 63; AND bitmask
widths 1 to 63; branch and call targets are labels; B.cond, CSEL and CSET take only EQ, NE, CS or CC (the
conditions the interpreter executes); x18 is never used. Mode `algebra` probes each
guard with out-of-range values (`isa_negative`).

## 1. Per-record programs and units (mode `kernels`, production shapes)

| Program | Instructions |
| --- | ---: |
| table step | 379 |
| emit | 77 |
| walk step (two calls, table step, emit, counter, branch) | 460 |
| relay record | 19 |
| home record | 15 |
| radix count, pass 1 | 9 |
| radix scatter, pass 1 | 15 |
| radix count, pass 2 | 10 |
| radix scatter, pass 2 | 16 |
| radix count, pass 3 | 9 |
| radix scatter, pass 3 | 15 |
| scan step | 8 |
| **per record** | **576** |

Per record: 460 + 19 + 15 + (9 + 15) + (10 + 16) + (9 + 15) + 8.

| Unit | Instructions |
| --- | ---: |
| walk_init (12 ring traces) | 963 |
| walk_all, per full walk call (MOVZ, BL, SUBS, B.NE) | 4 |
| walk_all, fixed (entry, last call, exit) | 7 |
| callsend_w with send (full walker block) | 22 |
| relay message, fixed part | 31 |
| relay full-block extra (call, callsend_r, send) | 21 |
| relay END extra | 191 |
| relay flush, per home | 205 |
| home message, fixed part | 33 |
| home END extra | 187 |
| walker flush, per relay | 207 |
| walker flush, fixed | 6 |
| report send | 22 |
| report handler (last report) | 27 |
| report handler (other reports) | 25 |
| relay_up | 30 |
| output | 5 |
| wait (cap: two empty passes) | 11 |
| setup_send (15 setup messages) | 2,361 |
| h_setup | 142 |
| zero, per 4 counters | 5 |
| prefix, per counter | 6 |
| table build, per entry | 19 |
| scan exit, per call | 8 |
| fp_mul | 155 |
| fr_mul (the same listing modulo r) | 155 |
| Fp12 product | 13,843 |
| fp_mul calls in one Fp12 product | 54 |

## 2. What each program does

- **Table step** (`tstep`, 12 copies by ring phase). Traces are canonical residues in normal form. The
  accumulator starts at p - t_n. For i = 1 to 5 the pair sum t_{n+i} + t_{n+12-i} (below 2p, so below 2^255)
  and, for i = 6, t_{n+6} are cut into five 51-bit digits (LSR or EXTR, AND); each digit v indexes row v of
  table t = 5(i - 1) + d (MADD into the interleaved layout, two LDP) and the entry -c_i v 2^(51d) mod p is
  added with carries (ADDS/ADCS into 5 words). The 30 entries and p - t_n are each below p, so the
  accumulator stays below 31p < 2^259. Its top part v = acc >> 254 (below 32) selects v 2^254 mod p from the
  fold table; the result is below 2^254 + p < 3p, and two conditional subtractions make it canonical. It is
  then stored into the ring. Checked on 26 terms of two production walks (w = 0 and w = 15: start
  z^(w + 1) g^(beta_w), step D_w = g^-delta_w).
- **Table build** (`tabbuild`): entry (v, t) = entry (v - 1, t) - c mod p (SUBS/SBCS, conditional add of p
  by CSEL, STP), then the next row and the loop test. Checked against the step's entries for v up to 40.
- **Emit**: two splitmix64 hash chains over the 4 trace words, home (d1, d2), walk tag w 2^60, index (2-word
  ADDS/ADC), block and cursor addresses, 3 stores, cursor update; the K-th record tail-calls `callsend_w`.
- **Walk**: `bl tstep`, `bl emit`, SUBS, B.EQ per step; one B per 12 steps; 3 for entry and exit.
- **walk_all**: `walk` counts its steps in the one 64-bit register x22, so the walker's `walk_all` calls it
  q times with x22 = 3 x 2^60 (one MOVZ; a multiple of 12, so every call starts at ring phase 0) and once
  with the remaining 1 to 3 x 2^60 steps. Executed with q = 0, and with a 24-step chunk for q = 1 and 2,
  against the walk's traces.
- **send**: count and END stores, sequence increment, STLR of the sequence word (the completing store),
  STADD of the sent counter, RET: 11.
- **h_relay** (relay delivery processor 1): entry (5 STP, MOVZ of the fixed window base, LDAR of the
  sequence word, count, end pointer, state loads); per record 2 LDP and LDR from the window, d2 by LSL/LSR, block and cursor
  addresses, 3 stores, cursor update and full test, loop; on a full block BL `callsend_r`. END path: copy
  the 176 parameter words (z, every delta_w and beta_w; LDP/STP through two moving pointers that advance by
  512 bytes), count END messages, and after the 16th one flush every home: block address, count, parameter
  copy, END flag, `send`, cursor reset. Exit: overflow word, STADD of the handled counter,
  5 LDP, RET (x30 holds the handler's entry address, so the RET is the handler loop's branch back to
  its entry).
- **h_home** (home delivery processor): the same entry and exit; the home's one list cursor is loaded once
  per message and stored once; per record the cap test (CMP, CSET, ORR into the overflow flag; CSEL
  redirects the store to a scratch slot at the cap), slot address, 3 stores, cursor increment, loop. END
  path: parameter copy and the wake store. A planted overflow (cap 12 for 40 records) sets the flag and
  stores nothing past the cap; cap 64 does not.
- **report / h_rep / relay_up / output**: one-record messages (x and found flag); the handler keeps a
  reported x by CSEL, ORs the found flag, counts reports and writes the wake word after the last one.
- **wait** (on the node's wake word, alone in its 2^11-byte block): LDAXR, CBNZ, WFE, B per empty pass; LDAXR, CBNZ, RET on exit. With at most 3 passes: 11.
- **flush_w**: per relay: block address, count, the 176-word parameter copy, END flag, `send`, cursor reset.
- **setup_send / h_setup**: walker 0 copies the 128 coin words (every delta_w and beta_w) and the 2 salts
  into its setup block for walker w and sends it, for w = 1 to 15; the handler copies them into the
  walker's state and writes its wake word.
- **Radix passes** (`count`, `scatter`, three LSD digits of the 128-bit key: bits 0 to 42 by AND, bits 43 to
  85 by EXTR and AND, bits 86 to 127 by LSR), **zero** (2 STP per 4 counters) and **prefix** (exclusive
  running sum): sort 32 full-width random records, and in a second test with small digits the zero and
  prefix loops run in place of the Python prefix.
- **scan**: compares each record's 128-bit key with the next one's (EOR, EOR, ORR, CBZ to the equal exit),
  alternating two register pairs so no copy is needed, then advance and end test. A call costs at most 8
  per advance plus 8. Checked on a sorted list with runs of equal keys: the expanded pairs are exactly all
  equal-key pairs.

All transport programs run in a simulation of 36 nodes (16 walkers, 4 relays, 16 homes; mode `transport`)
with real BN254 traces; the executed total of every processor equals the cost functions of the ledger.

## 3. Field kernels (setup and rechecks)

- **fp_mul** (4-limb CIOS Montgomery, 155): inputs a, b below 2p are allowed. Each round adds a b_i
  (at most (2p-1)(2^64-1)) and m p (at most p(2^64-1)) to t and divides by 2^64. If t <= 3p-1 before the
  round, then (3p-1) + (2p-1)(2^64-1) + p(2^64-1) = (3p-1)2^64 gives t <= 3p-1 after it, below 2^256, so the
  5-word accumulator never overflows. At the end t < (4p^2 + 2^256 p)/2^256 < 2p since 4p < 2^256, and one
  conditional subtraction makes it canonical. Tested to 2p - 1.
- **fr_mul** (155): the same generated listing reading the modulus r and its Montgomery constant; the
  argument above holds with r in place of p (4r < 2^256). Tested to 2r - 1. The rechecks do their
  arithmetic mod r with it.
- **Fp12 product** (13,843): Karatsuba tower, 54 fp_mul calls, unreduced Fp2 sums; matches the schoolbook
  product in the pinned basis. An exponentiation below 2^254 is at most 508 products.
- **Output step** (`outstep.py`, executed in every mode on the x that reaches walker 0 in the simulation and
  on 1, r - 1 and a seeded x): Jacobian doubling for a = 0 (7 products, 14 sums) and mixed addition of P0
  (11 products, 14 sums), counted operation by operation; the affine result must equal an independent
  affine double-and-add, lie on y^2 = x^3 + 3, and its operation count times the `fp_mul` count must stay
  within the point-output bound of proof.md section 6. The affine reference also checks [r]P0 = infinity.
  The pairing check is priced from the executed Fp12 product and `fp_mul` counts; it is not executed.

## 4. Memory layout (bytes; every range below 2^63, pairwise disjoint, checked in every mode)

Each role's byte count includes 2^20 bytes for code, stacks, the state block and the wake block.

- **State block** (every node, at 2^33): fixed fields below 512, the parameter copy at 512 (1,408 bytes:
  z, then delta_0 to delta_15, then beta_0 to beta_15), port words at 2048, report fields at 2560, a home's
  list cursor at 3072 and scratch slot at 3200; the wake word at 2^33 + 4096, alone in its 2^11-byte block.
- **Walker:** ring at 2^26, constants at 2^32 (p, r and their Montgomery constants), state block, cursor
  table at 2^40 (2^36 words), stacks below 2^41 and 2^42, setup blocks at 2^45 (walker 0: one of 2,048 bytes
  per receiving walker), receive window at 2^46, outgoing blocks at 2^50 (2^36 blocks of 11,712 bytes: 64
  header, 256 x 40 records, 1,408 parameters), fold table at 2^57, the 30 tables at 2^58 (2^51 rows of 960
  bytes).
- **Relay:** state block, cursor table at 2^40 (2^36 words), three stacks, report block at 2^45, receive
  window at 2^46, outgoing blocks at 2^50 (2^36 blocks of 11,712 bytes).
- **Home:** state block, two stacks, report block at 2^45, receive window at 2^46, radix counters at 2^47
  (2^43 words), the list at 2^60 and the sort buffer at 2^61 (each 40 x cap bytes).

## 5. Coins

One 64-bit word of the algorithm's ordinary fresh independent coins costs 16 instructions: load, store,
loop bookkeeping. Walker 0 draws 32 values (delta_w and beta_w for each walk), each with at most 128
attempts of 4 words plus a 64-instruction range test, and 2 words for the salts: 524,320 in all, charged in
the setup row.
