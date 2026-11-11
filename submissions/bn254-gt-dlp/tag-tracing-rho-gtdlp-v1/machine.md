# Priced implementation and instruction bounds

Ported from e-man07's 73463878-0c19-4309-a991-749ea525d51c
(bn254-pairing-inversion). The listings and counts are unchanged; the
G2 output row of the fixed-work table is replaced by the scalar output,
the meter widths are stated, and the store probe price is raised to 64.

All counts are executed scalar AArch64 instructions (64-bit general
registers; arithmetic, logic, loads, stores, address formation, branches,
calls and returns count one each). They are static counts of supplied
listings, not timings. One processor; W = R/20,000.

## Evidence chain

* `ttstep.c` is the complete source of the step loop `tt_run`;
  `tower.c` is the complete source of the window re-base `rebase`
  (`fp12_mul` and two `fr_add`) with all its callees.
* `ttstep.disasm.txt` and `tower.disasm.txt` are their disassemblies
  (`objdump -d --no-show-raw-insn`) after compilation with Apple clang
  17.0.0, target arm64, flags `-O2 -fno-vectorize -fno-slp-vectorize`.
  The compiler only produces a concrete listing; the bound is a property
  of the listing, which is itself a valid implementation.
* `experiment.py` section 5 parses both listings inside the sandbox,
  asserts that no SIMD/FP register occurs, and recomputes every count
  below; it also asserts the ledger inequalities of the last section.

## The step loop: 273 instructions per step

`tt_run(cursor*, y)` receives the 48 Montgomery limbs of the window base
Y0 and a cursor (current entry e, a 64-bit trip budget). One trip:

1. A = sum over j<48 of UMULH(y_j, F[2j]) + y_j*F[2j+1] (mod 2^64): in the
   listing 48 `umulh`, 48 `madd` and the accumulating `add`s. The loop has
   121 word loads in all (the 96 table words, those limbs of y the
   compiler keeps on the stack, and the child pointer) and no store.
2. `add x25,x24,#0x60; and ...,#0x1fffffffffffffff; cmp ...,#0x61; b.lo`:
   the certificate test of proof.md section 2 (exit when not certified).
3. `lsr #61`, `add`, `ldr`: load the child pointer for s=A>>61; `cbz`:
   exit at a window end (null pointer).
4. `mov`, `subs`, `b.ne`: advance e, decrement the trip budget, loop.

The listing has exactly one backward branch (`b.ne` to the loop head); the
loop body it closes has **273** instructions, contains no call and no
store, and its only other branches (`b.lo`, `cbz`) are forward exits out
of the loop, not taken on a certified in-window step. A trip therefore
executes exactly 273 instructions. Function entry and exit (register
saves, loading y, writing the cursor back, restores, `ret`) are the
remaining **65** instructions of `tt_run`, executed once per call.

Each entry is 1,280 bytes: 96 words F, 8 child pointers, 48 words T, 8
words (alpha, beta). The kernel uses only F and the pointers.

## The driver: counters wider than 64 bits

The global step counter (initially S < 2^128) and the per-walk length
counter (initially 2^86) are two-word integers kept by the driver. Before
each call the driver sets the 64-bit trip budget to
min(steps left, walk length left, 2^62); after the call it subtracts the
number of trips made (budget before minus budget after) from both
two-word counters. `tt_run` returns
* REBASE when the entry has no child (window end): the driver calls
  `rebase` (Y0 <- Y0*T_e, labels updated), sets e to the root's child for
  s=A>>61 (A was certified), counts that step, and calls `tt_run` again;
* EXACT when A is not certified: the exact evaluation below;
* BUDGET when the trip budget is used up: the driver ends the run or the
  walk if a two-word counter is zero, otherwise refills the budget.
Every call ends within 256 trips (a window end), and the budget is reset
before each call, so a BUDGET return occurs only when a two-word counter
has reached zero, i.e. at most once per walk plus once per run. The
driver never calls `tt_run` with budget 0. A trip that leaves through an
exit executes fewer than 273 instructions and does not decrement the
budget: at a window end the driver then takes that step itself, so the
aborted trip is covered by that step's 273 in the 273*S term; at an
uncertified tag it is covered by the exact-evaluation allowance below.
All driver work per call (two-word subtractions and comparisons, the
min, pointer set-up, the call and the branch on the status) is below 128
instructions.

## Window end: at most 24,576 instructions

| function | executed | composition |
|---|---:|---|
| fp_add / fr_add / fp_sub | 35 / 34 / 35 | leaf, four limbs |
| fp_mul | 218 | leaf, coarsely integrated Montgomery product |
| fp2_mul | 873 | 3 fp_mul + 2 fp_add + 3 fp_sub + 44 |
| fp2_add / fp2_sub | 85 | 2 leaf + 15 |
| fp2_mul_xi | 399 | 9 fp_add + fp_sub + 49 |
| fp6_mul | 7,412 | 6 fp2_mul + 15 fp2 add/sub + 2 fp2_mul_xi + 101 |
| fp6_add / fp6_sub | 274 | 3 fp2 add/sub + 19 |
| fp12_mul | 24,044 | 3 fp6_mul + 4 fp6 add/sub + fp2_mul_xi + 3 fp2_add + 58 |
| rebase | 24,132 | fp12_mul + 2 fr_add + 20 |

All functions of `tower.c` are straight-line (no loop, no conditional
branch; asserted by the recount), so executed = static + callees.
`experiment.py` mirrors the Karatsuba tower formulas of `fp12_mul` in
Python and checks them against the protected polynomial basis. The
Montgomery product is the standard four-row coarsely integrated form: the
running value stays below 2p<2^255 after each row, so one conditional
subtraction returns the canonical residue. Inputs and outputs of `rebase`
are Montgomery residues, which is the form the tag kernel reads.

One window end = `rebase` 24,132 + one `tt_run` entry/exit 65 + driver
128 = 24,325 <= **24,576**.

## Meter widths

Every cap is enforced by a counter wide enough to hold it. The global step
counter (cap S < 2^127.05) and the per-walk length counter (2^86) are
two-word (128-bit) integers, as above. The exact-evaluation meter (cap
2^90) and the store-probe meter (cap 2^70) are also two-word counters,
initialised to their caps and decremented with `subs`/`sbc`; the run ends
as a failure when one reaches zero. The walk counter (cap 2^50) and the
64-bit trip budget (at most 2^62) are single words. No counter can wrap
before its cap. Updating a two-word meter takes under 8 instructions,
charged inside the exact-evaluation allowance and the per-probe price
below.

## Exact evaluation (metered, at most 2^90 times)

`rebase`; convert the two coordinates of a_0 out of Montgomery form (two
fp_mul by 1), t = Re(a_0) - 9*Im(a_0) mod p (five fp_add/fp_sub);
s = floor(8t/p) by three four-limb comparisons with stored multiples of
p/8; the distinguished test t < ceil(p/2^80) by one four-limb comparison;
the aborted loop trip (< 273), one `tt_run` exit and re-entry (65),
driver work and the meter. Below 2^15 instructions, so at most 2^105 in
total.

## Fixed work, all inside F = 2^110

| item | bound |
|---|---:|
| parse the 12 x 64-hex input, basis change, Montgomery conversion, z^r=1 (<= 508 fp12_mul) | < 2^25 |
| 16 multiplier scalars and 8 multipliers (two exponentiations each) | < 2^30 |
| table: 525,783,425,977,953 entries. Each: T_e = T_parent*M_s (one fp12_mul), the row v(T_e) (twelve small-constant combinations, < 2^14), 48 values u = v*2^(64i-256) mod p (one fp_mul by a stored constant each) and F = floor(2^128 u/p) (a 384-by-254-bit division yielding two quotient words, < 2^11 each), two label additions, eight child pointers (rank of the child multiset by the combinatorial number system: eight lookups in a binomial table and additions, < 2^9 each). Below 2^18 per entry; charged 2^20 | < 2^69 |
| walk starts: at most 2^50, each 508 fresh bits by rejection sampling (64 instructions per bit, at most 1,024 draws per scalar) and two exponentiations (<= 1,017 fp12_mul): < 2^26 each | < 2^76 |
| exact evaluations (meter 2^90) | < 2^105 |
| store: at most 2^50 inserts and lookups in an open-addressing table of 2^52 slots of 96 bytes keyed by a multiplicative hash of t, at most 2^70 probes in total (meter), at most 64 instructions per probe including the meter update | < 2^76 |
| false matches and the final hit: at most 2^50 recoveries (inverse mod r by exponentiation, < 2^18) and verifications g^x=z (< 2^25) | < 2^76 |
| output: x is already a canonical residue mod r; encode it as 64 lowercase hex characters (or 64 zeros when z=1) | < 2^12 |
| every cap/failure path, allocation, code loading | < 2^40 |
| **sum** | < 2^106 < **2^110** |

## Storage

Table 525,783,425,977,953 x 1,280 bytes = 2^59.223 bytes. Store 2^52 x 96
bytes = 2^58.585 bytes. Multipliers, constants, binomial table, code,
stacks: under 2^30 bytes. Peak below **2^60 bytes**; every address is an
ordinary 64-bit pointer below 2^63 (the child pointers are plain
addresses; entry index times 1,280 is formed only during the build).
Memory is reported, not scored. For orientation, the same algorithm with
a 16-step window (735,471 entries, 0.94 GB) costs 273+24,576/16 = 1,809
instructions per step, i.e. 2.3 bits more.

## Total

    S = floor(sqrt(1409*r/1000)) + 1 = 175614732316003428141458596702221752171
    R_total <= 273*S + 24,576*(floor(S/256)+1) + 2^110
    W = R_total/20,000,   log2 W <= 121.285447

valid for every run: success, step cap, walk caps, meters, degenerate or
false matches.
