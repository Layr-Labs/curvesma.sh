# Priced implementation and instruction bounds

Version 2 of e-man07's 73463878-0c19-4309-a991-749ea525d51c. The step
loop, its listing and the driver are unchanged; the meter-width section
and the 64-instruction store probe come from our gt-dlp port
76bccbcf-95f2-46da-bce2-481011ac5601. The window re-base is new: it is built from the
leaf kernels of e-man07's 7e58d165-ebf6-41aa-8003-2f0a42370b6a (copied
unchanged) and is priced by executing it, instruction by instruction.

All counts are executed scalar AArch64 instructions (64-bit general
registers; arithmetic, logic, loads, stores, address formation, branches,
calls and returns count one each), the unit of `tracks/docs/COST_MODEL.md`
(R in scalar AArch64 instruction equivalents, W = R/20,000). They are
counts, not timings. One processor.

## Evidence chain and units

Two kinds of listing are priced, and both count the same thing:

* The step loop `tt_run` is C (`ttstep.c`) compiled by Apple clang
  17.0.0, target arm64, `-O2 -fno-vectorize -fno-slp-vectorize`, and
  disassembled (`ttstep.disasm.txt`, `objdump -d --no-show-raw-insn`).
  The loop body has no internal branch or call, so one trip executes each
  of its instructions exactly once: the executed count is the static
  count, which `experiment.py` recounts inside the sandbox.
* The window re-base `rebase` and its leaves are AArch64 listings
  (`listings.txt`) produced by the generators in `kernels.py`.
  `experiment.py` regenerates them, asserts that they equal
  `listings.txt` byte for byte, checks that every line is an encodable
  A64 instruction (below), and executes `rebase` on a simulator of the
  instruction subset used (adapted from 7e58d165's, extended with `bl`,
  `ret` and the link register), counting every executed instruction,
  including each `bl` and `ret`. The listings are straight-line (no
  conditional branch), so the count is the same on every input
  (asserted) and equals the driver length plus the lengths of the called
  leaves (asserted).

So the 273 per step and the 11,730 per re-base are both numbers of
executed A64 instructions on 64-bit general registers; no SIMD or
floating-point register occurs in either (asserted for the clang listing;
the simulator and the encodability check accept only general registers).
The compiler only produces a concrete listing; each bound is a property
of its listing, which is itself a valid implementation.

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

Each entry is 2,048 bytes: 96 words F (offset 0), 8 child pointers (768),
48 words T (832), 8 words alpha, beta (1216) and 96 words of stored
Karatsuba sums of T (1280). `tt_run` reads only F and the child pointers
(offsets 0 to 831, through the entry pointer in x12; asserted by
`experiment.py`), so the fields appended after offset 1,280 in v2 leave
its listing unchanged. `ttstep.c` declares the first 1,280 bytes, which
is all that `tt_run` accesses.

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

## Window end: at most 11,923 instructions

Representation (the same as the tag kernel's input): an Fp coordinate is a
fully reduced 4-limb Montgomery residue y*2^256 mod p, least significant
limb first. Tower: Fp2 = Fp[u]/(u^2+1), Fp6 = Fp2[v]/(v^3 - xi) with
xi = 9+u, Fp12 = Fp6[w]/(w^2 - v). An Fp12 element is 384 bytes: the
Fp2 coefficient a_i of w^i sits at byte 64*(3*(i mod 2) + floor(i/2)),
i.e. X0 = (a_0, a_2, a_4) then X1 = (a_1, a_3, a_5), each Fp2 as (re, im).
This is the coordinate order of `reference.py` (`BN254GT.coords`) and of
the 48 limbs read by `tt_run`. Walk state: y (384 bytes), then a and b
mod r (32 bytes each).

Leaves (copied unchanged from 7e58d165; all inputs and outputs fully
reduced, so they compose freely):

| leaf | listing length = executed | what it computes | range argument (7e58d165) |
|---|---:|---|---|
| f2mul | 443 | (a0+a1 u)(b0+b1 u) in Montgomery form: T0=a0 b0, T1=a1 b1, T2=(a0+a1)(b0+b1) unreduced (8 limbs); c1=REDC(T2-T0-T1); c0=REDC(T0-T1, plus p*2^256 if negative) | a0+a1, b0+b1 < 2p fit 4 limbs; T2 < 4p^2 < p*2^256; T2-T0-T1 = a0 b1 + a1 b0 in [0,2p^2); T0-T1+[borrow]p*2^256 in [0,p*2^256); REDC of t < p*2^256 gives < 2p, one conditional subtraction gives the canonical residue |
| f2mulxi | 151 | (9+u)(a0+a1 u) = (9a0-a1) + (9a1+a0)u | three doublings and one addition, each followed by a conditional subtraction of p; one masked-borrow subtraction |
| fp_add | 22 | a+b mod p | a+b < 2^255, one conditional subtraction |
| fp_sub | 23 | a-b mod p | subtract, add p masked by the borrow |
| fr_add | 22 | a+b mod r | fp_add reading the constants of r |

The facts used: p, r < 2^254, so 2p < 2^255, 4p < 2^256 and
4p^2 < p*2^256. f2mul and f2mulxi save and restore the callee-saved
registers they use (x19..x28) in their own frame; fp_add, fp_sub and
fr_add write only x0..x17. Each leaf reads its constants (p or r, and
-m^-1 mod 2^64) with `adr` from a constant block placed next to the code.

The driver `rebase(out, in, entry)` (495 instructions, new) computes
out.y = in.y * T_e with the Karatsuba product of 7e58d165
(`Tower.mul_counted` with stored table sums): three Fp6 products
X0*Y0, X1*Y1 and (X0+X1)*(Y0+Y1), each of six f2mul; C1 = t2 - t0 - t1;
C0 = t0 + v*t1, with v*(c0,c1,c2) = (xi*c2, c0, c1). The operand Y = T_e
enters through Y0, Y1, Ys = Y0+Y1 and, for b in (Y0, Y1, Ys), the sums
(b1+b2, b0+b1, b0+b2); these twelve Fp2 values are stored in the entry
(offset 1,280, 768 bytes) when the table is built, so the window end only
loads them. The sums of the walk operand X are computed in the driver.
Then out.a = in.a + alpha_e and out.b = in.b + beta_e mod r. Every read of
in.y precedes the first write to out.y, so the re-base may run in place.
The driver keeps its three pointers in x19..x21, saves them and the link
register x30 in a 992-byte frame (temporaries for t0, t1, X0+X1 and six
Fp2 scratch values), and sets up each call with `add x0/x1/x2, base,
#offset` and `bl`.

| piece | calls | instructions |
|---|---:|---:|
| f2mul | 18 | 18 x 443 = 7,974 |
| f2mulxi | 7 | 7 x 151 = 1,057 |
| fp_add | 48 | 48 x 22 = 1,056 |
| fp_sub | 48 | 48 x 23 = 1,104 |
| fr_add | 2 | 2 x 22 = 44 |
| driver: prologue 6, 123 calls of 3 or 4 instructions (pointer set-up and `bl`) 485, epilogue 4 | | 495 |
| **`rebase`, executed** | | **11,730** |

One window end = `rebase` 11,730 + one `tt_run` entry/exit 65 + driver
128 = **11,923**. The driver allowance is the itemised bound of the next
section, unchanged from v1.

Encodability. `kernels.check_encodable` accepts only the forms used and
checks their immediates: ADD/SUB/ADDS/SUBS with imm12 in 0..4095
(register 31 means sp except as the destination of ADDS/SUBS, where it is
xzr); LDP/STP with a signed 7-bit offset scaled by 8, i.e. a multiple of
8 in [-512, 504], and two distinct registers for LDP; LDR (unsigned
offset) a multiple of 8 in [0, 32760]; register-form arithmetic, CSEL and
CSETM on general registers or xzr, never sp; `adr` and `bl` to a label of
the program. The largest offsets used are 2,016 (`add` into a table
entry) and 192 (`stp` in the f2mul frame). All 1,156 lines pass. Stack
frames (992, 208, 32 bytes) are multiples of 16.

Simulation checks (`experiment.py` section 5b). Leaves alone on 2,231
cases (edge values 0, 1, p-1, p>>1, 2^64-1, 2^192-1, random and all-limbs-
large inputs, output separate and aliased with operand a; fr_add on r):
0 failures, only the output and the leaf's own frame are written, sp and
x19..x28 restored. `rebase` on 32 inputs (16 edge pairs from the
identity, all-(p-1), zero and (p-1,0,...,0), 8 GT pairs, 8 arbitrary
Fp12 pairs; labels r-1 in four of them; in place and into a second
buffer): the 56 output words equal the Montgomery limbs of Y0*T_e,
computed in the protected basis by `reference.py`, followed by the label
sums mod r; only the output state and the stack are written; sp and
x19..x28 restored; every run executes 11,730 instructions; reading any
word that was not written first would stop the simulator. The tag of the
executed output certifies coordinate 0 of Y0*T_e on all 32 inputs.

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

`rebase` (11,730); convert the two coordinates of a_0 out of Montgomery
form by one f2mul with the integer pair (1,0) (447 with its call; checked
in section 5b), t = Re(a_0) - 9*Im(a_0) mod p (five fp_add/fp_sub);
s = floor(8t/p) by three four-limb comparisons with stored multiples of
p/8; the distinguished test t < ceil(p/2^80) by one four-limb comparison;
the aborted loop trip (< 273), one `tt_run` exit and re-entry (65),
driver work and the meter. Below 2^15 instructions, so at most 2^105 in
total.

## Fixed work, all inside F = 2^110

| item | bound |
|---|---:|
| parse the 12 x 64-hex input, basis change, Montgomery conversion, z^r=1 (<= 508 Fp12 products, each below 2^14: the product of `rebase` plus 24 fp_add calls for the second operand's sums) | < 2^25 |
| 16 multiplier scalars and 8 multipliers (two exponentiations each) | < 2^30 |
| table: 525,783,425,977,953 entries. Each: T_e = T_parent*M_s (one Fp12 product, < 2^14), the twelve stored Karatsuba sums of T_e (24 fp_add), the row v(T_e) (twelve small-constant combinations, < 2^14), 48 values u = v*2^(64i-256) mod p (Montgomery products by stored constants, two per f2mul call with operand (c,0)) and F = floor(2^128 u/p) (a 384-by-254-bit division yielding two quotient words, < 2^11 each), two label additions, eight child pointers (rank of the child multiset by the combinatorial number system: eight lookups in a binomial table and additions, < 2^9 each). Below 2^18 per entry; charged 2^20 | < 2^69 |
| walk starts: at most 2^50, each 508 fresh bits by rejection sampling (64 instructions per bit, at most 1,024 draws per scalar) and two exponentiations (<= 1,017 Fp12 products): < 2^26 each | < 2^76 |
| exact evaluations (meter 2^90) | < 2^105 |
| store: at most 2^50 inserts and lookups in an open-addressing table of 2^52 slots of 96 bytes keyed by a multiplicative hash of t, at most 2^70 probes in total (meter), at most 64 instructions per probe including the meter update | < 2^76 |
| false matches and the final hit: at most 2^50 recoveries (inverse mod r by exponentiation, < 2^18) and verifications g^x=z (< 2^25) | < 2^76 |
| output Q=[x]Q0: 254 affine doublings and additions over Fp2 with one inversion each, serialization | < 2^30 |
| every cap/failure path, allocation, code loading | < 2^40 |
| **sum** | < 2^106 < **2^110** |

## Storage

Table 525,783,425,977,953 x 2,048 bytes = 2^59.901 bytes. Store 2^52 x 96
bytes = 2^58.585 bytes. Multipliers, constants, binomial table, code,
stacks: under 2^30 bytes. Peak below **2^60.4 bytes**; every address is
an ordinary 64-bit pointer below 2^63 (the child pointers are plain
addresses; entry index times 2,048 is formed only during the build).
Memory is reported, not scored.

## Total

    S = floor(sqrt(1409*r/1000)) + 1 = 175614732316003428141458596702221752171
    R_total <= 273*S + 11,923*(floor(S/256)+1) + 2^110
    W = R_total/20,000,   log2 W <= 121.077977

valid for every run: success, step cap, walk caps, meters, degenerate or
false matches. `experiment.py` recomputes R_total exactly and rounds
log2 W up at the sixth decimal.
