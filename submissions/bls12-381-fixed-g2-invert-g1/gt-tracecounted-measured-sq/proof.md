# Frobenius-class parallel rho with a cyclotomic-squaring branch, priced by its measured instruction count

**Claim.** For uniform nonidentity Q in G2 and uniform z in GT the algorithm below
recovers P in G1 with e(P,Q) = z with probability at least 0.39, using at most
**2^125.682** dense-Fp12-multiplication equivalents (C = 20,000), including setup, the two
pairings, table construction, walk starts, distinguished-point (DP) handling, cycle
escapes, the DP detection delay, output and verification. Peak memory is below 2^63
bytes, preprocessing (included) at most 2^44.08 units (the clearing pass over the report store, below
2^44.07 units, plus every other setup item, at most 2^31 units; section 5), nonuniform advice 0 bytes. The success
bound rests on the usual random-mapping heuristic for rho walks, adjusted for the branch
mix (section 4), and is supported by toy-field experiments and by running the production
step rule in the real GT. No production-size run is claimed.

The class construction (Frobenius classes of GT, the 4-multiplication key, the twisted
multiplier table that never applies a Frobenius to the walk element) is the one published
in submission 7da5f23a-6633-46b1-9c03-11f6cfedeeec (mpjunior92), which this package builds
on. The cyclotomic-squaring branch, the collision-law adjustment for the branch mix, the
explicit report and memory caps, and the reading of the branch, multiplier and DP bits from
the free class invariant T = Re(c_0) (so a squaring step needs no key computation) come from
later revisions of that design, the latest being the live record's package (claim 2^125.72).
Section 0 lists, item by item, every difference of this package from that record. The main ones are:
the squaring step is charged at the instruction count measured
for blst's own `blst_fp12_cyclotomic_sqr` (the measurement record ships in `measure/` and
`experiment.py` re-derives the price from it), the work meter and the multiword counters now have
explicit widths and priced updates (section 0-i), the report store is a fixed-depth array with
explicit addressing in one 64-bit address space (the distinguished-point rate is 2^-71 so that at most 2^55
reports of 80 bytes fit in 2^62.32 bytes), a stated load factor and a charged truncation rule instead of a
chained multimap (section 0-ii), and one restart cap is derived and quoted everywhere (section 0-iii).

## 0. What this revision changes

This package is a changed copy of the live record's package (`gt-tracecounted-exact-v25`, claim
2^125.72). The list below comes from comparing each file of the two packages (README.md, proof.md,
experiment.py, submission.json); it lists the input or rule that changed and, where a record value
is given, that value as written in the record's own files. Quantities that follow from a changed input
are given only for this package. Lines that only restate one of these changes are not listed again;
item 12 covers the small remainder.

The record charged its squaring routine 11,369 instructions, an analytic count: 9 Fp2 squarings (901
instructions each) and 34 Fp2 additions (90 each) from the measured blst_fp_mul = 383 and a
45-instruction modular addition, plus 200 for moves; every step also pays 270 instructions of
bookkeeping. This package charges the measured instruction count of the routine that performs the same
map, blst's `blst_fp12_cyclotomic_sqr`:

| Work | record's charge | this package |
| --- | ---: | ---: |
| every step: itemised bookkeeping (section 5) | 270 | **270** (same count, itemised listing) |
| every step: 192-bit work meter, exact compare with the cap on every step, and the 2-word stuck-walk counter (section 0-i) | not stated | **12** |
| every step, in all | 270 = 0.0135 units | **282 = 0.0141 units** |
| multiply step, keys only | 2,119 | **2,119** (same count, itemised listing) |
| multiply step in all | 1.11945 units | **22,401 = 1.12005 units** |
| cyclotomic squaring routine | 11,369 (analytic) | **9,803** (measured, section 5) |
| squaring step in all | 0.58195 units | **10,085 = 0.50425 units** |

**Changes from the record** (record value -> this package, where the record states one):

1. *Squaring price.* The squaring routine is the Callgrind-measured 9,803 instructions, not the
   analytic 11,369; the squaring step is 10,085 = 0.50425 units, not 0.58195 units.
   The package ships the measurement record (`measure/`, new) and the ledger reads it.
2. *Multiply-step price.* 1.11945 -> **1.12005** units (22,401 instructions). The
   2,119 instructions of key work are the record's; the extra 12 instructions are the every-step work
   of difference 3, which every step pays, a multiply step included (12 / 20,000 = 0.0006 units).
3. *Every-step work.* 270 -> 282 instructions: a 192-bit work meter updated by add/carry with the
   exact 3-word compare against the cap on every step (8 instructions), and a 2-word stuck-walk
   counter (4 instructions), both new and both with stated widths and charged updates (section 0-i),
   including 40 instructions per walk start or restart.
4. *Average step cost, ledger and claim.* c 0.9178875 -> **0.889125** units; W_nominal 2^125.7075 ->
   **2^125.6616283**; claim **125.72 -> 125.682**. This cap is the least thousandth covering 2 points of
   squaring-fraction drift (drift cost 0.6926% per point, ledger with that drift
   2^125.6814751), 0.0203 bits (1.422%) over the ledger, covering 2.053 points; the drift table of
   section 3.4 is recomputed, and the cap is checked by an exact integer inequality.
5. *Distinguished-point rate and what follows from it.* theta 2^-60 -> **2^-71**
   (DP bits 13 to 83 of T, fingerprint bits 84 to 211, not bits 13 to 72 and 73 to 200). Here the
   report cap is D = 2^55 (S theta = 2^54.83), the stuck-walk threshold 20/theta is 2^75.32, the detection
   window E = M 20/theta is below 2^95.33 steps, the expected false fingerprint candidates are 2^-19,
   and the restart cap is **2^36** with expected restarts 2^26.91 (counted per report gap) against the
   record's cap of 2^50 restarts, and the fresh-start work of the restarts is 2^46.1 units. The record's
   ledger had the term c M 20/theta; this one has c (M + 2^36) 20/theta (below 2^111.33 steps),
   so the restart allowance is charged as steps, and the non-walk bucket is 2^52.7 units (the record's
   bucket is 2^64.3).
6. *Report store.* The record kept a chained multimap at load factor 2/3 (peak memory below 2^73 bytes) and
   6,000 instructions per report. This package keeps
   a flat array of 2^52 buckets of 16 eighty-byte records with no pointer field, addressed in one 64-bit
   space (2^62.32 bytes, peak memory below 2^63 bytes), load factor 0.5 at the report cap, probe bound
   16, 4,000 instructions per report itemised (section 5), and a truncation rule: a report whose bucket
   holds 16 other records is dropped (section 0-ii). Before any walk starts the whole array is zeroed by an
   explicit clearing pass, charged in setup at below 2^44.07 units, and a stored record carries a used bit
   in a field whose true values never set it, which is what the free-slot test reads (section 0-ii). Section
   0-iv states the initial value of every counter, table and history, including a one-word restart counter
   charged below 2^23.30 units.
7. *Success chain.* The record's multimap lost nothing to the store. This package adds the
   store's charged truncation loss, at most 0.43 x 2q = 0.00709 with q = 0.00824 (the exact bucket-full
   tail P(Bin(2^55, 2^-52) >= 16) = 0.0082310), and its losses sum to 0.0170798, below 0.018: success
   after every loss at least 0.412, 0.022 over 0.39. The model value 0.43, the 0.01 charge for H1+,
   N, N', S = 2^125.8302 and eps = 6.5e-4 are the record's.
8. *Representation text.* Section 0a states the stored-limb rule directly (the record's text refers to a
   parent revision's rule on ordinary integers), section 3 reads T as the stored limbs, and a new
   paragraph of section 0a argues from the blst source that the squaring and multiply branches
   return the same stored limbs (an argument, not a measurement; section 6).
9. *Section 6 and 7 text.* The record's pessimistic-squaring sentence is replaced by a scope-of-prices
   item with sensitivities from part D; the record's sentence on an independent rerun of the script and
   its remark on the last control row are removed; the storage-size caveat and the part D description are
   new. The part-C readings quoted in section 3.4 are those of this package's declared run.
10. *experiment.py.* The code of parts A, B, C and E is byte-identical to the record's. The docstring
    and the `os` import changed, and part D is rewritten: it locates and reads `measure/`,
    asserts the three measured marginals and each raw record's total, recomputes three price vectors
    through one code path (asserting the first two), the integer cap inequality, the drift table, the
    store's load factor and truncation figure, the word widths and the per-report and meter charges, and
    reports the sensitivities.
11. *submission.json.* The claim is 125682/1000 (record 12572/100); the title and the check id are
    renamed; the evidence list has 5 entries (record 2: it adds `measure/README.md`,
    `measure/results.tsv` and `measure/disassembly-audit.txt`); the statement restates the bound, the
    squaring and every-step prices, the multiply step, c, W_nominal, the abort setting, the counter and restart
    description, the memory bound (2^73 -> 2^63) and the success floor 0.412, and no longer
    mentions the record's price-vector comparison; the check's purpose text adds what part D now reads.
12. *Small remainder.* The model rate 2^-13.6 per step is quoted as below 2^-13.58; the long-cycle
    allowance is restated at theta = 2^-71 (about 2^25.84 events, 2^101.2 steps in all); the
    squaring-fraction drift table lists 2, 2.5, 3, 4 and 5 points where the record's listed 1, 1.5, 2, 3
    and 4; and the claim line, the peak-memory bound and the section 7 description of part D are restated
    for the figures above.

`experiment.py` part D recomputes three price vectors through one code path and asserts the first two
before it prints this package's: the earliest rounded vector (11,400 / 300 / 2,200) gives c = 0.9225 and
W = 2^125.71472902377629, the record's vector (11,369 / 270 / 2,119) gives c = 0.9178875 and
W = 2^125.7074975, and this package's vector (9,803 / 282 / 2,119) gives the figures above. The two earlier
vectors are priced with part D's own bucket input 2^63.8 and window M 20 2^60 steps; the record's bucket
input is 2^64.3 (its experiment.py line 608), and the two differ by a term below 2^-60 of the ledger,
beyond the digits compared.

The measured count is of the stand-alone routine in a loop. Section 5 states exactly what it
covers, what the 282 instructions charge on top of it, and that the harness loop's own control
(about 5 instructions per call) is not deducted, so the squaring step is, if anything, over-priced.

### 0-i Counters and the work meter: widths and charged updates

Every counter of section 3.4 has a stated width and a priced update. The stuck-walk counter must
reach 20/theta = 20 * 2^71 = 2^75.32, which needs more than one word, so it is a
**2-word (128-bit) counter** in a register pair, incremented with `adds x, x, #1; adc y, y, xzr`
so that it cannot wrap. The report count must reach the cap D = 2^55 only: it is **one 64-bit
word**, incremented with `add x, x, #1` and compared with the cap by `cmp x10, x11; b.hs` (3 instructions, charged as 5 inside the report item of section 5), with no carry. The work meter must hold a whole run in instructions,
2^125.682 * 20,000 < 2^139.97 instructions, so it is a **3-word (192-bit) counter** in three
registers, with a headroom factor of at least 2^52.03 (192 - 139.97 bits).

Per step, both branches, the meter is charged **8 instructions** and the counter 4 in the table of section 5.
The step's charge (22,401 instructions for a multiply step, 10,085 for a squaring step: exactly the prices
the ledger of section 5 charges) is materialised in a register (`mov x4, #imm`, 1 instruction) and added
with `adds x0, x0, x4; adcs x1, x1, xzr; adc x2, x2, xzr` (3 instructions). The meter is then compared
with the stored 192-bit cap, which is held in three registers, by the full three-word compare
`cmp x0, x5; sbcs xzr, x1, x6; sbcs xzr, x2, x7; b.hs abort` (4 instructions). **This exact compare runs
on every step of every walk**, with no cheaper pre-test in front of it. The 2-word stuck-walk
counter (low word x8, high word x9) costs 4 instructions per step:
`adds x8, x8, #1; adc x9, x9, xzr; cmp x9, #2560; b.hs restart`. Only the high word is compared
because the threshold is an exact multiple of 2^64: 20/theta = 20 * 2^71 = 2560 * 2^64, so the counter
has reached it exactly when x9 >= 2560 (2560 fits the 12-bit immediate of `cmp`).

The stuck-walk counter is **reset to zero at every walk start, at every restart and at every
distinguished-point report**, so it holds the number of steps since the walk's last report (or since its
start), and the stuck-walk rule fires once per gap between consecutive reports (or from a start or restart to the first report), not once per walk. The
reset at a report is two moves into the register pair (`mov x8, xzr; mov x9, xzr`, 2 instructions); it is charged 4 instructions inside the
4,000-instruction report item of section 5 (a listed line of the item; the total stays 4,000),
at all 2^55 reports. At a walk start or restart the reset is part of the 20-instruction meter update
below, and the fresh start itself is charged in section 5 (at most 2^10.1 units each, 2^46.1 for the 2^36
restarts of the cap).

The meter is one 192-bit register triple for the whole run, not one per walk: every step of every walk
adds its charge to that same triple, so the abort test sees the total over all walks at every step. It is
shared by all M = 2^20 walks, so it is also updated outside the per-step path at
every walk start and every restart: 20 instructions per update (a full 3-word compare against the
cap, a 2-word counter reset, the register swap), 2 updates per walk (on entry and on hand-back).
The worst case is one walk start per walk plus the restart cap, so 2^20 * 40 < 2^25.33
instructions, below **2^11.04 units**, and 2^36 * 40 < 2^41.33 instructions, below **2^27.04 units** for the
restarts. Both are charged in the non-walk bucket of section 5.

### 0-ii The record store: one 64-bit address space, fixed-depth buckets, a charged truncation rule

Reports (at most 2^55, section 3.4) go into one array of **80-byte records with no pointer field**: 16
bytes of fingerprint F, 32 bytes for the p^kappa-scaled alpha, 32 for beta. The array is **2^52
buckets of 16 records**, 2^56 slots, **2^56 * 80 = 2^62.32 bytes**, which is the whole memory of
the store. It is one flat array in one address space; there is no sharding, no second address word and
no address-extension register. Every address quantity is below 2^63, so it fits one native 64-bit
register with the sign bit clear (`experiment.py` part D asserts each bound):

* bucket index i = the low 52 bits of F: `and x1, f_lo, #(2^52 - 1)` (i < 2^52);
* bucket byte base = array_base + i * 1,280 with 1,280 = 2^10 + 2^8: `add x2, x0, x1, lsl #10;
  add x2, x2, x1, lsl #8`. The offset i * 1,280 is below 2^52 * 1,280 = 2^62.32 < 2^63, so it fits the
  register without overflow;
* slot j of a bucket sits at bucket_base + 80 j, j = 0..15, reached by `add` with an immediate; the
  largest byte offset in the array is 2^62.32 - 1 < 2^63;
* record layout: bytes 0 to 15 hold F, bytes 16 to 47 the alpha field and bytes 48 to 79 the beta field,
  each 32-byte field as four little-endian 64-bit words, least significant first. Stored exponents are
  residues in [0, r) (each of the report's two multiplications by p^kappa mod r ends with its conditional
  subtraction, inside its 800 instructions), and r < 2^255 (`experiment.py` asserts r has 255 bits), so bit 255 of the alpha field,
  bit 63 of the 8-byte word at byte offset 40, is 0 in every true value. That bit is the **used bit**: storing
  a record sets it (`orr x, x, #0x8000000000000000`, 1 instruction, before the ten 8-byte stores), and a
  match masks it off when alpha is read back (`and`, 1 instruction, inside the at most 600 units of the
  candidate's verification);
* the **free-slot test reads exactly those 8 bytes at offset 40** of the slot and tests bit 63 (`ldr`, `tbz`):
  0 means free. The clearing pass below leaves it 0 in every slot, and every stored record has it 1, so no
  stored record reads as free and no free slot reads as used, with no probability to charge;
* slots of a bucket fill in order 0 to 15 and are never emptied, so the used slots of a bucket come before
  its first free slot. A probe of slot j is: offset add, load of the used word, `tbz` to the store if free,
  two 8-byte loads of the slot's F, `cmp`/`ccmp` with F, branch (8 instructions). The first used slot whose
  F equals F is the match; the first free slot ends the probe and receives the record.

**Clearing pass.** Before any walk starts, every byte of the array is set to zero: for each of the 2^52
buckets, 80 `stp xzr, xzr, [x0], #16` (16 bytes each, 1,280 bytes) and the loop's `cmp x0, x1; b.lo`,
**82 instructions per bucket**, plus 4 instructions to set the base and end pointers: 2^52 * 82 + 4 <
2^58.358 instructions, **below 2^44.07 units**. It is a named setup item of section 5 and is part of the
preprocessing bound (at most 2^44.08 units with every other setup item). No slot is assumed to start at zero
without this pass.

The report count is at most 2^55 < 2^64, a one-word counter, and its byte count 2^55 * 80 = 2^61.32 is
also below 2^63. The distinguished-point rate theta = 2^-71 is what makes this fit: the report cap scales with
theta, so at theta = 2^-60 (the record's rate) it would be 2^55 * 2^11 = 2^66 reports and the same load factor would need
2^73.32 bytes, which no 64-bit offset can address. The bits of T are disjoint (section 3): branch bits 0 to 2, table index bits 3 to 12,
distinguished-point bits 13 to 83 (71 bits), fingerprint bits 84 to 211; the low 52 bits of F are the
bucket index and F itself (128 bits) is the compared key.

The load factor at the report cap is 2^55 / 2^56 = **0.50000**, and the probe bound is **16,
fixed**. A report whose bucket already holds 16 records with other fingerprints is **truncated**
(dropped) after its 16th probe; nothing is ever chained, so the worst-case chain length is 16 and
the worst case (16 probes, 128 instructions) is what section 5 charges per report. A report
is dropped only if its bucket already holds 16 records with other fingerprints. After j earlier
reports (j <= 2^55, each landing in this bucket with probability 2^-52 under the uniformity that H1
assumes for F) the bucket holds at most Bin(j, 2^-52) of them, a count that grows with j, so every
report is dropped with probability at most the bucket-full tail at the cap,
q = P(Bin(2^55, 2^-52) >= 16). The mean is 2^55 / 2^52 = **8 records per bucket**, so the tail is that
of a Poisson(8) count (the 16 slots are not independent coin flips):
`experiment.py` part D computes it exactly (80-digit decimals) as q = 0.0082310, equal to the
Poisson(8) tail P(N >= 16) = 0.0082310 (the two agree within 1e-15; asserted), **quoted as q <= 0.00824**. The charge
uses this cap value for every report, an upper bound for each of them because the tail grows with j.

Effect on success. A useful collision needs two reports of one class to meet in the store. If the
first of them finds a full bucket it is dropped and so is the second (the bucket never empties); if
the first is stored the second finds it before any drop rule applies. Charging both reports, the
relative loss of success is at most 2q = 0.016463, **quoted as 0.01648**, an absolute loss of at most
0.43 * 2q = 0.0070787, **quoted as 0.00709**, so the 0.43 model value falls to at least
0.42292 (**quoted as 0.42291**). This is a charged loss, summed with the others in section 4 (total
below 0.018, success after every loss at least 0.412, a margin of at least 0.022 over 0.39). The bucket depth
is 16 and each stored report is 80 bytes, so the stride is 1,280 bytes and the probe charge of section 5 is for 16 probes.

Memory: the array 2^62.32 bytes, the multiplier table 7.1 MB, the exponent tables 0.8 MB and 2^20
walk states at 2 KB each (2^31 bytes) come to less than 2^62.33 bytes, so the claim states peak
memory below 2^63 bytes. Memory is reported, not scored (section 6).

Counters that exceed 64 bits are arithmetic counters, not addresses: the stuck-walk counter (it reaches
20/theta = 2^75.32, two words) and the work meter (2^139.97, three words), both as priced in section 0-i.
No address, index or report count needs more than one word.

### 0-iii One restart cap

The restart cap is **2^36 restarts**, derived once here and quoted with these digits everywhere
(sections 3.3, 3.4, 4 and 5, the README and the submission note). Derivation. A walk is restarted
when its stuck-walk counter reaches 20/theta = 2^75.32 steps, and the counter is reset at every report
(section 0-i). The rule therefore fires per **gap**, a run of steps between two consecutive reports (or
between a start or restart and the first report), and each gap is a separate chance. A gap ends at the
first distinguished point; each step is one at density theta = 2^-71, so a gap exceeds
20/theta steps with probability (1 - theta)^(20/theta) <= e^-20 < 2.07e-9. The gaps over the run are the reports of the
steps S(1 + eps) (the number of gaps is the number of reports, expected theta times the steps), so the expected number of
stuck-walk restarts over the whole run is S(1 + eps) theta e^-20 < **2^25.98** (the
experiment's 2^25.977; each walk has about 2^35 gaps, so this is not a per-walk figure). The other source of restarts is a walk
caught in a fruitless cycle of more than 64 steps (section 3.3): such cycles need all their steps to be
multiply steps (rate below 2^-100 per step), so the expected number of these events is S(1 + eps) 2^-100 = **2^25.84**.
The expected number of restarts of both kinds is **2^26.91** (the experiment's 2^26.906), which is 2^-9.09 of the cap 2^36.
No tail bound is taken on these counts, and the count of restarts is charged no success loss of its own: more than 2^36 restarts is a cap failure, which fails the run and is charged in the 0.01 of H1+ (section 3.4) like every other cap. The
expected counts are stated for the reader and the cap is far above them. The loss sum of section 4 does subtract two restart
terms, both for steps lost to restarts and not for their count: 4.2e-8 for the stuck-walk restarts and 3.8e-8 for the
long-cycle events (the two terms listed in section 4 and in experiment.py part D); no other restart term is subtracted. The meter supports the cap: 2^36
restarts cost 2^36 * 40 < 2^41.33 instructions of meter bookkeeping (below 2^27.04 units, section 0-i), 2^36 fresh starts of at most 2^10.1
units each (2^46.1 units, section 5) and, at worst, one full stuck-walk allowance each, 2^36 * 20/theta < 2^111.33 steps,
which section 5 charges in the step bucket and whose cost is 2^111.16 units, 0.0044% of the ledger. The cap
therefore neither breaks the meter nor adds to the success allowance, and it is the only restart cap in this document.

### 0-iv Initial values

Nothing the design reads is assumed to start at zero or at a known value; each is set and charged:

* the report store: zeroed by the clearing pass of section 0-ii before any walk starts (below 2^44.07 units);
* the 3-word work meter (set to zero), the 3-word cap register triple (loaded with W_cap), the one-word report
  counter, the one-word restart counter and the one-word false-candidate counter (each set to zero): at
  most 32 instructions once, before the walks, inside the at most 2^31 units of the other setup items (section 5);
* the restart counter: incremented and compared with the cap 2^36 at each restart (`add`, `cmp`, `b.hs`,
  3 instructions), 2^36 * 3 / 20,000 **below 2^23.30 units** at the cap, a named term of the non-walk sum
  of section 5; the false-candidate counter: the same 3 instructions per verified candidate, inside its at
  most 600 units;
* the stuck-walk counter: reset at every walk start, restart and report (section 0-i);
* the multiplier table T[h][t] and the scalar tables a_h p^t, b_h p^t mod r: computed in setup, the table item
  of section 5 (1,024 * (1,021 + 11) units plus 2 * 12,288 products mod r at 800 instructions each, below
  2^20.1 units; `experiment.py` asserts it);
* each walk's state y, alpha, beta: a fresh start g^alpha z^beta from fresh uniform (alpha, beta) at every
  walk start and restart (at most 2^10.1 units each, section 5);
* each walk's histories, the remembered previous T, the checkpoint T and the step count mod 64 of section 3.3:
  set at every walk start and restart, both T fields to the six-limb all-ones value 2^384 - 1, which no stored T
  (a residue below p < 2^381) equals, so the first steps cannot read a false repeat, and the count to zero.
  These few dozen instructions sit inside the fresh start's 2^10.1 units, whose two exponentiations take at
  most 1,021 units and leave more than 76 units (1.5 million instructions; asserted).

## 0a. Representation: which limbs the rule reads

A walk element is held as the six-limb Montgomery vectors that the pinned blst kernels store,
x~ = x*2^384 mod p. The decision rule is stated on those stored limbs, which is what a limb-level
implementation of the design reads at zero cost (reading an integer residue instead would need a
Montgomery reduction per T read and twelve per multiply step, work no charge here contains):

* T is the six limbs of Re(c_0) **as stored**. The branch bits are (stored T) mod 8, the table index is
  (stored T >> 3) mod R, and the distinguished-point predicate is read from the same word.
* kappa is the minimum of the twelve **stored** key products a*zeta^j and d*zeta^j, and the twist index
  is the smallest index attaining it.

This is the rule on ordinary integers composed with a fixed bijection of Z/p, so every structural
property carries over verbatim: the twelve stored keys are distinct; their multiset is the same for all
twelve conjugates (Frobenius permutes it), so kappa is class-determined and the twist index is covariant;
stored T is exactly invariant across conjugates because Re(c_0) is. What changes is which multiplier index
is selected, and that index is a pseudorandom function of the class under either rule, which is all H1 uses.
No conversion is performed anywhere in the charged path.

**The squaring branch returns the same stored form as the multiply branch.** blst keeps every Fp
value fully reduced: the double-width reduction tail `__redc_tail_mont_384`
(`build/elf/mul_mont_384-armv8.S`, lines 1336-1372) and the end of `__mul_mont_384` (same file, the
`subs`/`sbcs` block and the `csel ..., lo` group ending at line 818) both finish with one
conditional subtraction of the modulus, and `__add_mod_384` and `__sub_mod_384`
(`build/elf/add_mod_384-armv8.S`; the add falls through to `__add_mod_384_ab_are_loaded`) end with the
matching conditional correction (`csel ..., lo` after a trial subtraction, or an `and`-masked add-back), so every stored coordinate is the unique
residue in [0, p) times 2^384. An element of Fp12 therefore has exactly one stored limb vector, and
`blst_fp12_cyclotomic_sqr(a)` and `blst_fp12_mul(a, a)` (both return a^2 for a in the cyclotomic
subgroup, which contains GT) return identical limbs. **This is an argument from the blst source
(`measure/disassembly-audit.txt` section 5), not a measurement**: the harness image could not be
rebuilt for this revision, so a limb-by-limb comparison of the two routines' outputs was not run. The
shipped Python walk computes in ordinary integers and applies the Montgomery factor only to the values
it reads (part C, stored rule), which gives the same bits by the uniqueness just stated.

Evidence, part C, both rules, 2,000 production-field steps each from the same seed: the stored rule
gives 0 class-check failures in the first 150 steps, the exponent invariant y = g^alpha z^beta holds at
the end, squaring fraction 0.3677, chi-square of the low 3 bits 8.46 (7 dof) and of the next 4 bits 8.25
(15 dof); the integer rule gives 0 failures, 0.3864, 4.63 and 14.14. The selected index differs in
1,163 of 1,264 multiply steps, so the two rules are different walks and the claim is about the one the
charges can pay for. `experiment.py` reports both blocks (`partC_production_walk`,
`partC_production_walk_integer_rule`).

## 1. Reduction

Compute g = e(P0, Q) with the pinned convention pairing(Q, P0). Q has prime order r and
the pairing is nondegenerate, so g != 1 generates GT (prime order r). z = g^x for a unique
x in [0, r). Find x, output P = [x]P0, check e(P, Q) = z. If z = 1 output infinity. The
answer is for the supplied Q.

## 2. Arithmetic facts (checked in `experiment.py`, parts A and C)

Tower: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 1+u (the blst tower). With
u = w^6 - 1 this is the profile field, w^12 - 2w^6 + 2 = 0, and conversion of coefficients
is free: c_m = a_m - b_m, c_{m+6} = b_m.

**F1 Frobenius classes.** (sum c_m w^m)^p = sum conj(c_m) d_m w^m with d_m = xi^{m(p-1)/6}.
p has order exactly 12 modulo r, p^6 = -1 and 1 + p^4 + p^8 = 0 (mod r), so the maps
y -> y^{p^k} form a group of order 12 acting on GT, containing inversion. Every y != 1
lies in a class [y] of exactly 12 elements; N = (r-1)/12 classes, log2 N = 251.2721.

**F2 Keys.** Put zeta = xi^{(p^2-1)/6}; it lies in Fp and is a primitive 6th root of unity
(zeta^2 = zeta - 1). Let c_1 = a + b u be the w^1 coefficient of y and B = Re(conj(c_1) d_1)
= a Re(d_1) + b Im(d_1). The w^1 coefficient of y^{p^(k+1)} is conj(c^{(k)}) d_1, and
conj(d_1) d_1 = zeta, hence
`key_k := Re(w^1 coeff of y^{p^k}) = a zeta^j (k = 2j), B zeta^j (k = 2j+1)`.
With zeta^j in {1, zeta, zeta-1, -1, -zeta, 1-zeta} all 12 keys cost **4 Fp multiplications**
(a zeta, a Re d_1, b Im d_1, B zeta) plus 8 subtractions/negations and one addition. The
12 keys are distinct unless two of them coincide by a linear accident (probability <= 66/p
per element; ties are charged as failure). K([y]) = min_k key_k is a class invariant,
kappa(y) = argmin_k key_k, and y^{p^kappa} is the same element for all 12 members of the
class (the canonical representative y^). K is the minimum of 12 values that behave like
independent uniform residues (they come in negation pairs, so K lies in [0, p/2]); its low 200
bits are uniform up to statistical distance about 2^200/p = 2^-181 (the density of the minimum is smooth).

**F2b Trace invariant.** The w^0 coefficient of y^{p^k} is conj^k(c_0) (no twist), so T([y]) := Re(c_0) is the same
for all 12 conjugates, and 12 Re(c_0) is the trace of y to Fp. T is read directly from the coefficients (no
multiplication) and, for distinct classes, is an essentially uniform residue (collisions of T between distinct classes
have probability about 1/p per pair). It plays the role of the class identity for the bits below; K is needed only to
find kappa. `experiment.py` part A checks that Re(c_0) is equal over all 12 conjugates of random GT elements.

**F3 Cyclotomic squaring (Granger-Scott).** GT lies in the cyclotomic subgroup, where
squaring is cheaper than a product. Write y = A + B w + C w^2 with A=(c0,c3), B=(c1,c4),
C=(c2,c5) in Fp4 = Fp2[v]/(v^2-xi), v = w^3, and (A0,A1) = A^2, (B0,B1) = B^2,
(C0,C1) = C^2 computed in Fp4 with 3 Fp2 squarings each. Then
`c0' = 3A0 - 2c0, c3' = 3A1 + 2c3, c1' = 3 xi C1 + 2c1, c4' = 3C0 - 2c4, c2' = 3B0 - 2c2, c5' = 3B1 + 2c5.`
`experiment.py` verifies y^2 = formula on real BLS12-381 GT elements and counts the
operations of the implementation: **9 Fp2 squarings and 34 Fp2 additions/subtractions**
(the multiplications by xi are counted as one Fp2 addition).

blst's `blst_fp12_cyclotomic_sqr` (`src/fp12_tower.c` lines 630-660) is this map: three `sqr_fp4` calls
and the 3t - 2u / 3t + 2u combinations, on blst's Fp12 = Fp6[w]/(w^2 - v) layout, in which our coefficient
c_m of w^m is the Fp2 element a[m mod 2][m div 2]. Section 5 prices it by measurement.

**F4 Orders.** r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 *
2529403 * 52437899 * 254760293^2 (verified by multiplication in `experiment.py`). The order
of 2 modulo r is (r-1)/16 (about 2^250.9), so 2^t lies in the order-12 group <p> only if t is a
multiple of about 2^247.3.

## 3. Algorithm

Parameters: R = 1024 table multipliers, squaring probability s = 3/8, theta = 2^-71 (DP
density), M = 2^20 parallel walks. Bits of the class invariant T = Re(c_0) (read as the six stored Montgomery limbs of section 0a):
T mod 8 selects
the branch (squaring iff T mod 8 < 3), h = (T >> 3) mod 1024 selects the multiplier, the class
is distinguished iff (T >> 13) mod 2^71 = 0 (bits 13 to 83 of T). The fingerprint of a class used for matching reports is
F = (T >> 84) mod 2^128 (bits 84 to 211 of T), the 128 bits above the branch, multiplier and DP bits, so they are independent;
T has 381 bits, so every bit used is a bit of the stored T and none is shared between the roles. The bucket index of
section 0-ii is the low 52 bits of F.

Precompute M_h = g^{a_h} z^{b_h} for h < R (fresh uniform a_h, b_h mod r) and
T[h][t] = M_h^{p^t}, t = 0..11, with the scalar tables a_h p^t, b_h p^t mod r
(12,288 entries of 576 bytes, 7.1 MB).

A walk keeps some member y of the current class together with (alpha, beta) mod r such that
y = g^alpha z^beta. Each walk starts at fresh uniform (alpha, beta). One step:

```
1  T <- Re(c_0) of y (no multiplication); branch, h and DP bits from T
2  if class distinguished: compute the 12 keys (4 Fp mult) -> kappa; report (F, p^kappa*alpha, p^kappa*beta)
   into the record store of section 0-ii, keyed by F
3  2-cycle / checkpoint tests on T (section 3.3)
4  if T mod 8 < 3 (squaring branch):  y <- y^2 (Granger-Scott);  alpha,beta <- 2*alpha, 2*beta
   [no keys, no twist]
   else (multiply branch): compute the 12 keys (4 Fp mult) -> kappa; h <- (T>>3) mod 1024;
        t <- (12 - kappa) mod 12;  y <- y * T[h][t]   (ONE dense Fp12 multiplication);
        alpha += a_h p^t;  beta += b_h p^t
```

**3.1 The step is a map on classes.** Write y^ = y^{p^kappa}. Multiply branch:
y * M_h^{p^{-kappa}} = (y^ * M_h)^{p^{-kappa}} lies in the class [y^ M_h]. Squaring branch:
y^2 = (y^^2)^{p^{-kappa}} lies in [y^^2], because squaring commutes with Frobenius. Both
depend only on the class (T is a class function; kappa and y^ enter only through the conjugation
absorbed by the twist), so walks that reach the same class stay together and reach the same DP.

**3.2 Collision.** The store keeps every report (records are appended, never replaced), and its
truncation rule (section 0-ii) drops only a report whose bucket already holds 16 records with other
fingerprints, with probability at most q <= 0.00824 per report (the bucket-full tail at the cap, 8 records
per bucket on average); a useful class can be shadowed only by that rule, whose effect on success
(a charged loss of at most 0.00709) is priced in section 0-ii. Two reports with equal F from
different walks are the same class, up to
a false match of probability 2^-128 per pair of distinct classes: with D <= 2^55 reports there are about
D^2/2^129 = 2^-19 false candidates in expectation; every candidate is verified (one modular inverse and
g^x = z, <= 600 units), which is negligible. Their canonical exponents (alpha^, beta^) = p^kappa (alpha, beta) satisfy
alpha^1 + x beta^1 = alpha^2 + x beta^2 mod r, so x = (alpha^1 - alpha^2)/(beta^2 - beta^1).
beta^1 = beta^2 happens with probability 1/r. Output P = [x]P0 after checking g^x = z.

**3.3 Fruitless cycles.** A closed walk of the class graph is *useless* if the canonical
exponent pair returns unchanged. After l steps the pair is Lambda (alpha^_0, beta^_0) + tau
where Lambda = 2^t p^m (t squaring steps, m a sum of twists) and tau = sum_h c_h (a_h, b_h)
with c_h in Z[p]. Because (alpha^_0, beta^_0) is uniform and independent of the table, the pair
returns unchanged with probability about 1/r unless Lambda = 1 and every c_h = 0.
* **No squaring step can occur in a useless cycle.** Lambda = 1 with t >= 1 needs
  2^t in <p> (F4), i.e. t >= 2^247.
* A useless cycle made of multiply steps needs every table index used to occur at least twice
  with coefficients p^{e_1} + p^{e_2} + ... = 0 in Z[p]/(Phi_12(p)): for two occurrences
  e_2 = e_1 + 6 (probability 1/12), for three {e, e+4, e+8}, and so on; for any number of occurrences
  the probability is at most 1/12 (the last twist must hit one value). For 2-cycles this gives the
  model rate (5/8)^2/(12 R) <= 1/(12 R) < 2^-13.58 per step (measured 4.0e-4 to 5.6e-4 per step across
  seeds versus 5.1e-4 predicted at R = 64, section 7). For l = 2k indices drawn independently and
  uniformly (heuristic H2) the probability is at most (2k)!/(2^k k!) (12R)^-k (pairings dominate), which
  is below 2^-288 for 65 <= l <= 147. For l >= 148 all l consecutive steps must be multiply steps, which
  has probability (5/8)^l < 2^-100 (the branch bits are bits of T, H2). Such cycles are never met in
  practice; if one is, the stuck-walk rule restarts the walk after 20/theta = 2^75.32 steps, which adds
  less than 2^-14 relative work even at the restart cap of section 0-iii (2^-14.5, section 5).
* **Detection.** The step before last is remembered; if its T equals the current T the two
  classes form a 2-cycle and the walk replaces the representative of the class with the smaller
  T by its square (y -> y^2, exponents doubled), cost <= 3 units. This depends on the cycle
  only, so merged walks stay merged. A checkpoint T is stored every 64 steps and compared at
  every step; a repeat means a cycle of length <= 64, which is traversed once more to its minimum-T
  class that is then squared (waste <= 200 steps). Every other walk whose stuck-walk
  counter reaches 20/theta = 2^75.32 steps is restarted from a fresh uniform pair (stuck-walk rule). The counter is reset to zero
  at every walk start, restart and report (section 0-i), so a walk is restarted only after 20/theta steps with no report, and the
  rule is checked once per report gap.
* Total overhead factor on the step count: eps = 5/(12R) + 250/R^2 = 6.5e-4 (charged).

**3.4 Caps and the total-work abort.** The algorithm meters its own charged work with the prices of
section 5 (1.12005 units per multiply step, 0.50425 per squaring step, the stated costs of escapes,
starts, report handling and verifications) in the 3-word (192-bit) meter of section 0-i, and **aborts
with failure as soon as the metered work exceeds W_cap = 2^125.682**, the submitted bound. The claimed
bound therefore holds for every run by construction;
what remains is the probability that the abort fires, which is part of the failure budget below.

*Heuristic H1+ (no drift of the walk statistics).* The cap is sized to cover a drop of the squaring
fraction by 2 percentage points. A multiply step costs 1.12005 units and a squaring step 0.50425, a
difference of 0.6158 units, so a drop of D points raises the average step price by 0.006158 D units, a
factor 1 + delta with delta = 0.006158 D / 0.889125 = 0.6926% per point (2 points: 1.3852%). The nominal
ledger with every step priced that much higher is 2^125.6814751 (`experiment.py` part D), the least
thousandth of a bit at or above it is 125.682, and the exact integer inequality
(numerator^1000 <= denominator^1000 * 2^125682 for that ledger rounded up to 40 digits) is asserted there.
The cap therefore exceeds the nominal ledger W_nominal = 2^125.6616283 by at least 1.422% (0.0203 bits), a slack of
0.012644 units per step, and is not exceeded while the squaring fraction stays within 2.053 percentage
points of 3/8. The part-C walk of 2,000 steps has a binomial standard deviation of 1.08 points and reads
0.3677 in the declared run (0.73 points below 3/8; the integer rule reads 0.3864). A larger drop does not break the bound either: it makes the meter halt a
corresponding number of steps early, and the success probability of the shorter walk is still given by the
section-4 collision law. The meter then buys n/(1 + delta) steps of the nominal n:

| drop of the squaring fraction | success value | probability spent |
| ---: | ---: | ---: |
| up to 2 points (the cap's design point, 2.053 points) | 0.4304 | 0 |
| 2.5 points | 0.4284 | 0.0016 |
| 3 points | 0.4262 | 0.0038 |
| 4 points | 0.4219 | 0.0081 |
| 5 points | 0.4177 | 0.0123 |

H1+ states that over the whole run the frequencies of the branch bits, of the 2-cycle events and
of the DPs stay inside a drift of this size — the 0.01 budget is not exhausted until about 4.46
percentage points — and that the number of reports, false fingerprint candidates and restarts stays within
the caps below, with total probability of any cap or the abort firing at most 0.01. This is a much weaker
statement than independence of the steps; it excludes only a drift by percentage points over 2^125 steps.
Support: (i) the marginal probability of the squaring branch is exactly 3/8 for every class (T mod 8 is
uniform, F2b), so the expected charge of a step is exactly the nominal one with no independence;
(ii) in the toy runs the fluctuation of the mean charge of a run scales as the independent-steps
prediction: rms of (relative deviation x sqrt(steps)) is 0.300, 0.277 and 0.308 observed against 0.283
predicted (300, 150 and 200 runs in q = 67, 163, 211, part B; the declared run of this package
reads 0.3043 against 0.2848 when its 6-second budget completes all 40 trials; a slower machine completes fewer
trials and reads a nearby value), so there is no excess correlation of branch labels along the walk at that
scale; (iii) the production-field walk of part C has squaring fraction 0.3677 (stored rule) and 0.3864
(integer rule) in the declared run, both within the binomial scatter of 3/8. The toy runs cannot test
concentration at 2^125 steps; H1+ extrapolates it, exactly as H1 extrapolates the collision law. The
success analysis (section 4) charges 0.01 for it.

*Caps.* reports D <= 1.1 S theta + 2^25 < 2^54.969 <= 2^55 (S theta = 2^54.83; 2^25 covers the at most
20 M repeated reports of the merged trails), counted in a one-word counter (section 0-i); at most 2^16 false
fingerprint candidates are verified (the expectation is 2^-19, so Markov gives probability 2^-35; the
verification work is at most 2^16 * 600 units < 2^25.24); **at most 2^36 restarts** (section 0-iii: expected 2^25.98 stuck-walk restarts, counted per
report gap, and 2^25.84 long-cycle events, 2^26.91 in all, which is 2^-9.09 of the cap; a run that exceeded the cap is a cap failure inside the
0.01 of H1+, with no separate loss for the count (the step-loss terms 4.2e-8 and 3.8e-8 of section 4 are the only restart terms in the loss sum); each restart a fresh start plus 40
instructions of meter bookkeeping, 2^27.04 units in total, and at worst one stuck-walk allowance each,
2^36 * 20/theta < 2^111.33 steps, charged in the step bucket of section 5); the total number of steps
<= S(1+eps) + (M + 2^36) 20/theta. Repeat steps after the first
collision are at most E < 2^111.33 steps (<= 1.12005 units each), inside the cap slack, which is at least
2^119.5 units (1.422% of 2^125.6616283). Fruitless cycles
longer than 64 steps need all their steps to be multiply steps (rate below 2^-100 per step, about 2^25.84
events in expectation, each wasting at most 20/theta = 2^75.32 steps before the stuck-walk restart: 2^101.2 steps in total in expectation).
These caps bound the work and memory stated below.

## 4. Success probability

**H1 (random mapping with branch mix).** The class sequence behaves like iterates of a
random mapping on N classes in which two distinct classes have the same image with
probability (1 - sum_b q_b^2)/N, where q_b are the branch probabilities. Justification: for
c1 != c2 the equation F(c1) = F(c2) has, for each branch pair (b1,b2) with b1 != b2, about one
solution among the N classes (12 candidate representatives, each canonical with probability
1/12); for b1 = b2 = squaring it has none (squaring is a bijection of GT of odd order that
commutes with Frobenius), and for b1 = b2 = a table index it has about 11/12 of one. We
use the conservative value sum q^2 = s^2 + (1-s)^2/R = 0.141006 for all pairs; the exact
value (s^2 + (1-s)^2/(12R)) is smaller by 3.5e-4. Part E of `experiment.py` tests this law directly (section 7). This is the same correction as the Brent-Pollard
factor R/(R-1) of the pure adding walk (s = 0), and for Pollard's original three-branch walk
it is the mechanism that makes Pollard's original square-and-add walk slower than an r-adding walk
(Teske 1998). Hence the effective class count is
N' = N / (1 - sum q^2), log2 N' = 251.4914.

Under H1 the probability that the first S steps contain a class collision is
1 - exp(-S^2/(2N')). Choose P_model = 0.43:
`S = sqrt(2 N' ln(1/0.57)) = 2^125.8302.`
Losses before the first useful collision is detected: detection failure e^-20 < 2.1e-9, a
collision inside a single walk (fraction about 1/M = 1e-6), beta^1 = beta^2 (2^-255), cap
and abort failures (at most 0.01 together, H1+, section 3.4), steps lost to restarts (a stuck-walk restart
is charged a full allowance of 20/theta = 2^75.32 steps and fires with probability e^-20 per report gap of mean length 1/theta steps, so the expected
fraction of steps lost to it is 20 e^-20 < 4.2e-8, subtracted as 4.2e-8, and the expected fraction lost to the long-cycle events, each charged
one full allowance, is below 3.8e-8, subtracted as 3.8e-8; expected restarts 2^26.91, counted
per report gap, far below the cap of **2^36** of section 0-iii, each costing a
fresh start of at most 2^10.1 units, 2^46.1 in total for the starts at the cap, and the 40 instructions of meter
bookkeeping, inside the 2^52.7 non-walk bucket; the worst-case 2^111.33 restart steps are charged in
the step bucket of section 5) and the record store's truncation (at most 0.43 * 2q = 0.00709 of
success, with q = 0.00824 the bucket-full tail at the cap, section 0-ii). Their sum is the 0.01 of cap
and abort failures, plus 0.00709, plus about 1e-6 for the rest: 0.0170798, **below 0.018**
(experiment.py part D). The success after every loss is therefore at least 0.43 - 0.018 = 0.412,
leaving **0.022 absolute slack** over the required 0.39 for the inaccuracy of H1. Toy measurements (section 7) of the fraction of
runs with a useful collision within S(0.43) are 0.41-0.475 per setting, pooled 0.447 +- 0.011 over 2,200
runs, at 0.437 for S(0.42); individual settings of a few hundred runs scatter by +-0.04.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Instruction counts are executed scalar AArch64 instructions, counted as the cost model counts them
(`tracks/docs/COST_MODEL.md`): Valgrind Callgrind on Debian 13.7 with GCC 14.2.0, blst at the pinned
revision f7a88067, blst's C compiled with `-O3 -fno-builtin -fPIC -fno-tree-vectorize
-mgeneral-regs-only`, collection toggled around the call loop only, and the per-call figure taken as the
marginal (Ir(65 calls) - Ir(1 call)) / 64. The record is in `measure/` (README.md, results.tsv, six raw
Callgrind records, harness/, disassembly-audit.txt) and `experiment.py` part D reads it.

| routine | per-call marginal | status |
| --- | ---: | --- |
| `blst_fp_mul` | 383 | measured; equals the cost model's reference count |
| `blst_fp12_mul` | 22,996 | measured; equals the cost model's reference count |
| `blst_fp12_cyclotomic_sqr` | 9,803 | measured, on real GT elements, identical at five slice offsets |

The first two lines reproduce the cost model's own two figures to the instruction, which is the
calibration of the instrument; the third is the price used below. The squaring routine is
straight-line (18 functions, 1,017 static instructions, no conditional, backward or indirect branch
besides the 16 `ret`s; `measure/disassembly-audit.txt`), so its count cannot depend on the element.

**Table multiply:** the product y*T[h][t] is the dense blst_fp12_mul reference operation,
charged exactly 1 (22,996 instructions, rounded to C = 20,000).

**Squaring branch:** the Granger-Scott map of F3 is what blst's `blst_fp12_cyclotomic_sqr`
(`src/fp12_tower.c` lines 630-660) computes. Its measured marginal, **9,803 instructions**, is charged
for it in full. A squaring step needs no key computation. The harness loop's own control (about 5 of the
9,803 instructions: increment, compare, branch, call) is not deducted, although the 270 below also
charges loop control, so the squaring step carries that control twice. For comparison, the analytic
count of the same formula was 9 Fp2 squarings x 901 + 34 Fp2 additions x 90 + 200 moves = 11,369
instructions, 1,566 more than the measured routine; the routine's composition is in
`measure/disassembly-audit.txt` section 2.

**Per-step work.** Every step (both branches), analytic, as itemised before:

| Work | Instructions |
| --- | ---: |
| read Re(c_0) (6 limbs), branch bits, DP test, table address | 50 |
| (alpha, beta) update mod r (two 4-limb add + conditional subtract, or doubling) | 100 |
| 2-cycle and checkpoint comparisons on T, history | 60 |
| loop control | 30 |
| margin for loads and moves | 30 |
| **subtotal, the itemised listing of the record** | **270** |
| 192-bit work meter: step charge in a register (1), adds/adcs/adc (3), exact 3-word compare with the cap, cmp/sbcs/sbcs/b.hs, on every step (4) (section 0-i) | 8 |
| 2-word stuck-walk counter: adds/adc increment (2), cmp of the high word with 2560 and b.hs (2), exact as 20/theta = 2560 * 2^64 (section 0-i) | 4 |
| **total counted 282, charged** | **282 = 0.0141 units** |

Multiply steps in addition (the keys give kappa), analytic, built on the measured blst_fp_mul = 383 and
a modular add/sub/neg of 6 limbs called as a routine at 45 instructions (6 add/sub, 6 conditional-subtract,
6 select, 12 loads, 6 stores, 9 for call, return and address set-up):

| Work | Instructions |
| --- | ---: |
| 4 Fp multiplications (a zeta, a Re d_1, b Im d_1, B zeta) | 1,532 |
| 1 Fp addition (B) + 8 Fp subtractions/negations, 45 each | 405 |
| minimum of the 12 keys (12 top-limb loads, 11 x 3 compare/select, tie test, index) | 150 |
| load of the w^1 coefficient, kappa, twist index | 32 |
| **total counted 2,119, charged** | **2,119 = 0.10595 units** |

So a multiply step costs 1 + 0.0141 + 0.10595 = 1.12005 units and a squaring step
(9,803 + 282) / 20,000 = **0.50425** units.
Average step cost c = (5/8) * 1.12005 + (3/8) * 0.50425 = **0.889125 units**.

Setup and output (at most 2^44.08 units in all, the preprocessing bound): the **clearing pass over the
report store** of section 0-ii (2^52 buckets x 82 instructions, plus 4, below 2^44.07 units), and the
other items, at most 2^31 units together: pairing g and verification pairing (<= 2^14 each),
check z^r = 1, the table (1024 multipliers x (2 exponentiations <= 1,021 units + 11 Frobenius
maps) and the scalar tables, <= 2^20.1), M = 2^20 starts g^alpha z^beta (<= 2^30.1), the initial values of the
meter, the cap triple and the three one-word counters (at most 32 instructions, section 0-iv), final [x]P0 (<= 510 G1 operations,
each < 1 unit). DP handling: at most 2^55 reports (the cap of section 3.4) x **4,000 instructions**,
itemised as: 12 keys for kappa 2,119; two 255-bit modular multiplications by p^kappa mod r 1,600; the
one-word report counter increment, compare against the cap and branch 5; bucket index, bucket base and slot addresses
(64-bit) 10; **16 slot probes** (offset add, load of the used word, `tbz`, two 8-byte loads of F, `cmp`/`ccmp`,
branch) 128; storing the 80-byte record (`orr` of the used bit, ten 8-byte stores and address set-up) 13;
truncation branch after the 16th probe 4; reset of the
2-word stuck-walk counter 4; margin for loads, moves, call and return 117 — the worst case of the truncation rule of section 0-ii, with no
chain term. That is 2^55 * 4,000 / 20,000 < 2^52.679 units. The shared meter's updates (section 0-i)
add 2^20 * 40 / 20,000 < 2^11.04 units for the walk starts and 2^36 * 40 / 20,000 < 2^27.04 for the
restarts, the restart counter adds below 2^23.30 units (section 0-iv), and the restarted walks' fresh starts add
2^36 * 2^10.1 = 2^46.1 units. All non-walk work is therefore below
2^52.679 + 2^46.1 + 2^44.07 + 2^31 + 2^27.04 + 2^23.30 + 2^11.04 < 2^52.7 units (`experiment.py` part D
sums the exact terms, below 2^52.69675, and asserts the sum below 2^52.7).

Detection delay: after the first useful collision all M walks continue for 20/theta steps:
E = M * 20 / theta < 2^95.33 steps. The restart cap of section 0-iii adds one stuck-walk allowance per
restarted walk, 2^36 * 20/theta < 2^111.33, so the step bucket outside the walk is (2^20 + 2^36) * 20/theta < 2^111.33 steps
(2^-14.5 of S).

`W = c * S * (1 + eps) + c * (2^20 + 2^36) * 20/theta + 2^52.7 `, evaluated by `experiment.py` part D with the exact S = 2^125.830177..., 1 + eps = 1.000645319... and (2^20 + 2^36) * 20/theta = 2^111.321...: **W <= 2^125.6616283** (S is quoted elsewhere as 2^125.8302, its log2 rounded up to four decimals).

Submitted bound: **time_log2 = 125.682** (the least thousandth at or above the ledger recomputed with the
average step price raised by the 2-point drift factor of section 3.4, 2^125.6814751). It is at least 0.0203 bits (1.422%)
above the nominal ledger 2^125.6616283, which is the cap headroom of section 3.4. `experiment.py` (part D)
recomputes this ledger, and the two earlier price vectors of section 0 (2^125.7147 and 2^125.7075, at their own
theta = 2^-60 and the bucket input 2^63.8 of section 0) through the same code path, with exact decimal arithmetic, and checks the 125.682 inequality on integers.

Memory: the record store of section 0-ii is 2^56 slots x 80 bytes (128-bit fingerprint, two 256-bit
exponents, no pointer field) = 2^62.32 bytes, plus the 7.1 MB table, the 0.8 MB exponent tables and 2^20
walk states (at most 2 KB each, including the remembered previous state, 2^31 bytes): below 2^62.33 bytes,
and the claim states below 2^63 bytes. Memory is reported, not scored.

## 6. Heuristics and limitations

* H1 above (random-mapping behaviour of the mixed walk on classes), H2 (pseudorandom branch and DP bits
  of T), independence of table indices (used for cycle probabilities). These are standard in
  rho analyses (Teske 1998/2001; van Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999;
  Bos-Kleinjung-Lenstra 2010), supported here by toy and production-field runs, not proved.
* Scope of the prices. Three prices are measured with the cost model's own instrument and toolchain:
  blst_fp_mul (383), blst_fp12_mul (22,996) and blst_fp12_cyclotomic_sqr (9,803), each as a marginal
  over contiguous calls on real GT elements. The other charges are analytic instruction bounds built
  on the measured 383 and a stated 45-instruction modular addition: the 282 instructions of every-step
  work (270 of them itemised here plus the 12 of the meter and counter lines of section 0-i), the 2,119 of
  key work on a multiply step, and the 4,000 per report of section 5. The walk itself was not compiled and
  run as a blst program. The measured loop includes about 5 instructions per call of harness control that
  are not deducted (they are also inside the 282), so the squaring step is, if anything, over-priced. The
  measurement files shipped in `measure/` are the records of one session of the harness; the harness image
  could not be rebuilt for this revision, so they were not regenerated, and the equality of the stored
  limbs returned by the squaring and multiply branches (section 0a) is argued from the blst source, not
  measured. Sensitivity (part D): a squaring step priced at a full dense multiplication (1.0 unit) would
  give 2^125.936; each +0.1 unit on the squaring price adds 0.061 bits, each +100 instructions on the
  squaring routine 0.0031 bits and each +100 instructions of per-step work 0.0081 bits.
* Collision statistics are validated in toy fields only. The measured mean of the mixed walk is within
  5% of the model in all settings (pooled -1.4%); the plain-walk controls are -3.7% and +4.3%.
* The record store is one flat array of 2^62.32 bytes addressed by 64-bit arithmetic (section 0-ii). The track
  reports memory without scoring it, and the claim is an instruction-count bound: no physical AArch64 machine
  maps an array that large, and none is claimed.
* This is a generic square-root attack. It does not use index-calculus (TNFS) methods in Fp12.
  It is a valid upper bound, not a security estimate.

## 7. Experiments (`experiment.py`, standard library only)

* **Part A (BLS12-381):** GT membership; Frobenius formula = p-th power; 4-multiplication key
  formula = keys of the 12 explicitly computed conjugates, all distinct; kappa picks one identical
  canonical element from every conjugate; p^6 inverts; ord_r(p) = 12; 1+p^4+p^8 = 0; r-1
  factorisation; ord_r(2); the Granger-Scott formula equals y*y with 9 Fp2 squarings and 34 Fp2
  additions.
* **Part C (BLS12-381, production step rule):** a real walk of thousands of steps in GT with the
  production branch rule (branch and multiplier from the bits of Re(c_0), twisted table, squaring,
  2-cycle escape). Checks: every one of the first 150 steps lands in the class of y^ M_h or
  y^^2 as computed with explicit Frobenius; the exponent invariant y = g^alpha z^beta at the
  end; the squaring fraction (0.375 expected) and chi-square of the low bits of Re(c_0).
* **Part D (ledger and measurement record):** reads `measure/results.tsv` and the six raw Callgrind
  records, asserts the three per-call marginals (383, 22,996, 9,803) and that each raw record's total equals its
  table row, recomputes the ledger of section 5 for three price vectors through one code path
  (asserting the first two: 2^125.71472902377629 and 2^125.7074975), the cap of section 3.4 with its integer
  inequality, the drift table, the record store's load factor and truncation figure (section 0-ii), the
  word widths of every address, index and counter (asserted against 2^63 and 64 bits), the
  per-report charge and the meter charges of sections 0-i and 5, and the store clearing pass, the restart
  counter, the preprocessing bound and the non-walk sum of sections 0-ii, 0-iv and 5.
* **Part E (pair-collision law, toy fields):** sample m distinct classes, apply one step F to each and
  count coincidences of the images; N times the pair-collision rate is the quantity H1 predicts as
  1 - s^2 - (1-s)^2/(12R). Results (60,000 sampled classes, relative standard error about 0.9% except the last row):

  | q, R, s | observed | exact prediction | conservative value used | random mapping |
  | --- | ---: | ---: | ---: | ---: |
  | 67, 8, 3/8 | 0.8552 | 0.8553 | 0.8105 | 1 |
  | 67, 64, 3/8 | 0.8621 | 0.8589 | 0.8533 | 1 |
  | 67, 8, 0 | 0.9936 | 0.9896 | 0.8750 | 1 |
  | 67, 64, 0 | 1.0013 | 0.9987 | 0.9844 | 1 |
  | 67, 16, 1/2 | 0.7416 | 0.7487 | 0.7344 | 1 |
  | 211, 16, 3/8 (2.1%) | 0.8872 | 0.8573 | 0.8350 | 1 |

  The squaring branch removes exactly its s^2 of the collision probability, and the table branches
  keep 11/12 of theirs (rows with s = 0 and R = 8 separate the exact from the conservative value).
* **Part B (toy fields of the same tower shape, q = 67, 139, 163, 211, 379):** the complete
  algorithm (keys with tie-breaking by full conjugates, kappa trick, twisted table, squaring,
  2-cycle and checkpoint escapes, DPs, collision solving, verification g^x = z). Toy keys take
  few values, so the branch/DP bits are taken from a hash of the canonical conjugate (a
  stand-in for the bits of the 381-bit K; Part C tests those bits in the real field). Part B also reports the
  fluctuation of the mean charge per run against the independent-steps prediction (section 3.4).

Offline aggregate of the same script (all discrete logs recovered and verified; `trials` runs per
setting, `dp_bits` 3 and `walks` 4 unless noted; S(x) is the step budget at model success x):

| q, N classes, R, seed | squaring | runs | mean T_c / sqrt N (model) | P[T_c <= S(0.40)] | P[T_c <= S(0.42)] | P[T_c <= S(0.43)] |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| 67, 129,145, 64, seed 11 | 3/8 | 400 | 1.295 (1.357) | 0.438 | 0.468 | 0.475 |
| 211, 721,279, 256, seed 21 | 3/8 | 300 | 1.316 (1.353) | 0.433 | 0.443 | 0.453 |
| 163, 4,524,905, 64, seed 31 | 3/8 | 800 | 1.333 (1.357) | 0.415 | 0.436 | 0.445 |
| 163, 4,524,905, 64, seed 22 | 3/8 | 100 | 1.378 (1.357) | 0.390 | 0.400 | 0.410 |
| 163, 4,524,905, 64, seeds 61-63, walks 8, dp_bits 4 | 3/8 | 600 | 1.380 (1.357) | 0.405 | 0.422 | 0.433 |
| **pooled mixed walk** | 3/8 | 2,200 | 1.338 (1.357) | 0.418 | 0.437 | 0.447 |
| 67, 129,145, 64, seed 12 (control) | 0 | 400 | 1.217 (1.263) | 0.420 | 0.445 | 0.455 |
| 163, 4,524,905, 64, seed 32 (control) | 0 | 600 | 1.318 (1.263) | 0.370 | 0.385 | 0.395 |

The control rows are plain adding walks (their model values use the plain-walk effective class count). Measured
2-cycle rate at R = 64 with squaring: 4.5e-4 to 5.4e-4 per step across these settings (predicted
(5/8)^2/(12R) = 5.1e-4, bound 1/(12R) = 1.3e-3); without squaring 1.28e-3 and 1.35e-3 (1.30e-3 predicted).
Reproduce a row with e.g.
`echo '{"q":163,"trials":800,"R":64,"sq8":3,"dp_bits":3,"walks":4,"seed":31,"production_checks":0,"seconds":150}' | python3 -I experiment.py`.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999); Teske, Math. Comp.
70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Granger & Scott, PKC 2010; Aranha et
al., EUROCRYPT 2011; Bos, Kleinjung & Lenstra, ANTS 2010; mpjunior92, submission 7da5f23a.
