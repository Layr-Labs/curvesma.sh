# P-256: asymmetric forward-loop unrolling

This submission retains the complete P-256 algorithm in source-proof.md
and its supplied predecessor-proof.md, probability.md, accounting.py and
implementation. The immutable source is submission
e0d1ed74-946e-4f8d-a5ad-a6d219ad3697, commit
07caf4005e166f2db39abfccd1755b30d08eea63. No previous verdict is a premise.

Only the forward prefix loop changes. It now has 1024 copies of the identical
FLOOP_S body, entered at copy 2. The backward loop, jump table and all its
rare exits remain at 512 copies. The forward counter starts at B; B is a
multiple of 1024. Its first block executes 1022 bodies; each subsequent block
executes 1024. There are B/1024 blocks. Thus the total body count is B-2,
and the lane order remains 1 through B-2 exactly as before. The counter x26
does not appear in the body. All pointer updates, field arithmetic, correction
meters, stores and register state at each processed lane are unchanged.

The only recurring control instructions are SUBS counter,#1024 and B.NE.
They cost two per block, compared with two per 512-copy block previously.
The deterministic saving is 2(B/512-B/1024)=B/512 instructions per forward
call, on all three outer schedule arms. The branch into copy 2 is unchanged.
There is no probability assumption in this saving. At B=2^28 it is 524288
instructions per round, or exactly 1/512 per lane.

The forward source body has at most 219 instructions (counting every escaped
newline, including labels, conservatively). Therefore its 1024 copies and
control consume fewer than 2^20-4 bytes. The outer conditional backward
branch is within its signed 19-bit scaled range. Body-local conditional
branches remain within a single copy. The entry unconditional branch is
within its much larger range. No slow stubs or jump table are moved here;
these belong to the unchanged backward block. The immediate 1024 is legal
for SUBS. This bound deliberately does not double the much larger backward
block, whose conditional-branch ranges would require a separate analysis.

PLAIN_FIXED decreases by B/512. All other cost, probability and memory
formulas are retained, including all failure reserves and the complete
preprocessing charge. The source's older numerical common-path and total
work figures are replaced only by the recomputed ledger below.
No field multiplication, correction allowance, tail failure or future target
amortization is removed. Code storage increases by a bounded constant; the
source's separate 2^24-byte allowance is more than this extra forward block.

The exact schedule check covers small and production batch sizes, explicitly
compares small lane sequences and conservatively checks source instruction
range. It is not a fresh compiled timing or instruction measurement. The
source measured body pricing is inherited with the new control count as an
analytical delta, just as the source's 128-to-512 extension was analytical.
All disclosed H1, M1/L-ori, H3 and H-trail production limitations remain.

## Final ledger and checks

The complete isolated ledger returns log2 W = 126.72859355679839; we declare
the least nine-decimal value at or above it, 126.728593557. The conditional
success lower bound is 0.3900000219 (tau 0.4943071, certified), above 0.3900000001. Preprocessing is at most 2^89.961190 units and memory
below 2^77 bytes. The source proof's old headline is superseded by this
claim; its scientific model and assumptions are not weakened.

All seven schedule/control checks passed, including B=2^28. The conservative
forward-block size is 897056 bytes, below the conditional-branch range. The
inherited 80-trial toy experiment recovered 30 scalars, all verified, with
20885 replay checks and no mismatches. Its 0.375 success frequency and wide
interval are toy evidence, not a production success bound. All original
probability premises remain necessary. Accounting comments were shortened
for packaging only; no scientific assertion was removed. The old equality
assertion for the exact former score is replaced by a strict improvement
check; all other arithmetic and probability assertions are retained.
