# p256-ecdlp: Pollard rho with distinguished points and the negation map

Generic parallel rho on P-256 with class size 2, batched affine steps and fruitless-cycle escapes
made a function of the cycle (its smallest-u member).
Total work 2^127.82 units, with each P-256 step charged 1193 AArch64 instructions (1.193 units).
Success at least 0.39 for any walk quality factor at or above 0.908 (the ideal-mapping collision
probability is 0.42); memory at most 2^41 bytes; no advice. See proof.md for the full ledger and
p256_mul_count.md for the instruction listing that fixes the field-operation counts
(M = 146 per Montgomery multiplication, A = 21 per modular add/sub, each counted at 139 and 20
instructions plus a 5 % reserve; check_mul_mont.py executes the listing under the AArch64 flag
rules, where only adds, adcs, subs and sbcs write the carry flag). experiment.py runs the same walk
and the same cycle rule on a toy curve in under a second; proof.md section 7 lists what it shows and
why that is toy-scale support only.

The cycle rule in one paragraph. Each walker keeps the full 256-bit canonical u of a checkpoint,
refreshed every 64 steps, and compares every new point with all 256 bits: a match is a return to
that exact class, hence a cycle, and a false match would need a random point to reproduce a given
256-bit u (2^-256 per step, 2^160 below the distinguished-point rate). The walker then walks the
cycle once, keeps the member with the smallest u with its full state, and adds one escape point
from a separate 256-entry table with independent coefficients, so the escape is a function of the
cycle and merged walks leave together. Re-entry into the cycle just left is bounded by the basin of
the cycle: at most 64 x 2^113.5648 classes within the walker's reach over 2^255, so 2^-135.4 per
cycle and 2^-32.9 over all cycles of a run. Cycles longer than 64 are never detected; the whole
cycle loss is bounded at 2^-9.19 of all steps (detected cycles 2^-16.87, re-entry 2^-32.87,
undetected long cycles 2^-9.22) and the step budget uses 2^-9.

Revision 5 (2026-10-23): closes the hostile-review gaps. The checkpoint now stores and compares
the full 256-bit canonical u instead of 64 bits (a 64-bit checkpoint matched at 2^-64 per step,
2^32 times the distinguished-point rate, so merged trails split); the compare is charged in the
ledger (cycle line 20 instructions against 15.05 itemised, per-step total 1189 -> 1193). Re-entry
into a cycle just escaped is now bounded by the basin argument, and the escape table is stated to
be independent of the main table. Long cycles are bounded over every length up to the walker's
step cap instead of stopping at 2^20, which gives 2^-9.22 of all steps for that term; the step
budget's cycle allowance rises from 2^-12 to 2^-9 and K falls from 2^127.5623 to 2^127.5648. The
walk-quality margin is now computed rather than asserted: 0.39 holds for every quality factor at
or above 0.9076, and the finite-r Brent-Pollard correction and Teske's measured deviations are
cited. A toy experiment is added. The identity target returns 0. The field-operation counts are
unchanged (139 and 20, charged 146 and 21); the score stays 127.82 (127.8194 before rounding, up
from 127.8121, still below the 128.2 ceiling).

Revision 4 (2026-10-23): carry-chain instructions in the listing are now the flag-setting adcs/sbcs
(the count stays 139) and the checker models the real flag rules; walkers start from 64-bit
coefficients so the (A, B) bound holds; the table cost now includes field additions
(preprocessing 2^32.9); the non-multiplication per-step lines are itemised as upper bounds. The
score is unchanged at 127.82.

Revision 3 (2026-10-23): the field-operation costs are now counted from an explicit AArch64
instruction listing rather than estimated. The listing is 139 instructions for a 4 x 64-limb
Montgomery multiplication (the earlier text estimated about 185 and charged 220) and 20 for a
modular addition (charged 21 with the same reserve, previously 24). The per-step charge falls
from 1664.03 to 1189 instructions and the claim from 128.30 to 127.82. Nothing in the attack,
its success argument or its memory changes.

Revision 2 (2026-10-22): resubmitted unchanged in substance after the challenge's evaluation
workflow could not be dispatched during a GitHub Actions outage; the claim, proof and score are
identical to revision 1.
