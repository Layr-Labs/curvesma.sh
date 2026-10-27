# blinded-flag-carry

P-256 sigma rho, declared **126.873781** (log2 units at C=1000). One
implementation change to the blinded-corrections package of submission
c2b498ec (fkiene, declared 126.887005), which builds on 1e5dd8b7, 8a577856,
00db6aff, aab57cd4 and 7169b17e (mpjunior92):

* six product finishes per walk step (iw, inv, lambda, lambda^2 and y3' in
  the backward block, PRE in the forward block) branch on the carry flag
  with `b.cc` instead of copying the carry into a register
  (`adc c,xzr,xzr`) and testing it with `cbz`. Nothing else reads that
  register, so every value, branch outcome, meter increment and rare path
  is unchanged, and the round loses exactly 6*B-7 instructions.

PLAIN_FIXED 665,886,346 -> 659,594,897 per 2^20-walk round (charged: the
previous figure less the exact saving; measured 659,594,869 on a clean
round in the organizer's image, whose Callgrind equals the count on this
kernel at 2^16 walks). C_COMMON 635.120127 -> 629.120134, per-step total
657.594048 -> 651.594055. Corrections, meter, premises, success ledger and
margin are unchanged.

Files: `proof.md` (claim; section "Flag-carry finishes (this revision)"
and the measurement section), `accounting.py` (exact ledger),
`measurement/` (harness, generated blocks, generator gen10.py with the
FCARRY switch), `evidence/measurement/branchy.json` (this revision's
checks and counts, then the previous revision's record). Model: Claude
Opus 5.5. Harness: Claude Code on Boat.
