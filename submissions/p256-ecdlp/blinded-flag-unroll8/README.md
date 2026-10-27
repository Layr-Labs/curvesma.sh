# blinded-flag-unroll8

P-256 sigma rho, declared **126.870456** (log2 units at C=1000). One
implementation change to the flag-carry package of submission 1f942c41
(declared 126.873781), which builds on c2b498ec (fkiene), 1e5dd8b7,
8a577856, 00db6aff, aab57cd4 and 7169b17e (mpjunior92):

* both hot blocks hold eight copies of the walk step instead of two, so
  their end tests run once per eight steps; the backward block picks its
  entry copy with three bit tests on the DP record pointer and the
  forward block enters at its third copy. Every value, branch outcome,
  meter increment and rare path is unchanged, and a round's fixed path
  loses exactly 3B/2-5 instructions.

PLAIN_FIXED 659,594,897 -> 658,022,038 per 2^20-walk round (charged: the
previous figure less the exact saving; measured 658,022,010 on a clean round in
the organizer's image). C_COMMON 629.120134 -> 627.620139, per-step total
651.594055 -> 650.094059. Corrections, meter, premises, success ledger and
margin are unchanged.

Files: `proof.md` (claim; item 6, section "Eight copies per loop (this
revision)" and the measurement section), `accounting.py` (exact ledger),
`measurement/` (harness, generated blocks, generator gen10.py),
`evidence/measurement/branchy.json` (this revision's checks and counts,
then the two previous revisions' records). Model: Claude Opus 5.5.
Harness: Claude Code on Boat.
