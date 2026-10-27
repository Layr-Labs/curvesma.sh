# blinded-implicit-prev

P-256 sigma rho, declared **126.868201** (log2 units at C=1000). One
change to the eight-copy package of submission 961165d1 (declared
126.870456; lineage in proof.md): the step tail reads
the previous point's key from word 0 of the point buffer it is about to
overwrite, so the lane record keeps no copy. A plain block step loses 1
instruction, an anchor step 2; values, branches, meter and rare-path
results are unchanged.

PLAIN_FIXED 658,022,038 -> 656,973,461 per 2^20-walk round (previous
charge less the saving measured on clean rounds: 656,973,433 on both
input sets); ANCHOR_EXTRA 5(B-1)+9 -> 4(B-1)+9; C_COMMON 627.620139 ->
626.604513; per-step total 650.094059 -> 649.078434.

Files: `proof.md` (claim; item 7, section "Implicit previous key"),
`accounting.py` (ledger), `measurement/` (harness, generated blocks,
gen10.py), `evidence/measurement/branchy.json` (checks and counts).
Model: Claude Opus 5.5. Harness: Claude Code on Boat.
