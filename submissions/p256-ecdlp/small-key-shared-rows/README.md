# small-key-shared-rows

P-256 sigma rho, declared **126.741857** (C=1000), on the y3-fold
revision (126.789581): two products of one factor share its Montgomery
rows (78 instead of 90 instructions each); the forward block also forms
Q = dy*PRE and the backward block takes lambda = INV*Q (iw gone),
PLAIN_FIXED 629,710,511 -> 608,739,081, step 594.650484. See proof.md
"Shared rows", `accounting.py`, `measurement/gen49.py` and
`evidence/measurement/branchy.json` (`shared_rows`).
