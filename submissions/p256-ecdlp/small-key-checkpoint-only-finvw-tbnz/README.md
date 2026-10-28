# small-key-checkpoint-only-finvw-tbnz

P-256 sigma rho, declared **126.852665** (C=1000). Base: ae463634
(126.852724, checkpoint-only cycle detection). Two fixed-path edits:
`finv_w`, a 4-bit window over p-2 (-11,011 per 2^16-walk plain round), and
one `tbnz` for the 2^60 cap test (-(2^16-1) per anchor round).
PLAIN_FIXED 654,876,307 -> 654,865,296;
ANCHOR_EXTRA 4(B-1)+9 -> 3(B-1)+9; C_COMMON 624.604511 -> 624.578384;
per-step 642.149317 -> 642.123191. The base build reproduces 40,967,797
and 41,229,944; patched rounds are 40,956,786 and 41,153,398 with equal
meters and identical dumps. See `proof.md`, `accounting.py`,
`evidence/measurement/finvw-tbnz.json`.
