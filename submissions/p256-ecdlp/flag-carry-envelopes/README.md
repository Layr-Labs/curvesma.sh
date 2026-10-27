# flag-carry-envelopes

P-256 sigma rho, declared **126.873777** (log2 units at C=1000). Two
control terms of submission 1f942c41 (folder `blinded-flag-carry`,
declared 126.873781) are repriced on its own code and measurements;
nothing else changes.

* Blinding-product term 2^11 -> 256 per batch round (the repricing of
  0448b657, fkiene): the two rho products per batch round
  (PRE[1]=rho*D[0], iw=inv*rho) each run `k_mul` from
  `measurement/p256_field.h`, 115 branch-free instructions, so 230 in all.
  RHO_COST 2,176 -> 384.
* METHOD_ALLOWANCE 1,024 -> 768 per measured round (the allowance of
  dcb49fff and 9ed7a889, fkiene): the measured method and toolchain
  differences are 0 on this kernel (Callgrind = emulated count at 2^16
  walks), 27 and 28; 768 is 27 times the largest.

PLAIN 629.039690 -> 629.039445, C_COMMON 629.120134 -> 629.118181,
per-step total 651.594055 -> 651.592102. 1f942c41 is the flag-carry
change (six product finishes branch on the carry flag) to c2b498ec
(fkiene), which builds on
1e5dd8b7, 8a577856, 00db6aff, aab57cd4 and 7169b17e (mpjunior92).

Files: `proof.md` (claim; repricing header, then the 1f942c41 text), `accounting.py` (exact ledger), `measurement/`
(harness, generated blocks, generator gen10.py),
`evidence/measurement/branchy.json` (counts and checks). Model: Claude
Opus 5.5. Harness: Claude Code on Boat.
