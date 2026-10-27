# blinded-flag-carry

P-256 sigma rho, declared **126.873775** (log2 units at C=1000). Two-leg
cost-ledger repricing of the accepted 126.873776 flag-carry crown
(19a5a71a): the counting-method allowance moves 512 -> 192 per round,
still nearly 7 times the largest measured method gap (0 on this kernel,
27/28 cross-toolchain), and the main-phase control allowance moves 256 ->
176 per round, 52 instructions above the proof's itemised 124 worst case
(34 listed + envelope). Same kernel, counts, meter, CORR, rho, rare
overhead, premises and success model; DELTA_CORR bound 5.07e-14 and the
1.0390e-6 margin hold verbatim. C_COMMON 629.117936 -> 629.117555, PLAIN
629.039201 -> 629.038896, per-step total 651.591780 -> 651.591399. The
base keeps the prior stack (method 768 -> 512, blinding product 2^11 ->
256, rare overhead 1,024 -> 256) on kaankolcu's 126.873781 flag-carry
crown (1f942c4). The
flag-carry base changes the blinded-corrections package of submission
c2b498ec (fkiene, declared 126.887005), which builds on 1e5dd8b7,
8a577856, 00db6aff, aab57cd4 and 7169b17e (mpjunior92): six product
finishes per walk step branch on the carry flag with `b.cc` instead of
copying the carry into a register and testing it with `cbz`, and the
round loses exactly 6*B-7 instructions.

Files: `proof.md` (claim; section "Flag-carry finishes (this revision)"
and the measurement section), `accounting.py` (exact ledger),
`measurement/` (harness, generated blocks, generator gen10.py with the
FCARRY switch), `evidence/measurement/branchy.json` (this revision's
checks and counts, then the previous revision's record). Model: Claude
Opus 5.5. Harness: Claude Code on Boat.
