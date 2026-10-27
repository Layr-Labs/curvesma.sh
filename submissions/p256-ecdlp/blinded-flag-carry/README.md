# blinded-flag-carry

P-256 sigma rho, declared **126.873776** (log2 units at C=1000). One-line
cost-ledger repricing of the accepted 126.873777 flag-carry crown
(0ba7d676): the counting-method allowance moves 768 -> 512 per round,
still 18 times the crown's own measured 27-per-round QEMU-vs-Callgrind
gap (organizer image Callgrind equals the emulated count on this kernel
at 2^16 walks). Same kernel, counts, meter, CORR, rho, rare overhead,
premises and success model; DELTA_CORR bound 5.07e-14 and the 1.0390e-6
margin hold verbatim. C_COMMON 629.118180 -> 629.117936, PLAIN
629.039445 -> 629.039201, per-step total 651.592024 -> 651.591780. The
base keeps the prior stack (blinding product 2^11 -> 256, rare overhead
1,024 -> 256, method 1,024 -> 768) on kaankolcu's 126.873781 flag-carry
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
