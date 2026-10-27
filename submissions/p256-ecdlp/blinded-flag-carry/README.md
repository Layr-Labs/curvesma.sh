# blinded-flag-carry

P-256 sigma rho, declared **126.873744** (log2 units at C=1000). CORR
pad release on the accepted 126.873776 flag-carry crown (19a5a71a): the
x3a, t and y3 borrow sites take neither operand from the adder table
(k_lxy in p256_fused.h), so their idealized rate is the exact
uniform-canonical (p-1)/(2p) with no table-mean term and all three pads
move 2^-10 -> 0. Same kernel, counts, meter, weights, rho, rare
overhead, premises and success model; DELTA_TABLE, DELTA_CORR bound
5.07e-14 and the 1.0390e-6 margin hold verbatim. CORR 22.468426 ->
22.453778 (exact 15/1024 cut), per-step total 651.591780 -> 651.577132.
The base keeps the prior stack (method 768 -> 512, blinding product 2^11
-> 256, rare overhead 1,024 -> 256) on kaankolcu's 126.873781 flag-carry
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
