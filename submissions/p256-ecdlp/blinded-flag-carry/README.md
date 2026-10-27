# blinded-flag-carry

P-256 sigma rho, declared **126.873777** (log2 units at C=1000). Pure
cost-ledger repricing of kaankolcu's 126.873781 flag-carry crown
(1f942c4): three allowances move to their listed or measured support
with no shared fate. Blinding-product term 2^11 -> 256 per batch round
(two listed 115-instruction branch-free k_mul rows plus the rho reads).
Rare-event and look-ahead overhead 1,024 -> 256 per step (forced-path
own control at most 191 per 2-cycle check, 107 + 59 per anchor member).
Counting-method allowance 1,024 -> 768 per round (28x the measured
27-per-round gap). Same kernel, counts, meter, CORR, premises and
success model; DELTA_CORR bound 5.07e-14 and the 1.0390e-6 margin hold
verbatim. C_COMMON 629.120134 -> 629.118180, per-step total 651.594055
-> 651.592024 (rare reserve 0.005494 -> 0.005418 inside it). The
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
