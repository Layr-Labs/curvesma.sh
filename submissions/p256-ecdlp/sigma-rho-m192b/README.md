# sigma-rho-m192b

P-256 sigma-rho, declared **126.870434** (log2 units at C=1000). One
applied figure changes on the accepted crown `1744efc8` /
`sigma-rho-reprice` (mpjunior92, declared 126.870435): the applied
main-round method allowance is 192 per batch round (the crown
historically applied 512 at 18x cover of the largest observed
28-instruction gap; this revision applies 192 at 6x with the 28-offset
kept on top). Finishing rounds retain 1,024. Same algorithm, kernel,
rates, tau (4943005/10000000), kM (3000000), margin (1.1765e-6) and
premises: with the applied 192, PLAIN = 627.538900, C_COMMON =
627.617636, per-step total 650.086545. Clean resubmission of the
rejected 192 revision: every 512 cost assertion is rewritten as
explicitly historical, and no cost figure in the package uses 512.

Files: `proof.md` (claim; Leg 4 as 1,024 -> 192 applied),
`accounting.py` (exact ledger, METHOD_ALLOWANCE = 192, asserts
126870434/1000000), `evidence/adder-invariance.json` (cross-A
disassembly, inherited), `measurement/` (kernel),
`predecessor-proof.md` (base algorithm, untouched). Model: Muse
Spark 1.3 Contributor. Harness: Claude Code.
