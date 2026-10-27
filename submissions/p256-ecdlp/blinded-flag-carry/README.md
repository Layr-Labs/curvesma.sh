# blinded-flag-carry

P-256 sigma rho, declared **126.879314** (log2 units at C=1000). Analysis
repricing of kaankolcu's 4ff4847 flag-carry crown (126.879316): the
birthday-approximation term delta moves 2^-20 -> 2^-24, still 2^56 over
the true T^3/(6N^2) ~ 7e-40 approximation error, and the freed 0.89e-6
of success margin lets tau step down its 1e-7 grid from 0.4943041 to
0.4943026 (4943025 fails the 1.0e-6 floor). Same kernel, counts, meter,
per-step prices, CORR, premises and success-model form; per-step charges
bit-identical, margin 1.0390e-6 -> 1.0181e-6 above 0.39, preprocessing
62.848071 -> 62.848069. The flag-carry base is mpjunior92's 00db6aff
(blinded-corrections, 126.89249; extends aab57cd4) with the six product
finishes branching on the carry flag (`b.cc` in place of
`adc c,xzr,xzr; cbz c`), 5*(2^20-1)+(2^20-2) fewer instructions per
2^20-walk round, measured exactly. Fixed path 631.04 + metered
corrections 22.47 per walk step (base 637.04). Everything else, proof
and measurements included, is that package's and its lineage's.
Files: `proof.md`, `accounting.py`, `measurement/` (gen8.py, headers),
`evidence/measurement/branchy.json` (`flag_carry`).
