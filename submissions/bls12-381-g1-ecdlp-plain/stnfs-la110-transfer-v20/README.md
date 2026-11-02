# Ordinary BLS12-381 G1 DLOG via revised STNFS core (LA 1.10, exact prices, line sieve 2^45.5 rerun, slack 60, UDUP 0.45): 124.617970069062

Read proof.md (G1 algorithm, distribution, success, work), REVISION.md (every changed core number), source-core.md with DESCENT.md, KERNELS.md, STORAGE.md (the GT STNFS core) and field-map.md (charged Fp12-to-tower map).

Initial experiment: `audit_g1.py` (empty input, about 2 s); follow-ups {"descent": true} (about 10 s) or {"unit_entries": [15, 1]}. Local full diagnostic (not declared: about 19 s): `evidence16.py` with {"samples": 16, "seed": 1}.

Analytic exploratory claim under H1-H5; not a full-size run.
