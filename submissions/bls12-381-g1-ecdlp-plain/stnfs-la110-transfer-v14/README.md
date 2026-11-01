# Ordinary BLS12-381 G1 DLOG via randomized-base GT transfer and a revised STNFS core (linear algebra with pad 1.10): 124.85609

Read proof.md (G1 algorithm, distribution, success, work), REVISION.md (every changed core number), source-core.md with DESCENT.md, KERNELS.md, STORAGE.md (the GT STNFS core) and field-map.md (charged Fp12-to-tower map).

Initial experiment: `audit_g1.py` (empty input, about 2 s); follow-ups {"descent": true} (about 10 s) or {"unit_entries": [15, 1]}. Local full diagnostic (not declared: about 19 s): `evidence16.py` with {"samples": 16, "seed": 1}.

Analytic exploratory claim under H1-H5; not a full-size run.
