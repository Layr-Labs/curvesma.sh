# Ordinary BLS12-381 G1 DLOG via randomized-base GT transfer and a revised STNFS core: 124.9650

Read proof.md (G1 claim), REVISION.md (changed core numbers), source-core.md with DESCENT.md, KERNELS.md, STORAGE.md (the GT STNFS core) and field-map.md (charged Fp12-to-tower map).

Initial experiment: `audit_g1.py` (empty input, about 2 s); follow-ups {"descent": true} (about 10 s) or {"unit_entries": [15, 1]}. Local full diagnostic: `evidence16.py` {"samples": 16, "seed": 1} (about 20 s).

Analytic exploratory claim under H1-H5; not a full-size run.
