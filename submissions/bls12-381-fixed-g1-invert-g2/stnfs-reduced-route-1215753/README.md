# STNFS with 99.9% ECM, a measured record route and fully reduced LA kernels: 121.5753

Track bls12-381-fixed-g1-invert-g2: given nonidentity P in G1 and z in GT, return Q in G2 with e(P,Q)=z.

04073dd1 (121.8936) plus: the record route measured end to end (10.1 instead of 36 per record; relation collection 2^121.2533), and linear algebra with every mod-r reduction priced by compiled kernels (`row_tail.c`), the row tail charged per row, six-limb sums valid for any row, and pad 1.10 on the measured N/U scale trend (2^117.7327). Ledger <= 2^121.57530, success >= 0.391918 (proof.md, STORAGE.md). Analytic exploratory claim, not a full-size run.

Entry points: audit_refinement.py (bounded audit), unit_entries.py (H4u, U/N), evidence16.py (full diagnostic, 21 checks).
