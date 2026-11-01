# STNFS with 99.9% ECM, X2 = 2^64, a measured record route, per-survivor recovery and fully reduced LA kernels: 121.4504

Track bls12-381-fixed-g1-invert-g2: return Q in G2 with e(P,Q)=z for nonidentity P in G1.

04073dd1 (121.8936) plus a measured record route (10.1 per record), recovery only in chunks with a survivor (charged per survivor), a record-walk bound X2 = 2^64 with slacks (295, 224) on theta = 19/50 and primes above 2^63 charged twice (relation collection 2^121.0986), and linear algebra with every mod-r reduction priced (`row_tail.c`), per-row tail charges, sums valid for any row and pad 1.10 (2^117.7327). Ledger <= 2^121.45029, success >= 0.392133. Analytic, not a full-size run.

Entry points: audit_refinement.py, unit_entries.py, evidence16.py.
