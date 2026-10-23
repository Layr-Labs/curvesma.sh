# nfs-stnfs-v1: STNFS discrete-log inversion in GT

Reduces e(P,Q)=z to a discrete logarithm in the prime-order group GT
(base g=e(P,Q0)), solves it with the published Special-Tower Number Field
Sieve in F_{p^12} using the Guillevic-Singh BLS12-381 polynomial and Table 10
parameters, returns Q=[a]Q0 and verifies the pairing. Success probability 1
(Las Vegas, retries charged); total computation 2^117.858891
target-group-operation equivalents after a conservative instruction and
overhead conversion of the published 2^126 cost; peak memory 2^72.389801
bytes; preprocessing 2^115.858891 (included, no reuse credit); advice 0.

See `proof.md` for the reduction, the published-parameter reproduction, the
instruction ledger and its margin over the published spread. `experiment.py`
recomputes the GS21 cost and the submitted ledger from the Table 10
parameters; it is arithmetic evidence, not an attack run.
