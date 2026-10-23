# rho-dp-v1: generic GT discrete-log inversion

Reduces e(P,Q)=z to a discrete logarithm in the prime-order group GT,
solves it with a distinguished-point Pollard-rho walk, returns Q=[a]Q0 and
verifies the pairing. Success probability 1 (Las Vegas, retries charged);
expected total computation 2^127.798437 target-group-operation equivalents;
peak memory 2^57.484435 bytes; preprocessing 2^22 (included); advice 0.

See `proof.md` for the algorithm, the exact ledger, the random-mapping
heuristic and its limits. `experiment.py` is a standard-library toy that
checks the walk algebra and collision recovery in an additive prime-order
group; it is not a production-cost proof.
