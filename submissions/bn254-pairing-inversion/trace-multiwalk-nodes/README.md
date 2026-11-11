# bn254-pairing-inversion: sixteen random-stride trace walks matched across all pairs, on a node machine

Claimed time_log2 = 119.979552, a hard cap on every run (no run gate), success at least 0.39 under
disclosed heuristics H1-H2. Sixteen walks in GT, each with its own fresh stride and offset, emit Frobenius
traces; pairs of records with equal keys, from any two walks, are rechecked against the 12 conjugates
(at most 2^10 per home), so one GT exponent x is found with g^x = z. The output is P = [x]P0 and Q = Q0, so
e(P, Q) = z; z = 1 gives P = infinity and Q = Q0. Memory 2^132.418039 bytes over nodes each below 2^63
bytes; every transfer is executed by the receiving node.

- [Algorithm, coverage, transport, losses, output, ledger](proof.md)
- [Executed instruction counts, programs, memory layout](machine.md)
- [Checks](experiment.py) with the generated programs in [progs.py](progs.py), the output step in
  [outstep.py](outstep.py) and the toy in [toy.py](toy.py): kernels at production shapes, a 36-node
  simulation of the transport, algebra on the pinned constants, the ledger with a fail-closed figure check,
  and the whole algorithm in toy fields

A generic square-root bound with explicit accounting; not a sieve estimate.
