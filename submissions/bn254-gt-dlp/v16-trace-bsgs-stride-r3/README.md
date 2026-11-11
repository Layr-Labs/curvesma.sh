# bn254-gt-dlp: trace baby-step giant-step with a random stride on a node machine

Claimed time_log2 = 120.534383, a hard cap on every run (no run gate), success at least 0.39 under
disclosed heuristics H1-H2. Memory 2^132.871418 bytes over nodes each below 2^63 bytes; every transfer is
executed by the receiving node.

- [Algorithm, coverage, transport, losses, ledger](proof.md)
- [Executed instruction counts, programs, memory layout](machine.md)
- [Checks](experiment.py) with the generated programs in [progs.py](progs.py): kernels at production shapes,
  a 22-node simulation of the transport, algebra on the pinned constants, the ledger with a fail-closed
  figure check, and the whole algorithm in toy fields

A generic square-root bound with explicit accounting; not a sieve estimate.
