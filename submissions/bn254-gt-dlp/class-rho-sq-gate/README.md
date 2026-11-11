# bn254-gt-dlp: gated distinguished-point rho on Frobenius classes with a squaring branch

Claimed time_log2 = 123.459421. This is the expected total work over uniform challenges and fresh
coins: one fresh coin opens the run or stops it at once. On every run the work is at most
2^133.287917 units. Success at least 0.39 under disclosed heuristics H1-H3. Peak memory at most
2^57.000001 bytes (one address space), no advice.

- [Algorithm, success analysis, ledger](proof.md)
- [Instruction counts and storage](machine.md): Fp12 product 13,843, cyclotomic squaring 5,238,
  class selector 1,067 and per-slot glue 366 executed instructions of generated AArch64 programs
- [Checks](experiment.py): kernels executed against big integers, algebra on the
  pinned constants, the full ledger, and the complete algorithm in a toy field
- [Toy walk](toyrho.py): the complete algorithm in BN-shaped toy fields
- [Offline toy runs](toy_offline.json): every instance of the runs in proof.md

A generic square-root bound with explicit accounting; not a sieve estimate and
not a practical computation.
