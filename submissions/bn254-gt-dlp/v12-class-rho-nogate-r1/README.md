# bn254-gt-dlp: distinguished-point rho on Frobenius classes

Claimed time_log2 = 124.763426, a hard cap on every run (no run gate), success at least 0.39 under
disclosed heuristics H1-H3. Peak memory below 2^47.000001 bytes (one 48-bit address space), no advice.

- [Algorithm, success analysis, ledger](proof.md)
- [Instruction counts and storage](machine.md): Fp12 product 13,843, class
  selector 1,067 and per-slot glue 357 executed instructions of generated AArch64 programs
- [Checks](experiment.py): kernels executed against big integers, algebra on the
  pinned constants, the full ledger, and the complete algorithm in a toy field
- [Offline toy runs](toy_offline.json): every instance of the runs in proof.md

A generic square-root bound with explicit accounting; not a sieve estimate and
not a practical computation.
