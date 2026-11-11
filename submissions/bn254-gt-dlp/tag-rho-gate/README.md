# bn254-gt-dlp: gated tag-walk rho

Claimed time_log2 = 118.712713. This is the expected total work over uniform challenges and fresh
coins: one fresh coin opens the run or stops it at once. On every run the work is at most
2^133.320386 units. Success at least 0.39 under disclosed heuristics H1, H2 and H3. Peak memory at most
2^62.977583 bytes (one address space), no advice.

- [Algorithm, tag lemma, success analysis, ledger](proof.md)
- [Instruction counts and storage](machine.md): 119 executed instructions per walk step, one
  re-base with an Fp12 product every 3,700 steps
- [Checks](experiment.py): the priced programs executed against big integers, the tag lemma on
  production-size values, the full ledger, and the complete algorithm in a toy field
- [Interpreter and field kernels](kern.py), [toy walk and tag arithmetic](toy.py)
- [Offline toy runs](toy_offline.json): every instance of the runs in proof.md

A generic square-root bound with explicit accounting; not a sieve estimate and
not a practical computation.
