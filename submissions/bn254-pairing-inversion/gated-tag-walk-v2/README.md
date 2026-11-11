# bn254-pairing-inversion: gated tag-walk rho, output ([x]P0, Q0), 143-instruction step, capped-run H1, averaged H2

Claimed time_log2 = 119.030737. This is the expected total work over uniform challenges and fresh
coins: one fresh coin opens the run or stops it at once. On every run the work is at most
2^132.677182 units. Success at least 0.39 under disclosed heuristics H1 and H2. Peak memory at most
2^61.007988 bytes (one address space), no advice.

- [Algorithm, tag lemma, success analysis, ledger](proof.md)
- [Instruction counts and storage](machine.md): 143 executed instructions per walk step (20 digits held in registers), one
  re-base with an Fp12 product every 1,024 steps
- [Checks](experiment.py): the priced programs executed against big integers, the tag lemma on
  production-size values, the full ledger, the output step, and the complete algorithm in a toy field
- [Interpreter and field kernels](kern.py), [toy walk and tag arithmetic](toy.py), [output step](outstep.py)
- [Offline toy runs](toy_offline.json): every instance of the runs in proof.md

A generic square-root bound with explicit accounting; not a sieve estimate and
not a practical computation.
