# bn254-pairing-inversion: gated tag-walk rho with a two-level signed tag, output ([x]P0, Q0)

Claimed time_log2 = 118.828325. This is the expected total work over uniform challenges and fresh
coins: one fresh coin opens the run or stops it at once. On every run the work is at most
2^131.353421 units. Success at least 0.39 under disclosed heuristics H1 and H2. Peak memory at most
2^61.353586 bytes (one address space), no advice.

A modification of kaankolcu's packages 70974bf7 and 271f78c0; credits for the borrowed parts are in proof.md.

- [Algorithm, two-level tag lemma, success analysis, ledger](proof.md)
- [Instruction counts and storage](machine.md): 120 executed instructions per walk step,
  177 more when the second tag level is needed, one re-base with an Fp12 product every 1,024 steps
- [Checks](experiment.py): the priced programs executed against big integers, both tag lemmas on
  production-size values, the full ledger, the output step, and the complete algorithm in a toy field
- [Interpreter and field kernels](kern.py), [toy walk and tag arithmetic](toy.py), [output step](outstep.py)
- [Offline toy runs and census](toy_offline.json): every instance of the runs in proof.md

A generic square-root bound with explicit accounting; not a sieve estimate and
not a practical computation.
