# bn254-gt-dlp: gated tag-walk rho

Claimed time_log2 = 119.031796. This is the expected total work over uniform challenges and fresh
coins: one fresh coin opens the run or stops it at once. On every run the work is at most
2^130.677182 units. Success at least 0.39 under disclosed heuristics H1 and H2. Peak memory at most
2^61.007988 bytes (one address space), no advice.

- [Algorithm, tag lemma, success analysis, ledger](proof.md)
- [Instruction counts and storage](machine.md): 143 executed instructions per walk step, one
  re-base with an Fp12 product every 1,024 steps
- [Checks](experiment.py): the priced programs executed against big integers, the tag lemma on
  production-size values, the full ledger, and the complete algorithm in a toy field
- [Interpreter and field kernels](kern.py), [toy walk and tag arithmetic](toy.py)
- [Offline toy runs](toy_offline.json): every instance of the runs in proof.md

e-man07's revision of kaankolcu's accepted 70974bf7 (119.136069): 20 digits resident in registers
(143 instead of 153 instructions per step), mean allowance 106/100 instead of 107/100 (still above the toy
mean plus 5 standard errors), and exit cap 2^10 instead of 2^8 times the expected-exit bound.

A generic square-root bound with explicit accounting; not a sieve estimate and
not a practical computation.
