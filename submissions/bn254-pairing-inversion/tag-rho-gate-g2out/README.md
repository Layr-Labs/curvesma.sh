# bn254-pairing-inversion: gated tag-walk rho with pairing output (P0, [x]Q0)

Claimed time_log2 = 119.136069. This is the expected total work over uniform challenges and fresh coins:
one fresh coin opens the run or stops it at once. On every run the work is at most 2^128.779227 units.
Success is at least 0.39 under the disclosed heuristics H1 and H2. Peak memory is at most 2^61.007988
bytes in one address space, and there is no advice.

The walk, store, ledger and evidence are those of the GT-DLP package 70974bf7 (source 6263134e). The
change is the output: after g^x = z, the algorithm outputs (P0, [x]Q0) with a counted G2 scalar
multiplication (proof.md section 0), charged inside the unchanged output row.

- [Algorithm, pairing output, tag lemma, success analysis, ledger](proof.md)
- [Instruction counts and storage](machine.md)
- [Checks](experiment.py) (declared initial run, the parent's), [interpreter and kernels](kern.py),
  [toy walk](toy.py), [offline toy runs](toy_offline.json)
- [G2 output check](g2out.py) (new, follow-up callable)
