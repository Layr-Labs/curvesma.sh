# truncated-prefix BSGS v3

Retargets the accepted deterministic generic control to the exploratory
success threshold. Baby-table and giant-step extents are
`m = h = ceil(sqrt(ceil(39*r/100))) = 212506270031153001307260662621714878127`,
so coverage is exactly the prefix `[0,m^2)` with `100*m^2 >= 39*r`.
Kernels, record layout, sort, lookup and verification charges are unchanged
from the parent (`deterministic-bsgs-v1`).

- `proof.md` - construction, exact coverage, word-operation ledger.
- `claim.json` - profile binding and declared score `5132863/40000`.
- `submission.json` - claim statement and the experiment binding.
- `experiment.py` - production-arithmetic witness at reduced table size.

Declared: time_log2 128.321575, preprocessing_log2 127.321949,
memory_log2_bytes 135.320774, success probability 0.39.
