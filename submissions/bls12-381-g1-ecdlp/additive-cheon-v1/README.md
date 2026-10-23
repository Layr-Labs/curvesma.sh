# additive-cheon-v1

Cheon's p−1 attack (d = 2^20) as a baby-step/giant-step search in the subgroup of d-th powers. Both sets
are small integers a·b lying in that subgroup, walked in increasing order at one batched affine G1
addition per element. Claimed log₂ W(2^20) = 116.982765, against 117.894 for the KKM-accelerated Cheon attack under the
same per-element budgets.

- `proof.md`: the algorithm, the success argument, the cost ledger, W(n) and the heuristics.
- `evidence.py`: one standard-library script that recomputes the ledger and score, runs the toy coverage
  test and runs the end-to-end toy attack.

Prepared with AI assistance: Claude Code, model claude-opus-5-5, which was used for research, experiments
and drafting. The author reviewed the submission.
