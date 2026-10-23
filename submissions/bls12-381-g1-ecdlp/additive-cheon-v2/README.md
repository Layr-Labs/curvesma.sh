# additive-cheon-v2

Cheon's p−1 attack (d = 2^20) as a baby-step/giant-step search in the subgroup of d-th powers. Both sets
are small integers a·b lying in that subgroup, walked in increasing order at one batched affine G1
addition per element. Claimed log₂ W(2^20) = 116.889745, against 117.862 for the KKM-accelerated Cheon
attack under the same per-element budgets.

This is a revision of the published package of a submission reviewed valid at log₂ W = 116.982765
(2026-10-23T14:42Z). The algorithm, the coverage argument, the heuristics and the toy experiments are
that package's, unchanged. The revision is in the cost accounting:

- the key store is one flat block array in a single extended-memory heap, addressed by a two-level
  radix page table whose entries are two words (16 bytes) at both levels, so the level-2 entry holds a
  123-bit block base and the level-1 entry a 111-bit table base; the two-word index build, the second
  word of each entry and the record's two-word address and advance are itemised in `proof.md` §4.3;
  this replaces a table sharded over nodes with store-and-forward;
- every access is charged at its **worst case**: eight probes in a fixed-depth window plus four
  comparisons against the block's bounded overflow area. This replaces the average probe count of
  linear probing, which did not bound a data-dependent probe sequence. A key that finds neither is
  discarded by a stated rule, and the discarded fraction is charged at 0.01 and paid for by raising the
  coverage exponent from 0.52 to 52/99;
- every per-element budget is the itemised instruction count instead of a padded envelope over it:
  26 + 29 + 106 + 80 = 241 instructions for a baby element and 26 + 29 + 253 + 80 = 388 for a giant one.

The six changes are listed in `proof.md` §4.6.

- `proof.md`: the algorithm, the success argument, the store layout and addressing (entry widths,
  bit budgets and every multiword instruction), the probe and discard rule, the cost ledger, W(n) and the
  heuristics.
- `evidence.py`: one standard-library script that recomputes the ledger, the per-element budgets, the
  exact score and an exact integer proof of its rounding, measures the store's window and discard
  fractions at the production load, runs the toy coverage test and runs the end-to-end toy attack.

Prepared with AI assistance. The source package was drafted with Claude Code (model claude-opus-5-5);
this revision's cost accounting and text were drafted with an AI coding agent (DeepSeek V4.1 Flash on
opencode). The author reviewed the submission.
