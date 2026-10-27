# additive-cheon-pairs-v1

Cheon's p−1 attack (d = 2^20) as a baby-step/giant-step search in the subgroup of d-th powers, with baby and
giant sets of small integers a·b in that subgroup (the construction of the accepted `additive-cheon-v1`,
submission cbadde77). New here: every character-class list is paired offline into three-term progressions
(x−g, x, x+g), and the walk computes P+T and P−T from one shared field inversion, so two points cost 9M
instead of 12M. Claimed log₂ W(2^20) = 116.580443, against 116.982765 for `additive-cheon-v1` and 116.9047
for the same accounting without the pair trick.

- `proof.md`: the algorithm, the planner, the non-degeneracy lemmas, the success argument, the instruction
  listings, the cost ledger, W(n) and the heuristics H1–H4.
- `evidence.py`: one standard-library script that recomputes the ledger and score, runs the planner on real
  class lists of a toy prime at the production degree parameter, runs the end-to-end toy attack with the
  shared-inversion pair formula, and runs the coverage census.

Prepared with AI assistance: Claude Code, model claude-fable-5-1, which was used for research, experiments
and drafting. The author reviewed the submission.
