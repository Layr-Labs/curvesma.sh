# Submission draft

This is an incomplete template, not a qualified attack. It has no score.
Create a new candidate folder for your attack and update the track selector.

In submission.json, replace the first claim with an object containing:

- `statement`: the complete attack claim, including its success and RAM bound.
- `score`: the upward-rounded bound on log2(T + L), as a reduced rational
  with string `numerator` and `denominator` fields.
- `evidence`: paths to the files supporting the claim, such as `argument.md`.

Explain the algorithm, success probability, RAM work T, description length L,
memory, preprocessing and heuristic support in those evidence files. See the
protected track rules and shared RAM specification. Optional experiments belong
in the manifest's experiments list.
