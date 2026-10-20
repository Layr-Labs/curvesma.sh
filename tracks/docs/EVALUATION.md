# Exploratory review

Each track uses Eigen's shared `llm-verifier`, with one configured AI reviewer
and at most one follow-up experiment batch. Protected inputs are the track's
rules, target profile, verifier configuration and the shared
[Bernstein–Lange RAM specification](COST_MODEL.md).

Intake snapshots the submission and protected inputs. Declared experiments run
in the shared verifier's isolated, credential-free container. The reviewer reads
the complete claim, argument, source and captured evidence. It assesses the
algorithm, success probability and RAM bound together. Concrete, supported
heuristic attacks can qualify; missing details and unsupported heuristics cannot.

After acceptance, the verifier derives `ram_log2` from the first claim's exact
rational score. Every claim must meet the exploratory standard. Rejection or a
missing supported bound produces no score. Reviewer confidence is not algorithmic
success probability. Small experiments support only what they actually test.

Acceptance means plausible_not_refuted, not mathematical certification, novelty
or a practical production break. The review explains material uncertainty and
objections. No Lean proof, separate review committee or signed organizer receipt
is required. Human promotion through Yukon remains a separate decision.
