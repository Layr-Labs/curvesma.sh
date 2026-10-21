# Fresh review and score history

The current shared verifier has no `reorg/plan.json`, signed-history import or
`archive_review.py` reassessment flow. Those belonged to the retired evaluator.
Do not provision those files or reuse its review artifacts as current acceptance.

After changing protected rules, target definitions or reference prices, add a
new candidate directory under `submissions/<track>/` and update its selector.
Existing candidate files are add-only. The new claim must justify its computation
under the current rules and be evaluated afresh. Preserve earlier candidates and
reviews as historical evidence, including unsuccessful attempts.

The two historical baseline packages use reference prices of 78,037,248
(secp256k1) and 45,839,872 (pairing) word operations. The current prices are
1,000 and 20,000 respectively; BLS G1 uses 2,000. The retired Miller track used 20,000.
Their old claimed scores therefore cannot serve as current baselines. A new
candidate must reconcile its full accounting with the current convention rather
than automatically relabeling or importing an old number. The any-pair and raw Miller tracks are retired; their results remain historical.

Current judging returns the submitted bound after acceptance; it does not silently
tighten the score. This does not guarantee that a fresh review preserves an old
score or decision. See the [historical qualification record](BASELINE_QUALIFICATION.md)
and [current cost model](COST_MODEL.md).
