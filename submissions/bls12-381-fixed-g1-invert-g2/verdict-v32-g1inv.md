# g1-invert v32 verdict -- ticket bca1e72a (REJECTED)

Submission 8b82dd99 and 3817eb37/376a6935 died at INTAKE (protected-config step).
bca1e72a (fire 4, after the title fix) reached the judge: REJECTED, 3256 chars.

Banked intake class: 'protected-config step failure with NO archive-size message =
run gh run view --log-failed -j <jobid> -R Layr-Labs/curvesma.sh | grep llm-verifier
-- bounded_text_required on title is a real intake-killer class.' Fire log: title
312 chars (non-bounded) killed 8b82dd99; 249 chars passed intake. Score schema also
banked: claims[0].score = {"numerator": "<str>", "denominator": "<str>"} (an int
pair gives canonical_rational_required; a bare "a/b" string gives object_fields).

## Full rejectionReason, verbatim

The complete claim does not meet exploratory-v1. The reduction z=e(P,Q0)^x => Q=[x]Q0 is correct for the supplied P and nonidentity generator Q0. The code also supports the local repairs: integrated DP-first matching, recovery verification, and a launch-count clamp. However, material work and production-success objections survive.

1. Launch credit does not enforce the claimed computation budget. run_protocol bounds launches, not completed steps. Its test `t >= cfg.S and not inflight` is unreachable after the preceding empty-inflight check, and launches may continue after t exceeds S. The supplied trace has S=20000 but 654518 completed steps. The direct bound from launch count and per-trail cap is approximately 128S, not S plus the final in-flight population. Charging S*c_max and a 128*cap drain allowance therefore needs a stopping or probabilistic work argument that is absent. Adding an actual step cutoff would change the protocol and require reassessing success. The fixture trace is not a production lower bound, but it demonstrates that the asserted budget invariant is not implemented.

2. The delta=1/128 argument misuses Markov's inequality. Even assuming next-report distance D has mean theta, Markov bounds P(D>=128*theta). A follower misses when D exceeds its remaining lifetime R, where R<=128*theta; the former probability does not upper-bound P(D>R). The recorded ages do not repair this: protocol.py measures global completed work between insertion and matching, not a follower's residual report distance. The measured mean is 30891.7 with cap=4096, giving a ratio about 7.54, not 1/128; individual recorded ages also greatly exceed that cap.

3. The charged retention 6/10 is the fraction of distinct stored report keys subsequently revisited during one run with stopping on success disabled. It is not an estimate of detection probability conditional on a merge, which the sizing law needs. It depends on run length and store saturation. The small seeded success sweeps and their binomial intervals provide fixture-scope evidence only. Neither these measurements nor the disclosed shape extrapolation supports the production retention, conditional tail bound, or quantitative walk-quality allowance mu=1.1. In particular, theta is comparable to sqrt(N) in the fixture but vastly smaller in production. Disclosure and a cost debit do not establish that the measured statistic belongs in the sizing formula. Consequently success >=0.39 over full-production challenges and fresh coins remains unsupported; this objection does not demand a full-scale run or an unconditional mixing proof.

4. The original claim's exact minimality statement has a normalization conflict. The six ledger components are already target-operation equivalents, including pairing_row=19458889/5000. If W is their sum as stated, comparing W with 20000*2^c applies the normalization factor again, and 126.9558 is not the smallest four-decimal c satisfying that written inequality. Reinterpreting W as an instruction total would change the declared definition.

The self-generated certificate checks arithmetic consistency, not these missing premises. Thus neither the complete computation bound nor its required production success is adequately supported.

## Surviving flaws (my read)

1. BUDGET INVARIANT NOT IMPLEMENTED (hard). Credit bounds launches, not completed
   steps; `t >= cfg.S and not inflight` is dead code after the empty-inflight
   check; trace shows 654518 completed steps against S=20000, so work is ~128S.
   FIX: enforce a completed-step cutoff inside the round loop, re-run the success
   sweep under it (success will fall), and charge the true bound S instead of S
   plus a drain allowance.
2. MARKOV MISUSE (hard). P(D >= 128 theta) does not bound P(D > R), R <= 128 theta.
   The recorded ages are global completed work between insert and match, not the
   follower's residual distance. FIX: measure the follower's residual R at launch
   and the next-report distance D from the SAME trail, tabulate P(D > R) directly
   (no inequality), or derive the density for the shipped scheduler.
3. RETENTION STATISTIC IS THE WRONG OBJECT (hard). 6/10 is 'distinct stored keys
   revisited with stop_on_success=False'. FIX: P(detect | merge) measured with
   stopping on, per merge, over merges -- and mu=1.1 must be justified or dropped.
4. NORMALIZATION DOUBLE-COUNT (hard, cheap to fix). Rows are already
   target-operation equivalents (pairing_row = 19458889/5000), so W vs
   20000*2^c double-applies C_REF. FIX: either state W in instruction equivalents
   and divide once, or compare the summed target-op total directly.
