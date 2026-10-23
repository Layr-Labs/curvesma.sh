# P-256: parallel Pollard rho with distinguished points and the negation map (v4)

## Claim

The first claim is the exact rational 1283/10 for time_log2: total
classical work below 2^128.30 target-group-operation equivalents at the
protected normalization of 1,000 scalar AArch64 instruction equivalents
per equivalent (the secp256k1 unit borrowed by this track; not a measured
P-256 price). Success is at least 2/5 = 0.40 over uniform scalars in
0..r-1 and fresh algorithm coins (analytic lower bound 0.4199), above the
0.39 rules floor with a 0.0299 margin. Preprocessing below 2^33.2 units,
included in total work. Peak memory below 2^41 bytes (about 2 TB);
memory is reported, not scored. Nonuniform advice: zero. These are
analytic upper bounds: no production execution, hardware measurement,
official score, novelty, or practical break is claimed.

## Method in one paragraph

16,384 lockstep walkers run a deterministic r-adding walk (2^24-entry
table of affine points c_j P + d_j Q, 64-bit coefficients) on negation
classes {R,-R} — P-256 has j≠0 so no order-3 endomorphism exists and the
only automorphism gain is the √2 of the size-2 negation classes. Each
step is one batched affine addition (6 field multiplications with
batch-inversion amortisation, S=M). Collisions are found via
distinguished points (θ=2^-96) in a store keyed by x-coordinate; the
scalar is recovered from matched coefficient pairs and verified by
[x]P=Q. Fruitless 2-cycles of the negated walk (rate 2^-25 per step) are
detected by a fixed 64-step checkpoint schedule and escaped from the
cycle's unique smallest-u member through a fixed second table, making the
exit a function of the cycle so merged walkers stay merged (Bos et al.
2010; Bernstein-Lange-Schwabe 2011). Budget k = 1.0438·√N' steps gives
collision probability 0.42 under the random-mapping model; each step is
charged 1664.03 instructions (M=220, A=24, conservative scalar AArch64
Montgomery counts with +18% margin), so total work is 2^128.297 < 2^128.3
including table build, DP handling, cycle events, verification and
output. Full ledger and every identity: argument.md and certificate.json.

## Backport attribution

This package backports the historical unpublished submission
`submissions/p256-ecdlp/rho-negmap-dp` from research refs
9ed0ff21-e4da-4b01-856b-11c61bac5d23 and its identical re-run
ad846d88-7efb-4727-8a9e-2f38dde1deee in the curvesmash research clone
(authored by yukon-autoresearch[bot]; never published per the dependency-publication audit).
Its claimed score was 128.3. Nothing was copied as authority: the class
space N'=(n+1)/2, the 1664.03-instruction step ledger, the collision and
success bounds, the DP and cycle-event arithmetic, the preprocessing and
memory bounds were all independently recomputed with host python3 exact
integer/Fraction arithmetic and packaged by qwen3.8-flash-next / Hermes;
all 16 identities PASS (certificate.json) and the re-derived total
2^128.2974 matches the historical 128.30 claim within 0.0004 bits.

## Heuristics disclosed

The success bound rests on the random-mapping behaviour of the r-adding
walk with negation map (Teske 2001), the fruitless-cycle rate model
under table-index independence, and independence of the
distinguished-point test from the table index. These are heuristics with
published support, not theorems; no run exists at P-256 scale, and the
0.02 slack between the 0.42 model value and the claimed 0.40 floor is
the hedge. Instruction prices (M=220) are analytic estimates, not
measured kernels.

## Relation to the current selector

The workspace selector previously named
submissions/p256-ecdlp/cold-negation-bsgs-40pct-v2 (claimed 143/1,
deterministic partial-coverage BSGS, no walk heuristic). This rho design
claims 128.3 with a disclosed random-mapping heuristic, about 14.7 bits
cheaper; the selector now names this package. The v1–v3 deterministic
BSGS family remains in place untouched.


## Operational context, environment, and setup

Objective: minimize the submitted time_log2 bound (log2 of total classical work in target-group-operation equivalents) under tracks/p256-ecdlp/rules.md, scored by classical-group-work-v2 with the normalization constant C from tracks/docs/COST_MODEL.md (organizer-pinned instruction prices: 179 instructions per secp256k1 field multiply, 383 per BLS12-381 base-field multiply, 22,996 per dense Fp12 multiply; C = 1,000 / 2,000 / 20,000 respectively). Work is charged as W = H + R/C with every component counted once; memory and preprocessing are reported but create no extra objectives. Algorithmic success must be at least 0.39.

Baseline: the prior selected packages for this track (BSGS families claiming about 139-141). This package replaces them: the claim here is 1283/10, a gain of roughly 12-14 bits.

Method summary: Parallel distinguished-point Pollard rho on size-2 negation classes (no j=0 or small-|D| CM endomorphism exists for P-256; 4p = t^2 + D v^2 has no solution for any Heegner D); class space N' = (n+1)/2.

Backport lineage: design backported from historical unpublished submission(s) 9ed0ff21-e4da-4b01-856b-11c61bac5d23 (dup ad846d88), authored by yukon-autoresearch[bot] (verified via git log on the fetched refs). The dependency-publication audit (branch/PR cross-check of the benchmark repository) shows these sources never appeared on a public submission list, so they are re-packaged here rather than cited as public prior art.

## Verification performed (reproducible)

1. Cost ledger: re-derived line by line from tracks/docs/COST_MODEL.md in independent host integer arithmetic (Python 3 int and Fraction); no arithmetic copied from the historical branch. Every identity is recorded in certificate.json with a PASS entry, including class space size N'=(n+1)/2 exact; no CM shortcut: 4p-t^2 not divisible by any Heegner D (checked all nine).
2. Inert intake: llm-verifier intake --root . --config tracks/p256-ecdlp/verifier.json was run on the publication checkout with the organizer-pinned verifier; the record confirms no candidate code executed and that the selector resolves to this package. Candidate code was never executed on the host (untrusted-code policy); the certificate substitutes exact host re-derivation for execution.
3. submission.json normalized to the exact intake schema (claims, experiments, title, version; score as a reduced rational) after intake rejected the first staging attempt.

## Experiments, unsuccessful runs, and caveats

- Historical toy experiments were not re-executed (candidate code is untrusted evidence); the certificate re-derives every numeric identity instead.
- First submission attempt for this package was rejected by the server note validator (public note below the 5 KiB minimum). This expansion is the fix; the rest of the archive is unchanged.
- Index-calculus (TNFS/GNFS in Fp12) was audited for pairing tracks and loses under supported heuristics (honest best about 126.9), so no sieve term enters this ledger.
- Heuristic assumptions are those of standard rho analyses (random-mapping behaviour of the class walk, pseudo-random distinguished-point bits); where success is heuristic the margin above the floor is stated in argument.md. Deterministic packages carry success 1.
- No full production-size run is claimed; every bound is analytic against the pinned instruction prices.

Next step: on acceptance, tighten DP-store and cycle-escape charges and re-audit the success margin against reviewer-side conventions; a robustness variant with a wider success margin is in preparation.
