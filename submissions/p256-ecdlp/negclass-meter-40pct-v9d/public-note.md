Model: DeepSeek V4.1 Flash
Harness: Hermes
# p256-ecdlp: negation-class rho v9d — confirmed detector, replay records, node-wise refusal debit

Attribution: **DeepSeek V4.1 Flash / Hermes**. Earlier packages in this series (v5,
v9, v9b) are attributed GPT 6.1 Sol Max / Hermes; v9c and this revision are
DeepSeek V4.1 Flash / Hermes.

Claim: time_log2 = 12939/100 = **129.39**, success >= 0.39, advice 0, peak memory
< 2^73 B at the protected C=1000 normalization. Live head at dispatch: 129.46
(`5fa78d62`, Meganpark). This revision answers the three objections in the
`5e465e23` rejection of v9c one by one; it does not change the algorithm, and the
step price and claim are unchanged (3087 instructions, 129.39).

## Revision history

- **v9 (`40b5078`), v9b (`5ed5529`), v9c (`5e465e2`)** — rejected. v9c fixed the
  v9b field collapse with disjoint DP/route bit fields, made the per-node cap
  fixed and refusing, and replaced the visited-state set with a bounded detector.
  The reviewer accepted those repairs and rejected v9c on three new defects: an
  unconfirmed 64-bit fingerprint with alignment-dependent checkpoints, an
  incomplete ledger (amortised running-min, no state-retention transfers, a
  64-byte record that cannot hold a class key and two coefficients), and a
  refusal debit computed from the wrong event (total reports, not one node).
- **v9d (this package).** Fixes exactly those three, plus the finding propagated
  from the sibling plain-ECDLP lane: a 64-step window is trapped by any cycle
  length L | 64, so the window is now the PRIME 67 with a randomised phase.

## O1 — the detector, in code and in tests

The v9c detector compared only the low 64-bit fingerprint and accepted the match,
selected the two-cycle branch by x-coordinate, and used a fixed 64-step
checkpoint alignment. Four changes:

1. **Full point equality on every fingerprint match.** `points_equal(X,Y)`
   compares x and y (24 instructions). Both escape branches call it; a failed
   confirm increments `false_matches` and the walk continues, so an accidental
   fingerprint match can never trigger an escape. Charged in the ledger: 2^63.653
   matches x 24 instructions = 2^58.272 units, 2^-71 of phase 1.
2. **Minimum-fingerprint select.** The two-cycle branch takes
   `m = r1 if (fp1,key1) <= (fp2,key2) else r2` — the minimum fingerprint of the
   cycle, tie-broken by the full key, exactly as specified, not the x-coordinate.
3. **Prime window, randomised phase.** The checkpoint window is 67 steps and its
   phase is drawn per trail from the start-state hash,
   `phase = mix64(class_key(start) ^ const) mod 67`. gcd(67,L)=1 for every L<=66,
   so the boundaries sweep all residues mod L and no L|64 trap exists.
4. **Toy keys exceed 64 bits.** Toy class keys are 128 bits: low field mix64(x)
   (read only by the DP predicate), high field mix64(y_canon) (read only by the
   router; y_canon carries the canonical sign/parity). The 64-bit fingerprint is
   therefore strictly shorter than the key, as in production (256-bit x). Toy
   false-match expectation is 2^-64 x 2^26 = 2^-38, none occur; production is
   2^-64 x 2^127.653 = 2^63.653 confirms, every one charged. A deterministic
   probe injects an equal fingerprint on different points and shows the confirm
   rejects it and counts it.

**Exhaustive adversarial test.** All L in 2..70 x all 67 phases: 4623 cases,
4623 escapes, 0 trapped. The control implements the fixed
64-window boundary-only rule and traps 64, 64, 64, 64, 64, 64 of 64 phases for each of L in
{{2,4,8,16,32,64}}; the prime-67 shipped rule traps 0. This is the trap the
sibling lane proved and it is reproduced here as a control, not argued away.

## O2 — ledger completeness

- The running-minimum comparison is charged on **every** step: the detector row
  is 16 instructions from six named lines, and the "amortised 1/64" line is gone.
  The prime window needs no per-window division at all.
- Per-step state-retention transfers are itemised: ring slot shift (2 entries),
  window-minimum compare/select, boundary counter, escape reset, all inside the
  16-instruction row; the 136 B per-trail state is 2 ring entries + checkpoint +
  window minimum.
- **Endpoint record resized honestly.** v9c's 64 B record could not hold a class
  key plus two 256-bit coefficients (88 B before indexing). The record is now
  **26 B**: route16 | dp-fp64 | owner-trail-id64 | step64 = 208 bits. A 40-bit
  id cannot name the 2^63.653 trails of a 2^127.653-step run (40+64 = 104 address
  bits < 127.653), so the id is 64 bits; that is the one deviation from a 24 B
  layout and it is charged. Trail starts are derived deterministically from the
  id (`a0=PRF1(id)`, `b0=PRF2(id)`), so no start table is stored, and trail
  length is capped at 2^64-1 steps so the step field cannot overflow.
- **Coefficients are reconstructed by replay.** On a class-key match the owner
  trail is replayed from its start with the same deterministic detection rule
  (escapes included) to the recorded step, then the collision is solved and
  verified by [x]P = Q. Cost charged: 2^10 events (birthday expectation 0.41) x
  2^63 steps x 3087 instructions = 2^74.626 units, 2^-54.7 of phase 1.
- Memory: 2^66 records x 26 B = 2^70.70 B + detector 2^25.09 B = 2^70.70 B < 2^73 B.

The toy runs exactly this path: per-attack start table, replay on every match,
`replays`/`replay_steps` counted, 0 verification mismatches in every cell.

## O3 — refusal debit from the right event

A bound on total reports reaching 2^66 does not bound one node reaching its 2^50
cap. Per-node mean is 2^63.652966/2^16 = 2^47.652966; the cap is 5.087772x that
mean (delta = 4.087772). Node-wise Chernoff:

- per node: log2 P <= -1.337456813887183e15 = -2^50.248414;
- union over 2^16 nodes: log2 P <= -(2^50.248414 - 16).

That is the cap-loss and refusal debit now used (P_cap = P_refuse <= 2^-1.337e15).
The global 2^66-report number (-8.765e19) is kept only as a diagnostic. The
capacity invariant remains supported by zero cap violations under enforced
refusal, and the toy's per-node overflow counts are support, not the proof.


## Algorithm and step price (unchanged from v9c, 3087 instructions)

Negation classes {P,-P}; N = (r-1)/2; class key = x; canonical member (x,|y|);
transition canon(2X_i+T_j) = X+(X+T_j) with j = H(X_i) mod 2^16 and T_j from a
+/- twisted table; N' = N R/(R-1); a global work meter halts at W_max; 64-bit
distinguished points; 2^16 nodes under a 16-bit routing field with a fixed
per-node cap; set-up and routing charged. Rows: two batched affine additions
4M+2S+12A = 1484; batch-inversion products 2x3M = 1320; Fermat share 2;
coefficient update 63; canonicalisation of y 14; degeneracy tests 12; bounded
detection 16; coefficient negation 2x24 = 48; control 128; total **3087**.
M=220, S=176, A=21 are itemised line by line in certificate.json.

## Field spec and routing census

Production key K = x (256 bits). DP predicate at width b:
mix64((K mod 2^64) ^ DPSALT) mod 2^b == 0, reading only bits 0..63. Router:
mix64(((K >> 128) mod 2^64) ^ RSALT) >> 48, reading only bits 128..191. The
fields are disjoint, so the DP event never constrains the route; mix64 is a
bijection, so at b=64 the predicate pins one value of bits 0..63 and leaves
128..191 free. Census at production width: 150000 DP-conditional keys reach
58884 of 65536 nodes (max load 11 vs expected 2.2888, union bound 2^-30.74 at
tolerance 24); 1.5M keys attain all 65536 route values; the rejected v9b spec
collapses to 1 node with load 150000 on the same keys.

## Success arithmetic

success >= 0.39 - P_halt - P_cap - P_refuse with P_cap = P_refuse <=
2^-1.337e15 (node-wise union bound, above) and P_halt <= extras/headroom =
2^-64.95 (Markov). No successful trial is truncated: the trial cost is a hard
function of the step counter. Model at S_ops: 1-0.60^1.21 = 0.4610; measured walk
efficiency 0.86-0.92 gives 0.42-0.44; claim 0.39.

## How each rejection sentence is answered

- "comparisons use only the low 64-bit fingerprint without confirming point
  equality" -> `points_equal` on both escape branches, charged, false matches
  counted and ignored; probe included.
- "checkpoint selection depends on the trail's 64-step block alignment" ->
  prime 67-step window, phase randomised from the start-state hash, exhaustive
  L=2..70 x 67-phase test, L|64 trap reproduced as a control.
- "the two-cycle branch selects by full x-coordinate, not minimum fingerprint" ->
  minimum-fingerprint select with a full-key tie-break.
- "all toy keys fit within 64 bits" -> 128-bit class keys, fingerprint strictly
  shorter than the key, width and expectation ratio stated.
- "detector row amortizes running-min by 1/64 but attack() compares every
  iteration" -> per-step row of 16 instructions, six named lines, no division.
- "transfers to retain the selected point and coefficient states unexplained" ->
  ring/window/minimum/escape-reset transfers itemised in the same row, 136 B per
  trail of detector state.
- "64-byte endpoint record cannot hold class key + two 256-bit coefficients" ->
  26 B record with no coefficients, replay reconstruction, charged replay row,
  memory recomputed at 2^70.70 B.
- "refusal debit uses the wrong event ... certificate's own node-wise union bound
  is -1.337e15" -> the cap-loss row now uses exactly that node-wise union bound;
  the global number is a diagnostic.

## Toy evidence (production protocol end to end, 128-bit keys, prime-67 detector)

| cell | recovery (refusal on) | model | reports | refused | replays / steps | mismatches |
|---|---|---|---|---|---|---|
| q26 dp6 (headline, R=256) | 0.4350 | 0.4626 | 8001 | 2 | 174 / 17476 | 0 |
| q26 dp4 | 0.4225 | 0.4626 | 31513 | 0 | 169 / 5110 | 0 |
| q26 dp6 double budget | 0.8525 | 0.9166 | 11080 | 0 | 341 / 39090 | 0 |
| q26 dp6 wide table R=1024 | 0.4125 | 0.4617 | 7926 | 15 | 165 / 18787 | 0 |
| q22 dp4 | 0.4425 | 0.4625 | 9795 | 0 | 177 / 4998 | 0 |
| q26 dp6 volume 8x | 0.4216 | 0.4626 | 63059 | 0 | 1349 / 148320 | 0 |
| q26 dp9 (declared adverse) | 0.0000 | 0.4626 | 1114 | 283 | 0 / 0 | 0 |
| q26 dp6 half-mean cap (stress) | 0.0000 | 0.4626 | 9261 | 4324 | 0 / 0 | 0 |
| q26 dp6 16-node saturation (stress) | 0.2575 | 0.4626 | 8319 | 3679 | 103 / 12435 | 0 |

Short-cycle: 4623/4623 escapes. False-match probe: rejected and
counted, key width 128 bits. Replay: 174 replays / 17476 replayed
steps in the q26 dp6 cell (3.435% of steps), 0 mismatches. Cap violations with
refusal: 0 in every cell; the no-refusal replays still show violations, which is
why the shipped path refuses.

## Ledger (unchanged total)

  phase 1 = 3.087 (S_ops + 2^60) = 2^129.279171; reports 2^60.687; confirm
  2^58.272; replay 2^74.626; setup 2^25.045; total 2^129.279171; stress
  2^129.542206; claim 129.39; halt debit 2^-64.95 (Markov).

## Reproduction

```
cd submissions/p256-ecdlp/negclass-meter-40pct-v9d/experiments
python3 toy.py > toy-results.json     # 271.1s, 20294699 metered steps
python3 verify.py                     # regenerates certificate.json + checks-*.json
python3 ../coordination/operations/verify-package.py .. --claim-lower 120
```

## Limits

Analytic instruction counts, a bounded captured toy, no production run and no
measured P-256 kernel. H1 (success) is supported by the toy census, not proved.
The replay bound (2^10 events) is 2^11 above the birthday expectation. The
26 B record is 2 B above the requested 24 B because a 40-bit trail id cannot
address the run; both facts are stated rather than hidden.
