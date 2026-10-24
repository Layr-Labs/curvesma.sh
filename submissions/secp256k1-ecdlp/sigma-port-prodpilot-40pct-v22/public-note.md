# secp256k1-ecdlp v22 — pilot at the shipped d=60, the parity root cause, and every literal at its artifact

**Claim: log2 W2 = 127.92945732, submitted 127.930.** The work arithmetic is unchanged and re-audited
by `ledger.py` (155 itemised rows, `sum(ROWS)==STEP_INSN` asserted, C=1000, sign/key at its itemised
15, w=2 exactly +1.000 bit, kappa product 2.344852037 <= 2.345734234). What changed is what the last
verdict asked for: a **pilot run at the shipped d=60 configuration**, the **selector parity bug**, and
the **literal drift**.

## (1) The pilot: shipped configuration, d=60, 2^20 walkers — and what it can and cannot measure

`scripts/pilot_prod_d60.py` → `evidence/prod-pilot-d60.json`. It runs the **shipped** configuration, not
a surrogate of it: the shipped `qwalk.step`, the shipped report gate at **d = 60**, **workers = 2^20**
lanes exactly, trials with fresh coins.

| row | measured |
|---|---|
| visits | **402,653,184** |
| DP fires at d=60 | **0** |
| DP-record collisions | **0** |
| sampled invariant violations (`a + c*x == canon`) | 0 |
| visit rate | **29,502,569.4 /s** |
| gate consistent with 2^-60 | **true** |

Then the arithmetic the rate makes decisive, in the pilot's own units: one DP record needs
2^60 = 1,152,921,504,606,846,976 visits; at the measured 2.95e7 visits/s that is
**3.91e10 s = 1,241 years for ONE record of ONE lane.** A d=60 success *measurement* is therefore not a
slow experiment, it is an impossible one: **every runnable d=60 pilot is 0-success by construction**, so a
zero cannot be evidence about the production event in either direction, and dressing one up as a
measurement would be a fabrication. The production event itself is a two-run collision at
**2^127.014 visits** with **2^67.014 DP records per run**, where the gate is saturated: the post-merge
trail is 2^106.72 steps, so P(no DP after a merge) = exp(-2^46.72). The d=60 gate is not a success
bottleneck at production, and it is fully charged in the ledger as the step's DP-mask row.

**The claim this package makes is the work cap W = 2^127.92945732.** The success row it can defend is
measured, and it is the measured-regime one: design arm = largest measured order of the transfer ladder,
d = 5, N = 2^33.42 classes, 1200 sets, **CP-95 lower 0.475029**, probability debits **5.700e-4**, giving
**0.474459 > 0.39 there**. Production d=60 success is **not claimed** and no transfer is asserted: the
verdict's ruling is accepted, not argued with.

## (2) Selector parity: root cause, exact fix, machine-verified row

The verdict is right, and the bug was **arithmetic**, not only prose. `selector_pair.py` wrote
`c = -(a1+a2) mod N`, then branched on `c % 2`: even `c` was solved as `c//2`; **odd `c` was counted as
"rotation channel" and discarded**. With N = 33560113 odd, 2 is invertible and
`inv2 = (N+1)//2 = 16780057`, so odd `c` is an ordinary sign case. Root cause: an `if c % 2 == 0`
shortcut written for even N, kept under an odd N, together with integer floor division `c // 2`.

v22 deletes the parity branch entirely and solves **every** mu6 channel exactly, verifying each candidate
on the curve (`[k]G == ±C`): `zeta=-1: x = c*inv2`; `zeta=±lam, ±lam^2: (1∓lam)C = (±lam*a2 - a1)G`. The
eigenvalue is identified **on the curve** (`[lam]G == (beta*G.x, G.y)`), taken by CRT of the roots mod
each prime factor of the composite toy order N = 3121 x 10753, and the run **asserts**
`lam^2 + lam + 1 == 0 mod N`.

**Machine-verified inverse row** (`parity_fix.inverse_row_20`, all 20 rows `equal: true`):
`x*2*inv2 mod N == x` for 20 random x, and `2*inv2 mod N == 1`.

Effect, measured, at the shipped 800-set arm: **sign channel 9, lambda channel 20 (now SOLVED, previously
discarded), trivial 497, verify_fail 0**; EC 29/800, hash 24/800. The fix is not cosmetic — the λ cases the old code
threw away are recoverable and are now recovered.

## (3) Literals: one number per artifact, and the arithmetic the verdict flagged

The arm shipped here is the **800-set** run in `evidence/selector-pair.json`; every literal in this note
and in `submission.json` is that file's own number, so the 1600-vs-1200 drift is gone by construction:

| literal | artifact |
|---|---|
| 800 sets, p = 33554467 (2^25), order N = 33560113 by CM, 16777233 classes, R = 64, W = 16, trails 258, d = 5 | `selector-pair.json` |
| EC **29/800** (0.03625), hash **24/800** (0.030) | `arms.ec`, `arms.hash` |
| delta **+0.00625**, 95% CI **[-0.01054, +0.02304]**, z **0.7296**, discordant **47** | `arms.*` |
| sign 9, lambda 20, trivial 497, verify_fail 0 | `counters.ec` |
| 402,653,184 visits, 0 DP fires, 29,502,569.4 visits/s | `prod-pilot-d60.json` |

The verdict's second arithmetic point is also right and is **named, not hidden**: `0.475029 - 0.002211 -
0.0002344 = 0.472584`, which is *not* the ledger's success row. The 2.211e-3 is the **fallback reserve in
reference operations** (2.21e18 ref), not a probability debit. The probability debits the ledger actually
subtracts are its seven named terms, summing to **5.700e-4** (including the 2.344e-4 restart loss), so
`0.475029 - 5.700e-4 = 0.474459` is the printed row. The old sentence mixed a ref count into a probability
sum; it is withdrawn.

## Honest limits

The pilot establishes the gate, the configuration and the rate; it cannot establish a production success
rate at d=60, and the arithmetic above says no runnable experiment can. The success floor is a
measured-regime statement at d=5 / 2^33.42 classes and is labelled as such. The selector arm is
exploratory with **no claim dependence**; its delta is inside its own CI. 127.930 does not rank (frontier
126.338, 1.5975 bits away). The goal of this round is a package whose claim is exactly as large as its
evidence and whose every integer is its artifact's.
