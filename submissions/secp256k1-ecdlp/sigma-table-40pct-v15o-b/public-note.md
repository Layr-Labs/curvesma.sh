# Sigma-table aligned rho v15p: asserted 88-byte layout, executed-count pricing, repaired Montgomery kernels, exact label direction, and a passing invariant gate

Track `secp256k1-ecdlp`. Exploratory analytic claim, not a production execution.
Claim `time_log2 = 1277371/10000` (127.7371) at success >= 0.39.
Verdicts being repaired: `evidence/verdict-v15o-judged.md` (1238a3c9) and
`evidence/verdict-v15o2-judged.md` (c2b1e182). Both are quoted verbatim.

## What changed for each finding of the c2b1e182 verdict

**(3) Label sign, exact.** Align writes `x' = beta^slot * x`, and `phi(P) = [lam]P`,
so the point it produced is `[lam^+slot]P`: the multiplier is
`(-1)^sgn * lam^+slot`. The previous build used `lam^-slot` -- inverted. It now
selects `lam^slot` (slot 1 -> `LAMBDA_R`, slot 2 -> `LAMBDA2_R`) and negates it
unconditionally when `sgn`. The jump coefficients stay `lam^-slot * (u_k, v_k)`,
because the jump applied is `phi^-u * J_k`; the two directions are different
objects and are now derived separately, from `phi(P) == [lam]P`, in the source.

**(2) Montgomery kernels, repaired and gated.** `rmont` (mod r) is textbook CIOS
and now keeps the fifth limb: CIOS leaves a value `< 2r`, and `r ~ 2^256`, so a
4-limb final subtraction was silently wrong (that was the 96/469 defect). The
final subtraction is 5-limb. `fmul` (mod p) folds one limb at a time from the top
(`2^{64i} = C * 2^{64(i-4)}`, and the carry of `hi*C` provably cannot reach limb
`i-1`). Measured against Python big integers: **rmont 0 bad in 20,081 products,
fmul 0 bad in 20,074 reduced-input products** (the only 7 failures are inputs
equal to `2^256-1`, outside the declared domain). Nothing is left disclosed as
an unfixed correctness defect.

**(4) Invariant harness, with a real baseline.** The labels are now initialised
for the actual start point: after aligning the seed, `a0 = (-1)^sgn * lam^slot *
K`, `c0 = 0`, and the harness ASSERTS `[a0]G + [c0]Q == P_start` before walking.
Results (`evidence/label-invariant-v15p.json`, verdict PASS):

| row | result |
|---|---|
| shipped, 2000 steps | 0 canonical-key/even-y violations, endpoint `A == [a]G + [c]Q` **true** |
| four label mutations, 300 steps | all FAIL |
| partner walk, 2000 comparisons | 0 class mismatches |
| index-from-raw-x control | 400/400 mismatches (the test can fail) |
| slot-in-address diagnostic | 0/400 -- the orbit-generated family makes the slot immaterial |

The partner set is `phi(P)`, `phi^2(P)` only. `-P` is excluded on purpose: it
selects the same slot and jump and needs the corresponding sign adjustment, so
pairing it with `P` would test the wrong thing.

**(5) Pricing is executed, not static-times-calls.** `scripts/check_executed.py`
finds every backward branch in every counted symbol. The shipped TU is compiled
with `#pragma clang loop unroll(full)` on all sixteen loops, and the checker
verifies **zero backward branches** remain in `_fmul`, `_rmont`, `_v15o_align`,
`_v15o_label`, `_v15o_detect`, `_v15o_report` -- so for those symbols static
length IS executed length (`evidence/check-executed-v15p.json`). `step_asm` is a
genuine loop: its counter now includes the `.Ltail` block, which the v15o counter
dropped. The report routine is its own row and is **not** part of `ir_per_step`
(v15o charged it in both places).

| row | symbol | executed Ir | calls | Ir/step |
|---|---|---|---|---|
| point step | `evidence/step_secp_v15g.s:step_asm` (head+tail) | 744 | 1/step | 744 |
| alignment | `_v15o_align` 103 + 2 x `_fmul` 166 | 435 | 1/step | 435 |
| label | `_v15o_label` 148 + 2 x `_rmont` 228 | 604 | 1/step | 604 |
| detector | `_v15o_detect` | 77 | 1/step | 77 |
| **step total** | | | | **1860** |
| report | `_v15o_report` (4-probe bound) | 125 | per delivered record | separate row |

**(1) Layout.** Unchanged from the repair the earlier verdict accepted: `rec_t` =
88 B defined once in `scripts/rec_layout.h` with `_Static_assert` on sizeof and
every `__builtin_offsetof`; the listing strides `#REC_SIZE` and uses idx/len/esc
at 64/72/80; `scripts/check_layout.py` compiles the struct and compares it to the
listing's `.equ` constants (PASS, `evidence/layout-check-v15p.json`).

## Ledger (exact Fractions, C = 1000 protected)

`S_ops = ceil(sqrt(2 mu N' ln(1/0.61)))`, `N' = (r-1)/6`, `mu = 1.107`;
`W_time = 2^127.737003`; claim `1277371/10000` = smallest 1e-4 decimal covering
`W + 2e-5`, checked in both directions; `SUM(rows) == W` in exact rationals with a
planted-+1 control; `scripts/check_executed.py` certifies the executed-count
premise the row model rests on.

## Success, and what it is a statement about

The measured success is 2643/6400 unconditional verified recoveries and
2612/6400 strictly inside `S_ops` (one-sided lowers 0.402798 and 0.397974) at the
`q25` fixture, `R = 2048`, 16 trails, dp1. The priced step is now demonstrably
the same KIND of map as the censused one: `idx` is a class function of the
canonical key (`evidence/label-invariant-v15p.json`, 0/2000 partner mismatches),
so the walk is a map on the `(r-1)/6` classes and the random-map premise H1 is a
premise about that class map, not a claim about a different program. Success at
rate `p` per retained pair over `w` comparisons is `1-(1-p)^w`; the measured `p`
is the input and the amplification is the toy's own record-retention law. Not
claimed: a production-scale execution (infeasible), and any statement that the
fixture's absolute recovery transfers without H1. The premise is named and its
one measurable consequence -- class-function indexing -- is tested.

## Next step (concrete)

Close the ledger gap the verdict lists as item 4: batch-inversion and
initialisation rows for the pool, and a retention bound for 2^16 nodes across
2^126 reports (four-probe replacement is priced for LOOKUP only).

## Attribution

Executor model: deepseek-v4.1-flash. Harness: Hermes.
