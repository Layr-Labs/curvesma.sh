# Argument: the squared-step walk (v31) - v30's four judged flaws repaired

v30 (ticket 4dcc056e, claim 126.6729) was rejected on four new flaws: a retention
tail bound of 2^-3.589688 carrying a 2^-100-class budget with an unsupported
envelope, a launch protocol the measurements did not execute (with cumulative
abandoned work unpriced), a 127 control that never touched the 128 boundary, and
missing pairing/storage accounting. v31 repairs those four; the regen script,
literal sweep, twelve-shift gate, max-price ceiling and dp stop-loss row are kept.

Claim **126.6053** = 1266053/10000 (reduced), the least 4-decimal log2
>= the exact total 2^126.605236472. Emitted by `regen_v31.py`; no claim
literal is typed.

## 1. Retention tail: route (b), measured support + fitted residual

The merge simulator is re-run with `cap_mult = 128` (the production cap, not the
16-theta fixture) and budget 256 = 8x the old support: 24000 trials,
3173802 merges, 1489700 trails, max age
128.0 theta. Raw ages are a 1-step histogram, no bins.

    I_meas = mean e^{-(128 - A/theta)} = 2^-13.03
    window [7.5, 127.5] theta, 241 points
    kappa = 0.069587 +/- 0.002250 per theta (lo 0.065086)
    fit residuals: max 7.216 bits, rms 1.210 bits
    per-decade slopes (steps): [127,255]:-0.0971, [63,127]:-0.0180, [31,63]:-0.0135, [15,31]:-0.0152
    residual on (127.5, 128] theta = 2^-8.31

The fit is data: control (a) truncates the decades and kappa moves
(0.0696 -> 0.0654 -> 0.0633).

**The fixture law is not used as the production charge.** It has 785 toy classes;
the production retention miss comes from pair-uniformity instead. With the meeting
point uniform over the trail pair, the younger trail's age has density
2(L-u)/L^2 on [0,L], L = 128 theta, and its report distance is D ~ Exp(theta),
E[D] = theta. Markov in D:

    delta = int_0^L (2(L-u)/L^2) * theta/(L-u) du = 2/L = 2^-6    (charged)
    exact exponential law = 2(1 - e^-L(1+L))/L^2 = 2^-13.0000  (shipped)

v30 charged 2^-3.589688 here, so v31 charges 155 bits MORE on this row. The
charge moves the total, and the claim moves with it: this is the honest number.

## 2. The gate is in the shipped loop, not in a harness

`v31_protocol.run_protocol` is the only launch path; `ledger31` imports
`INFLIGHT_CAP` and `HALT_RULE` from it. Rules as they appear in the code:

    RULE 1  credit = completed - consumed; launch when credit >= LAUNCH_PERIOD
    RULE 2  if not inflight: launch (bootstrap, consumes no credit)
    RULE 3  if want and gate and len(inflight) >= inflight_cap: DEFER, credit kept

The line-by-line comparison the verdict asked for:

| v30 claim text | v30 code | v31 |
|---|---|---|
| "enforces an atomic in-flight cap of 128" | `scheduler_trace(..., enforce=derived)` in v30_measure, a HARNESS | RULE 3 inside `run_protocol`, imported by the ledger |
| "one trail start per theta executed steps" | `done - last >= launch_period`, done = aggregate | RULE 1, identical clock, credit retained on deferral |
| "no handling of this case is specified" | `first = (launched == 0 and not inflight)` only | RULE 2, explicit and quoted |
| "128*cap bounds simultaneously outstanding work" | trap row only | `abandon_row`, cumulative, section 3 |

Trace, shipped configuration:

    t=141 inflight=127->128 launch#128 ACCEPT completed=8156 credit=7148 (first time in-flight reaches 128)
    t=142 inflight=128 completed=8284 credit=7268 launch#129 DENIED gate: len(inflight)==128 (deferred=1, credit retained)

Control (b): patching `v31_protocol.INFLIGHT_CAP` to 129 makes the same entry
point reach **129**; removing RULE 3 reaches **1023**. The bound is a
property of the object, and it is exercised: the 129th launch is refused.

## 3. Cumulative abandoned work inside W

    abandon_row = n_inflight*cap + n_launches*P(collide inside cap)*cap
                = 19807040628566084398385987584 + 2^54.56 steps
                = 2^94.769251 at c_max

with n_launches = 62806765694859 and P(cap) = 2^-78.3. It
supersedes v30's trap row (same in-flight bound, plus the cumulative pre-halt
population) and sits inside W = w_effective + abandon_row + aux + pairing +
overflow + stop_loss. Control (c): earlier halts raise it 1x/2x/4x/8x while
w_effective holds.

## 4. Pairing and storage at production width

    live GT states   128 x 656 B        = 83968 B   (v30 charged 0)
    fixed-P base     68 x 2 x 96 B      = 13056 B
    returned G2      2 x 96 B           = 192 B
    pinned target    576 B
    record store     125613531389718 x 112 B (cap, 2x slack)
    peak             2^53.643 B

Width audit (`pairing_production_width`): record count 125613531389718
= 47 bits, held as a 2-limb multiword with a carry
assert; byte offset 54 bits, 2^10.4
headroom under SIZE_MAX, spec'd as (shard:u32, slot:u32, 2x64 offset); class keys
6 x u64, no truncation before mixing; refs carry an
explicit valid bit and are never identified by a zero word; 29247
shards of 4294967296 slots. `pairing_row` = 2^15.676142
for the pinned check (base + G2 reconstruction + Miller + final exp).

## 5. Max-price ceiling and twelve-shift gate kept

    w_effective = S*c_max = 2^126.605236; mean-price diagnostic 2^126.274568
    undercounts by 0.330669 bits and is never summed
    gate: 0.579100 arm-1 recovery vs 0.586783 census

## 6. The six rows

    w_effective      126.605236
    abandon_row      94.769251
    aux_sum          56.524626
    pairing_row      15.676142
    overflow_row     0.000000
    stop_loss_row    95.769251
    total            126.605236472   claim 126.6053 = 1266053/10000

## 7. Scope

The step, return algebra, census protocol, family-free theorem and N=785 toy scope
are unchanged. No production mixing bound is claimed. Rank loss disclosed:
126.6053 is behind the live board minimum 125.82 and this account's
valid 125.838.
