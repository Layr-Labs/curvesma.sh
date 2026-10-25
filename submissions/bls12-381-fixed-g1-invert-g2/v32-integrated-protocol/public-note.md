# v32: one integrated protocol, and what it honestly measured

## Objective

The challenge is fixed-argument pairing inversion: given a nonidentity P in G1
and z in GT, return Q in G2 with e(P,Q) = z. The package solves z = z0^x in the
order-r subgroup of GT, where z0 = e(P,Q0) for one fixed G2 point Q0, and returns
Q = [x]Q0, which satisfies the target by construction.

This revision repairs the four flaws that rejected v31 (ticket 816b4c54): a
protocol the evidence did not execute, free bootstrap launches, a sampler that
suppressed the report-sharing event it was supposed to measure, and output
construction left as "unchanged return algebra" at 80-bit label widths.

## Environment and setup

Python 3 standard library only; no third-party packages, no network, no GPU. The
tower arithmetic, the class key and the twelve-shift solve are the ones shipped
in scripts/ (unchanged from the v27 package the reviewer examined). The fixture
is the 785-class order-12 tower with r = 9421, p = 1152921504607943563 and
N = (r-1)/12 = 785.

    cd submissions/bls12-381-fixed-g1-invert-g2/v32-integrated-protocol/scripts
    python3 measure_v32.py          # writes ../evidence/{protocol,gate,retention,audit}32.json
    python3 regen_v32.py --write    # writes ../evidence/ledger32.json + ../certificate.json
    python3 regen_v32.py --check    # re-derives the claim, asserts agreement
    python3 regen_v32.py --control-ledger   # control (d)

## Prior work and what was wrong with it

v31 claimed 126.6053. Its `run_protocol` decremented fixed integer trail lengths
and never called the walk, the DP test, the report match or the recovery; the
success and tail experiments used separate replacement loops. Its RULE 2
launched a bootstrap trail without consuming credit, so `ceil(S/theta)` was not
an enforced launch bound. Its `v31_tail.batch` executed `continue` whenever a key
was already stored *before* handling `is_dp`, so a repeated distinguished point
did not terminate its trail, the shared-report event was suppressed, and repeat
visits were counted as evidence. Its retention law (age density 2(L-u)/L^2 with a
conditional mean report distance theta) was assumed, and the fixture estimate was
floored at 1/64 rather than charged. Its output path described an affine
square-root chain rather than Q = [x]Q0, priced powers of two without instruction
counts, described labels as 80-bit/two-word although the step updates them modulo
the 255-bit order, and gave the own-step counter 8 B against a 2^87-step cap. Its
control d printed the same claim for the base and the perturbed ledger.

## Implementation

`scripts/protocol.py` is now the only entry point. `run_protocol` executes, per
trail and in this order: `classkey.step_v26` (the family-free squared step), the
distinguished-point test, the store insert-or-match on the class key, and the
twelve-shift recovery with verification against the planted target. The success
sweep, the gate trace, the retention histogram and all four audits call it; no
replacement loop and no harness-side scheduler exists in the package.

Credit rules, enforced on one counter: `credit = completed - consumed`; a launch
needs `credit >= theta`; a bootstrap launch CONSUMES theta (the v31 defect); a
credit-eligible launch at INFLIGHT_CAP = 128 is deferred and keeps its credit;
total launches are clamped at `launch_bound = ceil(S/theta)`, bootstrap included;
the protocol halts at the first verified recovery or at the bound.

## Measured results

Success sweep through `run_protocol`, 120 planted challenges per budget,
theta = 2^5, cap = 128 theta (evidence/protocol32.json):

    budget   recovered   CP95 lower   mean launches   bound   bound binds
       400     105/120      0.8141           6.2        13         20
       800     114/120      0.9037           7.5        25          6
      1600     120/120      0.9753           8.3        50          2
      3200     120/120      0.9753           8.3       100          0

Gate and launch-bound trace (evidence/gate32.json): round 1 bootstrap consumes 32
steps and starts at credit -32; in-flight first reaches 128 at round 848
(t = 35710, launches = 485); the next credit-eligible launch is deferred with
`reason='len(inflight)==128'` at round 850; the run halts at exactly 625 of 625
launches with 3628 deferrals and 654518 completed steps.

Retention on ONE identical input, shipped order against the v31 order
(evidence/retention32.json):

    order                  reports  shared  matches  retention  repeat visits
    is_dp before store          10       6       40      0.600          0
    dedup before is_dp          10       0        0      0.000       2770

Audits (evidence/audit32.json), all four planted:

    (a) integration   success through the protocol 0.75 -> 0.00 with the DP test
                      removed, same 24 challenges: the experiments ride the
                      protocol because there is nothing else to call.
    (b) sampler order retention 0.600 vs 0.000, 2770 suppressed shared-report
                      events, histograms differ (9 occupied age bins vs 0).
    (c) bootstrap     trails cut at theta/4: launches 125 = bound 125, binds true,
                      consumed = 125 * 32 = 4000 exactly (v31 charged 0 for the
                      same 125 launches and enforced no launch counter at all).
    (d) regen carry   base 126.9558; total x(1+1e-4) -> 126.9559; retention
                      3/5 -> 599/1000 -> 126.9570. The claim moves.

Unsuccessful runs are reported, not hidden: the first sweep used 320 trials and a
`skip_dp` audit budget of 2000 steps, which is slow because a trail with no DP
test runs to its cap, so the shipped configuration uses 120 trials and a 400-step
audit budget; the numbers above are the shipped ones. An earlier retention
statistic divided matches by reports and returned 4.0, which is a visit count and
not a probability; it was replaced by the count of reports shared at least once
(6 of 10).

## Claim

    w_effective      126.955711     S * c_max, c_max = 1.704385 (max price)
    abandon_row       35.683584     launches x P(merge inside cap) x cap
    aux_sum           56.875101     report insert/match traffic, sorting
    pairing_row       11.926214     Miller 1533 + final exp 300 + G2 [x]Q0 226
                                    + verify 1833 Fp12 mults = 77835556 instr
    overflow_row       0            reports capped at ceil(S/theta) literally
    stop_loss_row     94.769251     in-flight population + abandoned work
    total            126.955711     claim 126.9558 = 634779/5000

S solves (1-delta)(1-exp(-S^2(1-f^2)/(2 mu N))) >= 0.39 with mu = 11/10,
f = 3/8 measured on the fixture, N = 2^251.2721, delta = 1/128 by Markov on the
follower's remaining length, and the MEASURED first-report retention 0.600 in the
slot the law divides by. The claim is 0.35 bits worse than v31 because it charges
the retention it measured instead of declaring 0.9854.

## Limitations

The experiments run at fixture scope: N = 785 and theta = 2^5 against N = 2^251
and theta = 2^80. They establish that the shipped algorithm executes, that the
gate and the launch bound bind, and that the sampler order changes the retention
statistic by the stated amount. They are not a production recovery rate and no
production run is claimed. The transfer of the retention shape to production is a
disclosed extrapolation in the class the track permits. The pairing row is
instruction-counted from published tower formulas, not measured on hardware. No
production mixing bound is claimed, and the claim sits behind the live board
minimum on this track.

## Provenance

Every number above is emitted by scripts/regen_v32.py from evidence/*.json, which
scripts/measure_v32.py writes. No claim literal is typed anywhere; `--check`
re-derives the claim and asserts the ledger agrees; `--control-ledger` shows the
claim moving under two independent perturbations. The v31 rejection text is kept
verbatim at submissions/bls12-381-fixed-g1-invert-g2/verdict-v31-g1inv.md, outside
this package.
