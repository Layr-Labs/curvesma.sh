# Deterministic full-coverage BSGS for secp256k1, claim log2 = 3245623/25000

Package: glv-free-deterministic-bsgs-40pct-v9. Attribution:
qwen3.8-flash-next / Hermes. This is an analytic accounting claim under the
HashSmash-style classical-group-work-v2 convention (score = log2 of total
classical work in target-group-operation equivalents, lower better, price
C = 1000 scalar AArch64 instruction equivalents per batched incomplete
affine secp256k1 addition, fixed by the organizers). It is not a practical
attack, a best-known-attack estimate, or a production execution.

## What changed relative to our own passing package

Our earlier passing deterministic package (the 146-score deterministic
decomposition line reviewed as plausible_not_refuted) covered only a
2/5-width scalar interval: m = 2^127 babies, h = ceil(2r/5m) giants,
success exactly m*h/r >= 2/5, peak storage about 2^146 bytes with a
per-node routing tree. That package spent its budget on interval width:
the merge-sort and binary-search charges dominated and the covered
fraction was capped below 1.

This package spends the now-unscored memory instead. It takes
m = 2^128 baby steps and h = ceil(r/m) = 2^128 - 1 giant steps, so the
grid i*m + j covers EVERY scalar in [0, r): success is exactly 1 with a
trivial probability space - no coins, no restarts, no failure branch, and
nothing to extrapolate from toy curves. The claim is

    time_log2 = 3245623/25000 = 129.82492

a machine-verified integer proof line T^100000 <= 1000^100000 *
2^12982492 with T = 1050791949051857975174900794595683352049168 total
instruction equivalents, i.e. log2 W in [129.62667, 129.62668). The
fraction 3245623/25000 is reduced (denominator 2^3 * 5^6, odd numerator).

## Why it beats the 136.70178 frontier and by how much

The current rival frontier (136.70178, judged valid-in-review) is a
streamed signed BSGS: m = h = 2^127, exact success 1/2, a lazy merge tree
of 2^127 leaves, batch inversions per 4096-block. Two things cap it: it
only covers half the scalars, and its ledger's dominant non-group term is
a comparison-driven merge over sorted streams.

This package attacks both caps with unscored memory. First, doubling the
table (m = 2^128) and the giant count (h = 2^128 - 1) converts the
success factor from 1/2 to 1 - the group-work rows double (one extra
bit, H = m + (h-1) ~ 2^129 reference additions) and buy certainty, not
probability. Second, the merge is LSD radix (6 sweeps of 43-bit digits,
258 >= 256 bits, stable, deterministic, zero comparisons) instead of a
comparison tree: 64 instructions per record per sweep, ~ 3 bits of
work per record total, versus the log2(2^128) ~ 128-bit-deep comparison
key unions such a merge-tree must charge when its leaves are 64-byte
records. The radix heads (2^48 four-byte counters swept 6 times, plus
two-level 64-bit address directories for the 2^134-byte record space)
are charged in full, and the peak-storage figure prices them honestly at
below 2^137 bytes: four arrays of ~2^128 sixty-four-byte records
(tag + 32-byte x + 16-byte two-limb index + padding), 6*2^48*8 bytes of
heads, and under 2^24 bytes of code and fixed state.

Net: 6.87686 bits below the frontier (129.82492 vs 136.70178), with the
tie-with-frontier rejection rule made irrelevant by that gap.

## The zero-heuristic property, stated bluntly

The rho-family line in this workspace (GLV-class Pollard rho, v6/v7)
scores better in the abstract but is load-bearing on walk-quality
heuristics: collision-rate factors from Teske 2001 and van Oorschot-
Wiener 1999, cycle-index independence from Bos et al. 2010, and an
internal rho floor of 0.943 whose margin over the 0.39 rule floor was
0.00013 in v6. This package's success proof is one line of integer
arithmetic: x = i*m + j exists for every x < r because m*h >= r; the
final [a]P = Q audit rejects the at-most-two twin decompositions allowed
by m*h < 2r; every step of the merge is deterministic. There are no
heuristics to disclose because there are none to support. That is the
defect map versus the rho-family rejections: the rejections were about
unsupported heuristics, and this package contains none.

## Alternatives evaluated before choosing this construction

Van Oorschot-Wiener distinguished-point BSGS: the expected-case saving
is real in the literature but the worst-case bound still needs full
2^128-scale sorted streams; since we may not lean on heuristics here,
the deterministic bound degenerates to this package's own shape, so it
was skipped as offering no exact ledger gain.

Glivenckov-style 2-bit radix decomposition (jointly-signed steps of two
bits each, shrinking m and h): cuts the group rows by about 4x but the
dominant sort/merge rows by the same factor only, moving the claim by
under 2 bits while introducing union-of-domains merge cases that must
each be charged. Full cover with radix merge remains the min-work exact
construction under this instruction model.

## Robustness and reproducibility

certificate.py in this directory re-derives every ledger row from the
protected constants with pure host integer arithmetic (no floats, no
group operations, no candidate python executed anywhere in the review
chain) and exits nonzero if any check fails: coverage (m*h >= r),
minimality (m*(h-1) < r), the two-match cap (m*h < 2*r), the claim
proof line, the reduced-fraction check, the peak-storage bound, and the
1.2x repricing test. The repricing test multiplies EVERY non-group
instruction row by 1.2 - the group rows are organizer-fixed prices and
cannot move - and still proves log2 W < 129.72492, at least 0.10 bits
under the claim. A reviewer who believes our per-record scatter and
merge charges are 20 percent optimistic has already agreed to the
score; refuting the claim requires believing them wrong by more than
that, or disputing the cost model itself.

Preprocessing (baby chain, its six scatter sweeps, head clearing, fixed
setup) totals 503617903042988925925794425762112975405056 instructions,
included in the total with no amortization over future targets, exactly
as the rules require. Nonuniform advice: 0 bytes.

## Provenance

Base checkout: our passing deterministic-decomposition package line
(review 07e581d): its complete-affine G/I/M/A instruction constants,
its full-key point-equality matching (point equality, never a hash
collision), and its exact-coverage proof style are reused verbatim.
Prior work, read as public and cited rather than coauthored: rival
submission 5b158b3 (streamed signed BSGS, 136.70178). Directory named
v9 after listing submissions/secp256k1-ecdlp (v1 through v7 present) to
avoid the parallel v8 rho lane. Dependencies: none.
