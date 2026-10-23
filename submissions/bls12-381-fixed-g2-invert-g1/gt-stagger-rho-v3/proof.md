# Fixed-argument pairing inversion on BLS12-381 G2: staggered-checkpoint rho in GT

Track `bls12-381-fixed-g2-invert-g1`. Candidate `gt-stagger-rho-v3`.
Claim **127.68** (`12768/100`) = log2 of total work in reference dense-`Fp12`
multiplication equivalents (`C = 20,000` instructions). Success at least **0.39**
over uniform nonidentity `Q`, uniform `z` and fresh algorithm coins, as
`rules.md` requires. Peak memory **227,416 B**.

The predecessor ids cited below (`59149222`, `a7ac78ea`, `4c234c5c`) are submissions on the
sibling track `bls12-381-fixed-g1-invert-g2`. This package is their role-swapped port (the
challenge input is `Q` in G2 and the returned point is `P` in G1); the reviews of those
predecessors apply to it verbatim, and the port changes no cost, no schedule and no claim.

`v3` is `v2` (`a7ac78ea`, rejected) with the review's three items repaired, and
with each repair *measured* rather than asserted: the toy now runs the production
coefficient regime (`cap*2^cw < q`, zero gap-sum wraps), the storage peak is the
honest 13 records (12 compared), and the filter carries a hard `2^70`-comparison
allowance that aborts instead of an expected cost. §8.1 maps the three
objections to the three fixes.

Every number below is recomputed by `experiment.py` (standard library only) from
the declared inputs; `results.json` is the captured run of that file. The
rejected predecessor `59149222` (`gt-stagger-rho-v1`) is superseded, not patched
silently: §8 lists each objection it earned and where the fix lives.

## 1. Reduction

Instance: nonidentity `Q` in G2 and `z` in GT, `r = |G1| = |G2| = |GT|` prime.
Fix a public nonidentity `P0` in G1 (the pinned generator), and compute

    g = e(P0, Q) = pairing(Q, P0)        (the pinned library takes (Q, P) order)

`e` is nondegenerate and bilinear, `r` is prime and `P != O`, so `g != 1` and
`<g> = GT` (a prime-order group has no proper nontrivial subgroup). Hence a
unique `a` in `Z/r` has `z = g^a`, and `e([a]P0, Q) = e(P0,Q)^a = z`, which
solves the instance, including `z = 1` (`a = 0`, `P = O`). What is left to price
is one discrete logarithm in the target group, with no auxiliary powers.

`experiment.py` §A re-checks the production facts on BLS12-381 itself with
integer arithmetic (`p`, `r` prime; `r | Phi_12(p)`; `ord_r(p) = 12` with
`p^6 = -1`; the blst tower is a field; a random GT element has order `r`; the
collision identity of §5 recovers a real 381-bit logarithm).

## 2. The walk

Public jump table `Th[h] = g^(A_h) z^(B_h)`, `h = 0..R-1`, `R = 256`, built once
per instance from `(g, z)`. State `(X, a, b)` with `X = g^a z^b`, started from a
fresh uniform `(a0, b0)`:

    h  = H(X)           one 64-bit limb of the canonical Fp12 encoding, fixed odd
                        multiplier, xor-shift, top 8 bits
    X <- X * Th[h]      one dense Fp12 multiplication = 1 reference unit
    a <- a + A_h mod r,  b <- b + B_h mod r        (A_h, B_h in [1, 2^64))

The next state is a function of the group element alone, so the trajectory is a
deterministic map `X -> X * Th[H(X)]` on a set of size `r`, modelled as a random
map on a cyclic group of prime order `r`. There is no history, no value lookup,
no distinguished-point set and no dictionary; the only sequence data kept is the
twelve-entry comparison window of §3 (thirteen stored at peak).

`A_h, B_h in [1, 2^64)` is not decoration: it makes every gap exponent a
strictly positive integer below `2^192 < r`, which is what removes the `den = 0`
degeneracy in §5 unconditionally rather than with probability `1 - 1/r`.

The toy run of §5 runs the **same restriction**, not the whole of `Z_q`: its
coefficients are drawn from `[1, 2^cw)` with `cw = bitlen(q) - bitlen(cap) - 1`,
so `cap*2^cw < q` and every toy gap sum is likewise a positive integer below the
group order. `cw = 12` at the declared inputs, and the toy's `wrap_violations`
counter (true gap exponent outside `(0, q)`, or disagreeing with the reduced
difference) is **0** over the 600 declared trials.

The table is checked at build time for `Th[h] != 1` (possible only if
`A_h + B_h*a = 0 mod r`, probability about `2^-255` per entry); a hit would leave
the walk in one step and the same solve formula still applies, so it is charged
and not relied on. Cost of the table: `2*R = 512` GT scalar multiplications, each
about 383 dense multiplications (`(1<<13.6)` units), plus one pairing for `g`.

## 3. Detection: staggered checkpoints, and the working set is *counted*

No predecessor text survives here. The rejected candidates charged `2*rho+3` or
`2*S0` evaluations; a tail-0 cycle of length 9 is detected at evaluation 24, so
that sentence was false and it is deleted. This candidate charges **no
cycle-length bound at all**: it charges a hard cap (§4), and the cap is a
property of the run, not of the trajectory.

The detector is a staggered checkpoint schedule

    P_k = { ceil(2^j (k+i)/k) : j >= 0, 0 <= i < k },      k = 4

so `P_4` carries exactly 4 points in every dyadic level `[2^j, 2^(j+1))`. A
checkpoint created at position `p` stores `(X, fp(X), a, b, p)` and is compared
with the hare only at hare positions in the **time window** `[p+1, W*p]`,
`W = 8`, and is then dropped. Expiry is one test on the oldest entry, because the
live list is in increasing `p`.

**Working-set law.** A window of multiplicative span `W = 2^b` straddles exactly
`b` dyadic levels and each level carries `k` points, so

    n_compared(k, W) = k * log2(W)        (W a power of two);   (4,8) -> 12
    n_stored(k, W)   = k * log2(W) + 1    (the record being inserted)  -> 13

This is the correction that matters, and it is the predecessor's fatal error in
its 192-versus-56 form: `k*(1-1/W)` is not merely a bad constant, it has the
wrong form. The `+1` is the review's second item on `a7ac78ea`: the record
created at the current position is appended *after* the comparison, so at hare
position 32 the twelve live records and the new one coexist — thirteen stored
records, of which twelve are ever compared. `experiment.py` §B1 does not assert
either number, it **enumerates** them: it walks every hare position in
`[1, 2^20]`, records the maximum compared (12) and the maximum **stored** (13,
first at `h = 32`), and prints the compared count at the probe positions 1024,
4096 and 2^20 (12, 12, 12). The old `k*(1-1/W)` value is printed beside the
enumerated one so the two cannot be confused again, and the toy's own
`observed_max_live` (13) and `observed_max_compared` (12) are measured on the
running walk, not on the schedule in the abstract.

**Detection index.** For a trajectory with tail `mu` and cycle `lam`, positions
`c < m` carry the same element iff `m = c`, or `c >= mu` and `m = c mod lam`. A
checkpoint stored at `p` and compared with the hare at `p+1 .. W*p` is therefore
matched first at

    m = p* + lam,  p* = min{ p in P_k : p >= max(mu, ceil(lam/(W-1))) },

because the checkpoint must lie on the cycle (`p >= mu`) and its window must
still be open when the hare reaches `p + lam` (`lam <= (W-1)p`). `m` is the
number of map evaluations, hence the evaluation count of the run.

Verified at two independent scopes: inside `experiment.py` §B2, against a
**direct simulation** of the comparison decisions (build the grid, advance the
hare, test the live checkpoints) for every `(mu, lam)` with `mu + lam <= 60`
(1830 pairs, **0 mismatches**); and offline, exhaustively over every pair with
`mu <= 600, mu+lam <= 600` for `k=1,W=2` (180,300 pairs) and over 9 `(k,W)`
schedules x 4,095 pairs — **0 mismatches**.

## 4. Success probability: a computed quantile, and a hard cap

Model, disclosed as a heuristic: `rho = mu + lam` satisfies
`Pr[rho >= t] = exp(-t^2/(2r))`, and given `rho`, `mu` is uniform on `[0, rho)`.
This is the standard random-map (Flajolet-Odlyzko) model for rho. Pushing the
**exact** detection law of §3 through it by deterministic quadrature in
`(rho, mu)` — no sampling error — gives the quantiles of `m/sqrt(r)` for
`(k,W) = (4,8)`:

| quantile | 0.39 | 0.42 | 0.45 | 0.50 | mean |
| --- | --- | --- | --- | --- | --- |
| `c = Q_q(m)/sqrt(r)` | 1.0507 | 1.1046 | 1.1569 | 1.2432 | 1.3255 |

`experiment.py` **§D1** is that quadrature (600 points in `rho/sqrt(r)` x 200 in
`mu/rho`, truncated mass 0.99998, renormalised), and it is implemented in the
declared experiment, not quoted from elsewhere. **§D2** is an independent
Monte-Carlo cross-check of the same model, 20,000 draws: `c(0.39) = 1.0453`,
`c(0.45) = 1.1524`, i.e. agreement to `0.005` in the ratio, which is the
sampling error of the cross-check itself.

**Dyadic alignment (formerly "phase").** The grid is dyadic, so where
`sqrt(r)` sits between two powers of two shifts the rounding cost. `§D3` sweeps
16 alignments `sqrt(N_t) = sqrt(r)*2^(t/16)`, one full period:

    c(0.39) in [1.0463, 1.0527]     spread 0.61%
    c(0.45) max 1.1591              success at the cap, min over alignments 0.4501

The submission does **not** claim a per-alignment statement — `rules.md` asks for
success over uniform `Q`, uniform `z` and fresh coins. The sweep is here as
supporting evidence that the declared cap is not an artefact of the lucky
alignment of the production scale: the worst alignment still has its 0.39
quantile at 1.0527, far below the cap.

**Hard cap.** `B = ceil(1.16*sqrt(r)) = 265627019777073510555230149274057194865`
(exact integer square root; `B-1` verifiably fails). The run stops after `B`
evaluations and reports FAILURE, so the work of a *failed* run is bounded by the
same `B` as a successful one. With `c = 1.16`,

    P[success] = P[m <= B] = 0.4518 at the nominal alignment, 0.4501 worst case
    claimed: success >= 39/100  (rules.md threshold, 0.06 absolute margin)

No amplification or restart is charged, because none is needed: one trial above
0.39 is claimed and the trial above is at 0.45.

**Toy validation, at the production table size.** `§B3` runs the complete
algorithm — same grid, same window, same expiry, same solve, `R = 256` as
declared, and in the production coefficient regime of §2 — on a 26-bit
prime-order subgroup, 600 trials, cap `ceil(1.16*sqrt(q)) = 6720`:

| quantity | value |
| --- | --- |
| solved inside the cap | 272/600 = **0.4533** (Wilson 95% [0.4139, 0.4933]) |
| solved + capped = trials | 600 = 600 |
| recovered logarithm verified correct | 272/272 |
| bookkeeping failures | 0 |
| live-bound violations / max stored / max compared | 0 / **13** / **12** |
| degenerate relations | 0 |
| gap-sum wrap violations (true gap outside `(0,q)`) | 0 |

The model predicts 0.4518 at this cap; the toy gives 0.4533, inside the interval
and closer to the model than the old whole-of-`Z_q` toy's 0.4967.
`§B4` repeats 120 trials at `R = 1024`: 65/120 = 0.5417 (Wilson [0.4526, 0.6281]).
So the success fraction is not an artefact of the declared table size — raising
`R` costs memory only, it does not change the per-step cost or the detection law.
This is the check the earlier "experiments use R=64 rather than production R=1024"
objection asked for: the declared `R` **is** the `R` used for the main trial set.
Toy statistics are evidence at toy scope only.

## 5. Solve, degenerate draws, and the withdrawn argument

On a detected match the two states are equal group elements with exponent pairs
`(a,b)` and `(ta,tb)`, so `(ta-a) + (tb-b)a_full = 0 mod r`:

    den = (b - tb) mod r;   if den = 0 the draw is restarted inside the same cap
    a_full = (ta - a) * den^-1 mod r

exact because the elements are equal. `den = 0` requires the product of the
jumps over the gap to be the identity: `sum A_h = 0 AND sum B_h = 0 mod r`.
Because every `A_h, B_h` is in `[1, 2^64)` and the gap is shorter than
`B < 2^128`, both sums over any gap lie in `(0, 2^192)` with `2^192 < r`, so
**neither can vanish** — `den != 0` unconditionally, and the `1/r` degeneracy
argument is not needed. `experiment.py` §A6 checks the related single-step case
(`Th[h] = 1`) on real BLS12-381 elements.

The toy does not take this on faith: every solve asserts the *true* gap exponent
(unreduced `bi - tbi`) lies in `(0, q)` and agrees with `(b - tb) mod q`, and
counts violations. Over the 600 declared trials: **0 violations, 0 degenerate
relations**. The lemma that `den != 0` is therefore measured in the same
coefficient regime the production claim uses, which is exactly what the review of
`a7ac78ea` asked for.

The toy's two unreduced accumulators (`ai`, `bi`) are *measurement scaffolding*:
they exist only to test the no-wrap property, they are never part of the priced
680-byte record, and production does not carry them — there the lemma is a
proof, here it is an experiment. The priced record remains
`(X, fp(X), a, b, p)` = 576+8+32+32+32 = 680 bytes.

The earlier argument of this account (`T_i^2 = 1`, "two-step returns are
impossible") is **withdrawn**: a return needs `T_i * T_j = 1`, and elements of
odd prime order can be inverses of one another. The positivity lemma above
replaces it.

No quotient by the order-12 Frobenius class and no canonicalisation is used, so
the structured two-cycle family raised against quotient walks on this track does
not exist here.

## 6. Cost ledger (exact)

Units: one dense `Fp12` multiplication at the organizer's fixed `C = 20,000`
instructions. `r = 2^254.857089`, `sqrt(r) = 2^127.428545`.

Per evaluation, charged `512` instructions on an itemised sum of `116`
(safety factor 4.41):

| item | instructions |
| --- | --- |
| jump index: load limb, multiply by the fixed odd constant, xor-shift, mask | 8 |
| table addressing and coefficient loads (`A_h`, `B_h`) | 4 |
| `a += A_h mod r` (4-limb add, conditional subtract) | 11 |
| `b += B_h mod r` | 11 |
| fingerprint of the new element | 2 |
| 12 live fingerprint comparisons (load, cmp, branch, pointer) | 48 |
| ring insertion and expiry, amortised | 24 |
| loop control and cap test | 8 |
| **itemised sum** | **116** |
| **charged per evaluation** | **512** |

so `step_units = 1 + 512/20000 = 641/625 = 2^0.036468`. The reference dense
`Fp12` multiplication is 1.0 unit by definition and is not re-priced.

**Fingerprint discipline, stated once and used everywhere.** The stored
fingerprint is `fp(X)`, the low 64 bits of the canonical encoding of `X`. It is
an exact *necessary* condition for equality (equal elements have equal
encodings), so the filter cannot miss a true match; it can only produce false
positives, and each is settled by one **full 576-byte element comparison** of the
stored `X` — never by a reconstruction from exponents, and never by trusting the
fingerprint. There is no 16-bit mask anywhere in this candidate: the width is 64
bits in the walk, in the ledger and in the experiment. The run performs at most
`12*B = 2^131.3` fingerprint comparisons. Treating a 64-bit fingerprint as
uniform gives an *expected* `2^67.3` matches, each settled by one full 576-byte
element comparison at 216 instructions (72 load+cmp+branch pairs) —
`2^74.98` instructions `= 2^60.69` units.

The **number** of matches is a random variable, so an expected cost is not a
work bound (the review's third item on `a7ac78ea`). The run instead carries a
**hard allowance of `2^70` full element comparisons and aborts** if it is
exceeded:

    charge = 2^70 * 216 / 20000 = 2^63.46 units, deterministic for every run

The allowance is `2^(70-67.3) = 6.5x` the expected number of matches, so the
abort is a bound and never an event on any realistic run; and it is 64 bits below
the main term, so it moves the claim by less than 0.001 bits. The walk in
`experiment.py` implements the same structure: it counts `full_compares` and
returns `fp_abort` when the allowance is passed, and `experiment.py` §B3 reports
`fp_aborts = 0` at the declared inputs. The toy group of §B3 is 26 bits wide, so
its fingerprint is injective there and the toy cannot even produce a false
positive: the allowance is carried because production elements are 576 bytes and
the filter is 8 of them.

Additive work:

| item | units |
| --- | --- |
| jump table, `2*R = 512` GT scalar multiplications (196,096) | 2^17.58 |
| one pairing for `g = e(P0,Q)` | 131,072 = 2^17 |
| verification `e([a]P0,Q) = z` and `[a]P0` | 135,168 ≈ 2^17.04 |
| hard full-comparison allowance (`2^70` compares x 216 instructions, abort) | 2^63.47 |
| **total additive** | **2^63.47** |

Total work

    W = B * 641/625 + additive = 2^127.679138,   claim = 12768/100 = 127.68

The additive term is `2^63.47` against a main term of `2^127.68`; it moves the
claim by less than `2^-64` bits and is carried anyway.

**Peak memory in bytes**, itemised by the structure that holds them
(`experiment.py` §C computes exactly these terms):

| item | bytes |
| --- | --- |
| jump table: 256 x (576-byte dense `Fp12` element + 8-byte `A_h` + 8-byte `B_h`) | 151,552 |
| stored checkpoints: 13 x (576 `X` + 8 fingerprint + 32 `a` + 32 `b` + 32 position) | 8,840 |
| hare state: 576 element + 8 fingerprint + 32 `a` + 32 `b` + 8 step | 656 |
| `g = e(P0,Q)`, the table's base | 576 |
| control: power, `lam`, counters, RNG state, pointers | 256 |
| code, stack and scratch for the walk and the detection loop | 65,536 |
| **peak total** | **227,416** |

Notes on the three record widths that the predecessor got wrong: `a` and `b` are
32 bytes (four 64-bit limbs, exact for a value mod `r`, 255 bits); the grid
position is 32 bytes, because the schedule runs to `B > 2^127` and an 8-byte
field could not address it; and `X` is **stored** at 576 bytes per checkpoint, so
no comparison ever needs a reconstruction. Nothing grows with walk length, the
number of recorded points, or the number of trials.

## 7. Heuristics disclosed

1. The random-map model for `X -> X * Th[H(X)]` on GT, with the Rayleigh
   collision time. Standard; supported by §B3/§B4 at production `R`, in the
   production coefficient regime: 0.4533 (272/600, Wilson [0.4139, 0.4933]) and
   0.5417 (65/120, Wilson [0.4526, 0.6281]) at `R = 1024`, against a predicted
   0.4518, and by the exactness of the detection index given `(mu, lam)`.
2. The uniform tail `mu`. The toy's mean detection index over solved trials
   (0.73 sqrt(q)) is where the model puts the conditional mean.
3. The per-item instruction counts of §6. The claim is insensitive to them:
   dropping the 12-comparison line entirely moves the claim by 0.0034 bits
   (to 127.6758), and charging a full 576-byte comparison against all 12 states
   every step instead of one fingerprint moves it to 127.8478 (0.1687 bits).
   With the non-reference overhead priced at zero the claim would be 127.6427,
   which is the floor of this schedule.
4. The success probability is a model quantile, not a theorem. That is why the
   cap sits at the model's 0.45 quantile while only the 0.39 floor is claimed,
   and why the 0.39 quantile is checked at every dyadic alignment.

## 8. Answers to the rejection of `59149222`

| objection | where the fix lives |
| --- | --- |
| checkpoints store `X` in the experiment but the claimed 112-byte record omits it, and a 120-instruction comparison on a prefix hit did not pay for reconstructing it | §3 and §6: `X` is **stored** (576 B), the record is 680 B, and the full comparison is charged at 216 instructions on a fingerprint hit. No reconstruction is used anywhere, so nothing uncharged is needed. |
| the 589,824-byte jump table held only the group elements, not `A_h`/`B_h`, and four limbs per exponent consumed the whole 64 KiB allowance | §6: the table is `R = 256` entries of `(element, 8-byte A_h, 8-byte B_h)` = 151,552 B. Coefficients are below `2^64` by §2, so one limb each is exact and the positivity lemma still holds (`2^192 < r`). |
| the 8-byte position field cannot address a `2^128`-step schedule | §6: 32 bytes, and the ledger lists it. |
| `proof.md` specified a 16-bit mask while the ledger used `2^-64`, understating expected hits by `2^48` | §6 and §3: there is no 16-bit mask. Fingerprints are 64 bits in the walk, in the ledger and in `experiment.py`; the expected false-positive count, its per-hit cost and the charge all use 64. |
| the quadrature, Monte-Carlo cross-check and phase sweep were asserted, not implemented in the declared experiment | §4: §D1 (quadrature), §D2 (20,000-draw Monte-Carlo), §D3 (16-alignment sweep), §D4 (success at the cap) are all in `experiment.py`. |
| "success 0.45 for EVERY alignment" was only supported by 20 sampled phases | §4: no per-alignment statement is claimed. `rules.md` asks for success over uniform `Q`, uniform `z` and fresh coins; the sweep is supporting evidence, and the 0.39 quantile is at most 1.0527 at every alignment against a cap of 1.16. |

## 8.1 Answers to the rejection of `a7ac78ea` (v2)

The review accepted the reduction, the collision-recovery equation and the §D
model block ("the completed initial run reproduces the reported toy successes,
schedule comparisons and model calculations"). It refused three items of support.
Each is repaired in the code, not in prose.

| objection | fix |
| --- | --- |
| the toy drew `A_h, B_h = rng % q`, so it ran a walk with wrapping gap sums while production never wraps — "quadrature, Monte Carlo and alignment sweeps all evaluate the assumed model; they do not independently validate that applicability" | §2/§5: the toy now draws `A_h, B_h` from `[1, 2^cw)` with `cw = bitlen(q) - bitlen(cap) - 1`, asserts `cap*2^cw < q`, and at every solve asserts the *true* gap exponent `bi - tbi` is in `(0, q)` and equals `(b - tb) mod q`. Reported: `coefficient_bits = 12`, `cap_times_2_to_cw_lt_q = true`, `wrap_violations = 0`, `degenerate_relations = 0` over the 600 declared trials. The positivity lemma is now a measured property of the walk that is actually run. |
| peak storage is 13 records, not 12; `max_live` was updated before insertion | §3: the law is `n_stored = k*log2(W) + 1 = 13` (12 compared, 13 stored), `max_live` is updated *after* insertion, §B1 enumerates the stored peak (13, first at `h = 32`) beside the compared peak (12), and the ledger prices `13 x 680 = 8,840 B` for a peak total of 227,416 B. |
| a hard evaluation cap plus an *expected* filter cost is not a hard work bound | §6: the filter now carries a hard allowance of `2^70` full element comparisons and aborts on exceed, charged deterministically as `2^70 * 216 / 20000 = 2^63.46` units. The walk implements the abort; §B3 reports `fp_aborts = 0` and `max_full_compares = 1` at the declared inputs. |

The §D model block was not attacked and is unchanged: `c(0.39) = 1.0507` offline,
1.0507 in `experiment.py`; Monte-Carlo 1.0453; 16-alignment band
`[1.0463, 1.0527]`; success at the cap 0.4518 nominal, 0.4501 worst alignment.

## 9. Reproduce

    python3 experiment.py <<'EOF'
    {"production_checks":2,"table_checks":2,"trials":600,"toy_bits":26,
     "R":256,"cap_factor":1.16,"seconds":60,"seed":1,"model_draws":20000}
    EOF

Standard library only, about 12 s: §A production facts on BLS12-381, §B1 the
compared (12) and stored (13) enumeration, §B2 detection law versus direct
simulation, §B3 toy search at the declared `R` in the production coefficient
regime with the wrap counter, §B4 jump-count sensitivity, §C the exact ledger and
memory, §D the model quantiles. `results.json` is the captured stdout of exactly
that command; `tools/capture_v3.py` runs it *from* `submission.json`'s declared
inputs, so the two cannot drift.
The claim is `max(127.68, computed total)` and the file prints both, so the
submitted score can never be tighter than the computation.
