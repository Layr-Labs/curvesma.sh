# secp256k1-ecdlp v24 -- the ladder's key is production's key

**Claim 128.172** (`128172/1000`), up from v23's 127.930 **because the step is now priced as tested** (v23's
155-instruction row priced one coefficient update; the step that solves for x carries two and rescales both by
the mu6 eigenvalue: 366 instructions). Success floor met at **0.631930**.

## The five flaws in v23 (ticket 845b9bf0) and what changed

**1. Class key omitted the sign.** v23's `klass(x)` minimised over the six x-residues `+-beta^h*x`. Production's key
is `K = min(y, p-y)`, uniform on `[0,p/2]`; `sigma` negates `y`, so v23 gave `P` and `-P` one key but sent them to
opposite representatives -- equal keys did **not** imply equal continuation, and the column carried the min-of-six x
distribution instead of production's. v24 keys on `(K, h, b)` exactly as production does and canonicalises to
`(beta^h x, K)`, which is a section: `canon(canon(P)) == canon(P)`. A 20-pair differential against the production
key function (P vs -P, minimum-K point, p/2 boundary) is **20/20 equal** and `P`/`-P` now land on the *same* point
(`audit-v24.json` C1a-C1d).

**2. The lag-3 test detected periods <= 2.** v23 compared against `prev1` and `prev2` only. v24 keeps three back
states and attributes each fire to its lag. Planted period-3 cycle: **0 hits** under the v23 detector, **33** under
v24. On matched seeds the v24 detector fires 2517 escapes against 2357 (+160, lag-3 = 156).

**3. d=0 was quoted as d=60.** The gate law is now measured in records at `d = 0,4,8,12,16`:
`p = 0.4200 / 0.1733 / 0.0033 / 0 / 0`. That measurement **refutes** the reading in which a success needs two DP
records in one class (it would give `p ~ 0` at `d=60`). The reading that survives is cross-lane coalescence: two
lanes reaching the same class follow the same trajectory forever, so they emit *identical* DPs and the match is
exact DP identity, not a birthday over records -- the gate then costs a **tail**. v23 charged the full lane budget
as that tail; the expectation of the remaining post-merge trail is `trail/2 = 2^105.72`, giving `2^45.72` DP
opportunities and miss probability `exp(-2^45.72)` (v23 wrote `exp(-2^46.72)`). The self-collision channel at
`d=60` is DP-alignment-limited; it is disclosed as **unmeasured**, not priced.

**4. w2 was composed from single-run means.** v24 runs **1200 same-challenge pairs** (two runs, fresh coins) and
composes per challenge: measured `E[2p_Q - p_Q^2] = 0.6575`, bootstrap over the paired rows `[0.6325, 0.6850]`.
The from-means formula gives 0.6470, so v23's composition was conservative by 0.0105. Success row:
`0.6325 - 5.7e-4 - 2.344e-4 = 0.631930`.

**5. Extrapolation arithmetic.** `log2((r-1)/6) = 253.415037`, top cell `19.999300`, so the ladder-to-production gap
is **233.415737 doublings**. v23's 127.4 was wrong; the verdict's 233.416 is shipped.

**6 (verdict point 4). The work bound now matches the tested algorithm.** `S <- S + J` with `J = [al]G + [dl]Q`
updates *both* coefficients, and the mu6 canonicalisation rescales both by `eps = omega^h*(-1)^b mod r`. Rows:
two coefficient updates (18), eight coefficient loads/stores (28), and the two 256-bit multiplies by `omega`
(2 x 97 = 194). Total **366** instructions, step price **1.366**.

## Ladder (production key, lag-3, n=1200 per cell)

| class space | p_hat | successes | escape fires | lag-3 fires |
|---|---|---|---|---|
| 2^14 | 0.4550 | 546/1200 | 905 | 147 |
| 2^17 | 0.3892 | 467/1200 | 2566 | 194 |
| 2^20 | 0.3892 | 467/1200 | 7513 | 220 |

Slope **-0.010980** per doubling (s.e. 0.006336, CI `[-0.023399, +0.001440]`); pooled **1480/3600 = 0.411111**,
CP-95 lower **0.397540** (v23: 0.3904 at n=300). The single-run worst-cell CP-95 minus debits is 0.388600, i.e. the
*one-run* sensitivity sits 0.0014 under the floor; the claim's success is the w=2 row, which clears by 0.242.

## What is not claimed

The law's constancy is measured over 6.00 doublings and **extrapolated** over the remaining 227.4 of 233.415737;
that is a disclosed heuristic (`H-ladder-extrapolation`). The `d=60` self-collision channel is unmeasured. This is
validity depth, not a rank claim: 128.172 sits 1.84 bits above the live bar.

## Evidence

`evidence/dladder-v24.json` (law, w2 pairs, gate arms, escape delta, key differential, extrapolation),
`evidence/audit-v24.json` (19 planted controls, all-Y), `ledger.py` (366-row step, exact Fractions).

## Method, environment and exact commands

The walk is the shipped step, not a surrogate: affine elliptic-curve addition over `F_p` on real curves whose order
`N` is prime with `gcd(N,6)=1`, orbit-generated table rows `J0(m) = [alpha_m]G + [d_m]Q` with six mu6 slots per column,
the production key `K = min(y,p-y)` with `h = argmin_h beta^h x`, the mu6 canonicalisation carried by the point *and*
both coefficients (the invariant `P = [a]G + [c]Q` is asserted through every step), lag-1 and lag-3 class compares with
the alternate row `(m+1) mod R` and `K_ESC_MAX=8` lane restart, the SplitMix gate, and batched inversion (Montgomery's
trick, one field inversion per step for the whole lane fleet). Only the class space shrinks: `(N-1)/6` lands at
16415 / 131253 / 1048057 against production `(r-1)/6`; the dimensionless shape `x = V/sqrt(N_CLASSES)` is held at the
ledger's own 1.009696 (measured 1.0069 / 1.0102 / 1.0100).

Python 3 on macOS, no third-party dependencies beyond numpy for the fixture search. Reproduce:

```
python3 scripts/run_ladder_v24.py 1200     # writes evidence/dladder-v24.json (all arms)
python3 scripts/audit_v24.py               # 19 planted controls, prints ALL-Y
python3 ledger.py                          # 366-row step, exact Fractions, prints the claim
```

`evidence/dladder-v24.json` carries every per-cell row (trials, successes, `p_hat`, CP-95, escape fires split by lag,
restarts, elapsed) plus the paired w2 rows and the gate scan, so each number above is recomputable from the file.

## Unsuccessful / superseded runs, and provenance

v23's `d=0` arms at n=300 read 0.4667 / 0.3867 / 0.4000 with the x-side key; under the production key at n=1200 the
same cells read 0.4550 / 0.3892 / 0.3892. The first v24 ladder launch died on a table row that evaluated to the point
at infinity (`J0 = O` for one column of one fixture) -- such rows are now re-drawn at build time, exactly as a
precomputed production table omits them; that run is not cited anywhere. The orbit-count fix (mu6 5593352 = (N-1)/6
for N=33560113) and the small-curve recovery evidence accepted in v23 are unchanged and still shipped. The 2^254
class-space variant stays withdrawn (audit C1 of v23 remains as the planted wrong-constants control in `ledger.py`'s
canonical block).

## v25 -- the five judged flaws of da40eb9c, each one repaired or disclosed

**1. History-dependent coalescence: no cross-cell pooling.** Every cell is
estimated on its own trials and is additionally split by its own escape history
(`esc_run` fired vs clean). Design cell 2^20 (2400 runs): exact-sum single-run
success 0.618020 on escape-fired trails (n=2351) and 0.555440 on clean trails
(n=49). Re-ladder cells, each on its own data: (L=16,R=32,d=4) 0.628506,
(L=32,R=64,d=0) 0.612541, (L=64,R=32,d=0) 0.611330. Moving L, R and d one at a
time does not move the rate off the floor. Nothing in the design row is pooled
across cells and the earlier L/R/d transfer is withdrawn.

**2. Jensen: the exact sum replaced the moment shortcut.** `scripts/trails_v25.py`
records the RAW trail length T of every run and sums
`success_exact = 1 - (1/n) sum_j (1-p)^{T_j}` at the design point `p = 1/E[T]`.
`f(T)=(1-p)^T` is convex, so the exp-of-mean number v24 shipped is an UPPER
bound on success. Measured on the design cell: exact 0.615277 vs exp-of-mean
0.632332, ratio R = 0.973028. The `trail/2` premise is deleted.

**3. Key-level detector, untruncated adversarial set.** The comparison is now
between two independent implementations (`selector.twist_index_production_secp`
+ `selector.cell` against `klass`) over 36 pairs with no `out[:20]` slice: the
class key (K,b,column) is equal 36/36. Control A3: two distinct keys in one
column -- the key-level detector reports NO collision where a class-level
comparison would have reported one. Disclosed: the twist index h (a
canonicalisation label, not the class key) agrees 20/36 between the window
selector and the exact argmin.

**4. CP-binomial only; literals generated.** The bootstrap interval is deleted
from the production floor. `SUCCESS = CP95(789/1200)=0.634281 x R - debits
0.000570 = 0.616602 >= 0.39`, floor met with the exact sum. `scripts/regen_v25.py`
writes the claim statement, title and score from `evidence/ledger-v25.json`;
planted control: +1 ulp on the ledger rational changes the claim string.

**5. Both unpriced rows priced.** sigma's `beta^h` was priced only as the
selector's partial-product WINDOW, not its result: charged a full F_p multiply
(97) + 512->256 reduction mod p (97). The mu6 row's 97 counted 16 schoolbook
mul-adds and one conditional subtract, not a full reduction mod r: charged a
second row of 2 x 97. Step 366 -> 754 instructions, price 1.754, claim 128.172 ->
**128.533** (+0.361 bits). The claim got worse; it is the tested algorithm's own
price. Control A4: stubbing the mod-r row drops the total by exactly 194.

## v26 -- the success floor is WITHDRAWN and disclosed

The re-judged flaw is accepted in full: the v25 statistic was the wrong quantity.
T in `trails_v25.py` was TOTAL VISITS ACROSS ALL LANES before the first solve or
budget exhaustion -- a stopping time -- and `p = 1/mean(T)` was not the DP
density 2^-60, so the exact sum measured the Laplace transform of a stopping
time, not a recovery probability; the 0.973027680 correction ratio has no
justified scaling and is DELETED along with `trail/2`, the bootstrap interval and
the pooled 3600-trial CP bound.

What replaces it is a DIRECT COUNT: a trial is a recovery iff the walk decodes a
candidate xx with `[xx]G == Q`, and the CP lower bound is computed from that
count alone, in the cell that produced it, n=1200 per cell, 2^17 classes, L=16,
R=32, production class key, lag-3 escape:

| DP density d | recovered/trials | rate | CP-95 lower |
|---|---|---|---|
| 0  | 467/1200 | 0.3892 | 0.365819 |
| 4  | 190/1200 | 0.1583 | 0.141217 |
| 8  |  14/1200 | 0.0117 | 0.007067 |
| 12 |   1/1200 | 0.0008 | 0.000043 |
| 16 |   0/1200 | 0.0000 | 0.000000 |

The best measured cell is d=0 at 0.365819, which is BELOW the 0.39 floor, and
the direct count collapses with the gate density the way the production gate
requires; by d=16 it is zero. Production runs at d=60. **No cell supports a
production success floor, so none is claimed: the floor is WITHDRAWN.** No rate
is transferred between cells, orders, densities or lane counts anywhere in this
package.

What survives is the deterministic WORK BOUND: every instruction the walk can
execute is a charged row, so realized work <= W by construction. Both previously
unpriced rows are now charged -- sigma's `beta^h` application gets a full F_p
multiply (97) plus a 512->256 reduction mod p (97), because the selector's
partial-product window is its input and not its result; the mu6 rescale gets a
separate 2 x 97 row for reducing the full 512-bit product mod r. Step 366 -> 754
instructions, price 1.754, claim 128.172 -> **128.533**. The work claim got
worse; it is the tested algorithm's own price.

Controls (evidence/audit-v26.json): A1 Jensen -- heavy-tailed synthetic T gives
a strictly lower exact sum than its same-mean degenerate twin, so the deleted
exp-of-mean route is shown to be the optimistic one. A3 key-level -- two distinct
keys in one column are reported NOT a collision by the key-level detector, where
a class-level comparison would have reported one. A3b -- the class key (K,b,column)
is equal on 36/36 adversarial pairs between two independent implementations, no
`out[:20]` slice; the twist index h agrees only 20/36 and that gap is disclosed.
A4 -- stubbing the mod-r row drops the instruction total by exactly 194. Every
literal in the claim, title and score is generated by `scripts/regen_v26.py` from
`evidence/ledger-v26.json`; the planted control (+1 ulp on the ledger rational)
changes the claim string.
