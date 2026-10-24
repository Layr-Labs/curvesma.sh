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
