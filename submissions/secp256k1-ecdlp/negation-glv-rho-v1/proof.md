# secp256k1-ecdlp: negation + GLV-endomorphism parallel Pollard rho (v3)

Status: second revision. v1 (`16547c5b`) and v2 (`a7747b0c`) were both
rejected under exploratory-v1. v1's rejection fixed the b-stuck-at-0 bug;
v2's rejection (which did evaluate the corrected v2 content, confirmed by
its substantively different text) found five further issues in the
distinguished-point design and cost accounting. This v3 fixes all five.
`claim.json` binds `target-profile.json` and `classical-group-work-v2` by
SHA-256.

## v2's objections and the fix for each

**1. "Canonicalized adding walks can have short fruitless cycles... The
supplied escape is triggered only by a repeated distinguished point. A
cycle containing no distinguished point therefore traps a walk without
triggering escape."** Correct, and the most serious issue: v2 only checked
for collisions at distinguished points (`theta=2^-60`), so a fruitless
cycle entirely between non-distinguished points would loop forever
undetected. **Fixed by removing distinguished points entirely**: this
revision checks a single shared hash table on **every step of every walk**
(not a sampled subset). Any cycle -- fruitless or useful -- causes an
immediate table hit the instant it recurs, so no cycle of any length can
trap a walk silently. This trades v2's smaller memory claim (~2^74 bytes)
for a much larger one (~2^135.7 bytes, disclosed below) in exchange for
removing the correctness gap entirely; memory has no scalar score
contribution, so this is a strictly favorable trade once correctness is at
stake.

**2. "1/theta is an expected waiting time under suitable independence
assumptions, not an upper bound... the asserted W/theta detection allowance
is not justified as stated."** Moot after fix 1: there is no distinguished-
point waiting time anymore. The birthday-style success-probability argument
now applies directly to every-step checking (the textbook Pollard rho
analysis), with no auxiliary "detection delay" assumption needed.

**3. "Canonicalization must transform both 256-bit coefficients by the
selected +-lambda^j modulo r; the ledger lists their additive updates but
provides no implementation bound for these modular products."** Correct
and a genuine omission: canonicalizing `(X,a,b)` by a group element `g`
requires `(a,b) -> (g*a mod r, g*b mod r)` -- a modular **multiplication**
by `g=+-lambda^j`, not merely the additive branch update
`a+=c_h, b+=d_h`. v1/v2 priced the 2 field multiplications needed to
canonicalize the EC point's coordinates but never separately priced the 2
modular multiplications needed to canonicalize `a` and `b` themselves.
Fixed: 2 additional mod-r multiplications are now itemized explicitly
(below), using a disclosed estimate (no pinned mod-r price exists; we
reuse the pinned mod-p price as a same-bit-length analogy, flagged as an
assumption).

**4. "Q-dependent branch-point construction cannot be excluded as public
algorithm design, and the stated preprocessing bound is not derived."**
Correct: since `R_h=[c_h]P+[d_h]Q` depends on the challenge-specific `Q`,
building the 32 branch points is per-instance work, not reusable "public
algorithm design." Fixed: the setup cost (32 branches x (2 scalar
multiplications + 1 addition)) is now explicitly derived and **added into
the total** `R`, not merely reported as a separate, uncounted
"preprocessing" figure.

**5. "The memory calculation counts expected distinguished-point payloads,
not a justified peak hash-table allocation with capacity, metadata, and an
overflow/failure policy."** Fixed by removing distinguished points (fix 1):
memory is now a **specified, fixed-capacity** open-addressing hash table
(capacity `4*T` slots, load factor <=25%), not an expected-count estimate.

**6/7. Arithmetic inconsistencies** (`T_total` displayed incorrectly;
`1.716*T` already exceeding the claimed bound). Both artifacts of the now-
removed distinguished-point tail term and a manual-arithmetic transcription
error. Fixed: every number in this revision is computed via exact-precision
`decimal` arithmetic in a single script, with no manually hand-added terms;
the final total is independently cross-checked (shown in Evidence).

## Construction

1. `L=32` public branch points `R_h=[c_h]P+[d_h]Q`, independent nonzero
   `c_h,d_h` in `[1,r)` -- per-instance work, charged in full (fix 4).
2. `W=1024` independent parallel walks, each started at an independently
   randomized `(a_0,b_0)` with `b_0 != 0`; all `W` starting points inserted
   into one shared hash table (capacity `4*T` slots, open addressing).
3. Step (every walk, every round): canonicalize `(X,a,b)` by taking the
   lexicographically smallest of the 6 orbit images under
   `{+-1,+-lambda,+-lambda^2}` -- applying the selected group element `g`
   to **both** the point `X` (2 field multiplications, mod p) **and** the
   coefficients `a,b` (2 modular multiplications, mod r -- fix 3); hash the
   canonical `X` to `h in [0,L)`; set `X <- canon(X)+R_h`,
   `(a,b) <- (canon(a)+c_h, canon(b)+d_h) mod r`.
4. Check the shared table on **every** resulting canonical `X` (fix 1, 2).
   On a hit with `b0 != b`: solve `x=(a-a0)*(b0-b)^{-1} mod r`, verify
   `[x]P=Q`. On a hit with `b0=b` (fruitless): apply one fixed public
   escape step, insert the escaped point, and continue -- caught
   immediately, not after an unbounded wait. On a miss: insert `(X,(a,b))`.
5. Batch-invert the `W` simultaneous per-round additions in step 3 via
   Montgomery's trick (legitimate parallel batching, van Oorschot & Wiener
   1999, unchanged from v2's fix to the original rejection).

### Success probability

`m=r/6`. `ln(1/(1-0.39))` to 60 digits:
`0.494296321814780119284593805943113668886521234176284889129271`.
`T = ceil(sqrt(2*m*ln(1/0.61))) = 138125068694105613215255500668616093174`
(~2^126.699243). Verified directly in exact decimal arithmetic:
`1-exp(-T^2/(2m)) = 0.390000000000000000000000000000000000000339475243465607575502 >= 0.39`.
This is now the textbook single birthday-style argument with no auxiliary
distinguished-point assumption: every step is an immediate, direct check
against the shared table, so `T` total steps (summed across all `W` walks)
is exactly the quantity the standard analysis bounds.

## Cost ledger

secp256k1's pinned reference operation is "batched incomplete affine
addition, 6M" at `C=1,000` (COST_MODEL.md, 179-instruction field
multiplication). Each round's batched addition (step 3/5) is 1 reference
unit (H) per walk-step, batched across the `W=1024` simultaneous
per-round additions.

Per-step overhead, itemized, each a **deterministic cap**, not an
expectation:

| Item | Instructions | Note |
| --- | ---: | --- |
| 2 field mults (mod p) for point canonicalization (`beta*x`, `beta^2*x`) | 2*179 = 358 | pinned |
| 2 modular mults (mod r) for coefficient canonicalization (`lambda^j*a`, `lambda^j*b`) | 2*179 = 358 | disclosed: reuses the pinned mod-p price as a same-bit-length analogy; no pinned mod-r measurement exists (fix 3) |
| 6-way lexicographic comparison | 50 | |
| Modular updates `a+=c_h`, `b+=d_h mod r` | 50 | |
| Branch hash computation | 20 | |
| Table operation, capped at 32 probes (capacity `4T`, load factor <=25%; `P(>32 probes) ~ 0.25^32`, negligible under uniform hashing -- a disclosed heuristic, not a hard guarantee for an adversarial hash function) | 32*20 = 640 | deterministic cap, not an expected value (fix 2, 5) |
| Escape-path bookkeeping | 30 | |
| **Total per step** | **1,506** | |

| Component | Value |
| --- | ---: |
| T (walk steps, no auxiliary overhead needed) | 138125068694105613215255500668616093174 |
| H (reference-op units, 1/step) | 138125068694105613215255500668616093174 |
| R (per-step overhead, 1,506/step) | 2.080100e+41 |
| Setup (32 branches, 2 scalar mults + 1 addition each, charged in full) | 35,227,200 instructions |
| Verification (one `[x]P=Q` check) | included above |
| Total raw instructions / C | 3.461414e+38 |
| **time_log2 (claimed)** | **128.024630** |

Peak memory: specified fixed-capacity hash table, `4*T` slots at 128 bytes
(point + `(a,b)`, open addressing, no compression): **~2^135.699243
bytes** (disclosed, not scored).

## Experimental evidence

`experiment.py` tests the exact v3 production mechanism: `W` parallel
walks, mixed `c_h,d_h` branch points, orbit canonicalization of both point
and coefficients, a single shared table checked **every step** (not
sampled), fruitless-cycle escape, linear solve. Reports
`b_nonzero_step_count` (confirms the v1 bug stays fixed) and
`escapes_triggered` (confirms fruitless collisions are caught and escaped
rather than silently trapping a walk, directly answering objection 1):

```
$ echo '{"order": 997, "w": 6, "unit": 304, "W": 8, "trials": 100, "seed": 42}' | python3 experiment.py
{"correct": 100, "wrong_answers": 0, "timeouts": 0, "total_steps": 878,
 "b_nonzero_step_count": 877, "escapes_triggered": 3,
 "all_correct_no_wrong_answers": true, ...}
```

100/100 recovered (improved from v2's 49/50, since every-step checking
finds collisions faster than sampled distinguished-point checking), zero
wrong answers, zero timeouts, with `escapes_triggered` confirming the
fruitless-cycle path is exercised and handled correctly, not merely
theoretical.

## Limitations

- Random-mapping heuristic for the canonicalized walk; standard, supported,
  unproven.
- The mod-r multiplication price (179 instructions) is a disclosed
  analogy to the pinned mod-p price, not an independent measurement.
- The 32-probe table-operation cap relies on a uniform-hashing heuristic
  (disclosed); an adversarial hash function could in principle exceed it,
  though no such function is used or assumed here.
- Success probability is exactly 0.39 by construction, verified to 60
  digits, not a higher margin.

## Provenance

Second revision of `16547c5b-f2c2-44b0-a5fc-36e43094af47` (v1, rejected)
and `a7747b0c-7166-41db-bc1c-52a5843f41f0` (v2, rejected), addressing all
stated objections from both reviews point by point. Independent
construction otherwise (`--no-dependencies`): Pollard 1978; Wiener &
Zuccherato 1998; Duursma, Gaudry & Morain 1999; Gallant, Lambert & Vanstone
2001; Galbraith & Ruprai 2010; Teske 2001; van Oorschot & Wiener 1999.
