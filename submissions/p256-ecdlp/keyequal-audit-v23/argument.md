# p256-ecdlp v23 -- key-equality hazard measured on the shipped predicate, moment-based success, locator access priced

v22 (`57eef28`, `bf2faa9`) was rejected twice. This tree answers the four named
flaws with an executed predicate measurement and a re-priced ledger; it does not
re-argue any v22 number that the verdict did not accept.

## 1. Row occupancy is not key equality -- the predicate is now measured directly

v22 fire 2, flaw 2, verbatim: *"The new experiments measure collisions of 20-bit
selector rows, not collisions of complete P-256 x-coordinate keys. ... The planted
same-point control verifies recovery conditional on a deliberately constructed
collision, not its occurrence probability."*

The production collision predicate is **key equality inside a bucket**:

```
bucket = canonical_key & (R-1)              # the dictionary bucket
probe  = scan the bucket, compare the 32-byte canonical x
hit    = key equal  (a shared bucket with a different key is NOT a hit)
```

At production `R = 2^20` and `r = 2^255` sign classes, so same-bucket pairs occur
`r/R = 2^235` times more often than key-equal pairs. A row count therefore cannot
estimate the key hazard, and v23 does not use one.

**What is measured.** Key equality at `r = 2^255` is unobservable at any feasible
`T`. The exact predicate is measured on a **real elliptic curve of the same
shape** whose group is small enough that complete canonical-x keys collide:

| property | production | v23 fixture |
|---|---|---|
| curve | `y^2 = x^3 - 3x + b` over F_p, P-256 b | `y^2 = x^3 - 3x + 76` over F_65519 |
| group order | prime, ~2^256 | 65447, prime (counted in selftest) |
| sign classes r = (n-1)/2 | 2^255 | 32759 |
| bucket count R | 2^20 | 1024 (ladder 256..4096) |
| walk | `T <- T + [a_row]P + [b_row]Q`, `Q=[x]P` | identical |
| selector | `row = (x ^ t*2654435761 ^ trail*7919) mod R` | identical |
| record | 32 B canonical x || identity | identical |
| probe | bucket scan, key compare | identical (shipped `probe_insert`) |

`experiments/key_equal_v23.py` runs three launch mixes at the production key load
`lambda = C(204,2)/32759 = 0.6328`, 6000 seeded trials each, and counts **key-equal
pairs**, not row arrivals:

| quantity | measured | model | note |
|---|---|---|---|
| key hazard per pair | 2.857e-5 | 1/r = 3.056e-5 | ratio 0.935 |
| bucket hazard per pair | 9.924e-4 | 1/R = 9.766e-4 | ratio 1.016 |
| conditional P(key | same bucket) | 0.028917 | R/r = 0.031293 | ratio 0.924 |
| per-trial key collisions | mean 0.5915 | lambda 0.6328 | Poisson chi2 (3 dof) 1.04 |
| factorial ratio E[X(X-1)]/E[X]^2 | 0.9908 | 1.0 (Poisson) | 0.9851, 0.945 in the other cells |

The 6.5% deficit of the measured key hazard against the uniform `1/r` model is
**measured, reported and charged**: the success block uses the measured rate, not
the model. It is not explained away.

**The conditional decomposition** the verdict asked for is the measured ratio
`P(equal | same bucket) = 0.028917` against `R/r = 0.031293` (ratio 0.924), and the
same quantity at R = 256 / 2048 / 4096 measures 0.007451 / 0.060476 / 0.116598
against the model 0.007814 / 0.062513 / 0.125048. A shared bucket with a different
key is compared and **not** counted: the planted control reports the same-bucket
different-key pair as one comparison, zero collisions, while the old row predicate
would count it.

**q is re-measured on the key predicate**, not on row arrival: 3549/3549 nontrivial
key-equal pairs recover x through the probe path (yield 1.0). 1707 are same-point
(`T = T'`, difference formula) and 1842 are **opposite-point** (`T = -T'`, sum
formula); a difference-only replay would fail on 52% of the collisions, so the
shipped replay rule reconstructs the exact point from the identity and uses the
matching linear relation. Zero trivial duplicates (distinct launches `c + d*x`,
no duplicate starts).

## 2. Cardinality is not a union bound -- measured moments, explicit inequalities

v22 fire 2, flaw 3: *"A union bound needs bounds on the probabilities of pairwise
state coincidences; state-space cardinality qR does not imply those probabilities."*

v23 states every inequality it uses and never treats a row or state count as an
event probability:

```
X = number of key-equal pairs among N = 2^127.67 records, r = 2^255
E[X]      = N(N-1)/(2r)                   = 0.632878      (arithmetic)
E[X(X-1)] = E[X]^2                        = 0.400535      (uniform-key model)
E[X^2]    = E[X(X-1)] + E[X]              = 1.033413
Paley-Zygmund: P(X>=1) >= E[X]^2/E[X^2]   = 0.387584      (distribution-free floor)
Poisson:       P(X>=1) = 1 - exp(-E[X])   = 0.468939
union bound:   P(X>=1) <= E[X]            = 0.632878      (upper only, not used)
```

The **measured** factorial moments are Poisson-tight (0.9908 / 0.9851 / 0.945 vs
1.0), and the 3-dof Poisson chi-square of the per-trial counts is 1.04 (critical
7.815). The charged success is the **minimum Clopper-Pearson 95% lower bound over
the three seeded key-collision cells**: 0.435029 (2686/6000), 0.448990 (2770/6000),
0.448325 (2766/6000) -> charge **0.435029 >= 0.39**, attempts 1. Paley-Zygmund is
reported at 0.387584 as the distribution-free floor; it sits below the measured
lower bound because it uses moments only, and the charge does not rest on it.
At the measured key rate the Poisson value is 0.446562; the CP lower bound is the
smaller and is what is charged.

## 3. Locator access is priced as its own row

v22 fire 2, flaw 4: *"That control does not implement access to the identified
page or node. The seven-instruction translation price still assumes that loading a
multiword identifier and one address add suffice for the next dereference."*

`experiments/store_v23.py::anchor_chain_counters` **executes** the 128-bit anchor
chain for the 122-bit page number (6 levels, fanout 2^21) and counts every
instruction class: 2 word loads, low add, high adc, carry select, digit shift,
digit mask, node-base + slot add, tag compare, branch = **10 Ir/level = 60 Ir per
access**. The shipped probe path performs **4 bucket accesses per record**, so the
charged row is **4 x 6 x 10 = 240 Ir/record**, added to the ledger on top of the
unchanged v22 sort price (no double charge: the v22 per-pass page walk is retained
inside the sort price). The locator row including that page-walk share is
**246.518 Ir/record**; doubling the page width (4096 -> 8192 B) moves it to
**243.259 Ir/record**, so the row responds to the layout (control (c)).

## 4. Work ledger (machine-copied from `evidence/ledger_v23.json`)

| term | Ir/visit |
|---|---|
| affine addition | 606 |
| batch-inversion share | 444 |
| batch chain share | 0.054340 |
| table read/index | 10 |
| canonical sign rule | 20 |
| record append (48 B) | 12 |
| global radix sort + witness + directories (v22 price) | 1258.880229 |
| **locator access, 128-bit x 4 probes** | **240** |
| **per visit** | **2590.934569** |

Units 2.590934569 at C = 1000; log2 rows 129.0434725825; replay 1.643e-04
relative, returns 2.745e-19 (unchanged first-moment bound), table build folded
in, attempts 1. **log2 work = 129.0437061946**, claim **129.0438 = 645219/5000**
(ceiling at 1e-4; one ulp lower would understate -- control (d)).

Movement from v22: locator access +240 Ir/record; claim 128.9035 -> 129.0438
(+0.1403 bits). The success predicate changed from row arrival to key equality
(CP 0.424555 -> 0.435029, both above the floor). Nothing else was re-priced.

## 5. Four planted controls (`evidence/audit_v23.json`)

* **(a) predicate** -- a same-bucket different-key pair is compared by the probe
  path and produces zero collisions; the same runs count 3549 key pairs against
  122732 bucket pairs, and the old row predicate would count the planted pair.
* **(b) moment** -- a synthetic overdispersed X with the carrier mean gives
  Paley-Zygmund 0.010 vs 0.372 Poisson at the same mean; the bound responds.
* **(c) locator page width** -- 4096 -> 8192 B moves the locator row by -3.259
  Ir/record (control cannot be vacuous).
* **(d) regen 1-ulp carry** -- claim bounds the exact log2 work and one ulp lower
  would fall below it; the reduced rational equals the claim exactly.

The curve self-test counts the order (65447, prime), checks `n*G = O`, checks the
base point on the curve, and replays 24 walk steps against an independent scalar
multiple with 0 mismatches.

## 6. What is claimed and what is not

Claimed: total classical work **2^129.0438** target-group-operation equivalents
(reduced rational 645219/5000) with success **>= 0.435029** charged from 6000
seeded key-collision trials on the shipped predicate, at the production key load
lambda = 0.6328.

Not claimed: that the 16-bit fixture group is the P-256 subgroup; that the key law
measured at r = 32759 is proven at r = 2^255; that a production-size attack was
executed. The step from 32759 to 2^255 sign classes is the disclosed class-scale
equidistribution conjecture; the measured key hazard, the measured conditional
law and the measured deficit are the evidence at the fixture scale, and the
production row hazard on real P-256 (v22 evidence, retained) shows the selector's
source map does not carry the row hazard.

## 7. Reproduction

```
python3 experiments/key_equal_v23.py selftest   # curve order, nG=O, walk check
python3 experiments/key_equal_v23.py all        # ~18 s, 6000+1500-trial cells
python3 experiments/store_v23.py                # locator row + page-width control
python3 experiments/ledger_v23.py               # writes claim_v23.json
python3 experiments/audit_v23.py                # four planted controls
```

`experiments/production_selector_v23.py` is the retained real-P-256 row-hazard
runner whose output is `evidence/production_selector_v23.json` (v22 measurement,
same code); `returns_v23.py` regenerates the first-moment returns bound.
