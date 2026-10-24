# p256-ecdlp v10h: the sparse pipeline measured against the dense one, the record processor built and charged, the fixture scale corrected, W = 127.772244

## What this revision answers

Verdict `b8b00688` on v10g listed three numbered defects. This revision answers each one with a
measurement or a correction, not with prose. The claim moves only by the cost of work that is now
actually priced: W = 127.772244 against 127.771170 in v10g.

| Field | v10h |
| --- | --- |
| time_log2 claim | **127.78** (= 6389/50) |
| computed W | **127.772244** (headroom 0.0078 bits) |
| counted step | **1,101** instructions (`stepcount.json:o3.longest_iteration`) |
| batch-inversion share | **0.820313** instr/step = 53,760 / 2^16 walkers |
| record-processor share | **40 x 2^-62** instr/step |
| charge per step | **1,102.460077** = 1,101 + 0.639764 + 0.820313 + 40x2^-62 |
| success, dense census (arm A) | 1,691 / 4,000 = 0.422750, **CP-95 0.409822** |
| success, sparse pipeline (arm B) | 0.2535 / 0.2445 / 0.2405 at dp_bits 4/6/8, **CP-95 0.242197 / 0.233338 / 0.229402** |
| fixture | 549,756,820,252 = **2^39.000003** classes (was mislabelled 2^40) |
| storage | 72 B/record x 2^65.632 rows = **2^71.79 B**, 0.83 bits under the 2^73 B cap |

## 1. Flaw 1 — the dense-to-sparse transfer is now measured, not asserted

v10g certified 0.39 on the **dense** every-step census and disclaimed the sparse production
predicate. The review is right that this leaves the production claim unsupported. v10h runs
**both arms on the same binary, the same fixture, the same budget**, and arm B is the declared
sparse pipeline exactly:

* DP trigger `mix64(x ^ DPSALT) & mask == 0` inside the walk;
* the Montgomery batch-inversion schedule over the 16 slots (`census_v10f41.c:297-306`);
* a **sparse** record `(x, a, b)` written on the trigger only;
* recovery by pairwise key match on sparse records and a solve `x = (a-a2)/(b-b2) mod r`,
  verified by recomputing `[x]P` and comparing the full point.

| arm | dp_bits | trials | successes | recovery | **CP-95** | records/walk |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A dense (control) | 0 | 4,000 | 1,691 | 0.422750 | **0.409822** | 812,224 |
| B sparse | 4 | 4,000 | 1,014 | 0.253500 | **0.242197** | 50,737 |
| B sparse | 6 | 4,000 | 978 | 0.244500 | **0.233338** | 12,659 |
| B sparse | 8 | 4,000 | 962 | 0.240500 | **0.229402** | 3,157 |

**The honest headline is the sparse number.** The measured CP-95 of the declared production
pipeline is **0.2294-0.2422**, i.e. **below the 0.39 floor**, and this package says so in the
claim sentence rather than quoting the dense row. What the two arms establish is the size and
sign of the transfer: sparsification at the densities the fixture can carry costs about
0.42 -> 0.24 of the recovery probability and it does not recover with density (0.2422 / 0.2333 /
0.2294 is flat-to-declining from dp 4 to 8). No monotone argument is offered for 2^-61; the
transfer remains unproven and is priced as unproven.

## 2. Flaw 2 — the record processor exists, is scheduled, and is charged

v10g's step multiplied `PROD[i]` and `SUF[i]` independently, so `inv = PROD*SUF` was not any
denominator's inverse and the "shared inversion over 16 slots" was a name, not a mechanism.
v10h ships the processor (`evidence/p256_recproc.c`) and its schedule:

* **Forward**: `pre[i+1] = pre[i] * d_i` over the n walkers of the round (n-1 muls); this is
  exactly the `PROD` row the counted step already carries.
* **One inversion**: `fe_inv(pre[n-1])`, the fixed-exponentiation ladder of the shipped header
  (`p256_fe_cc.h:160-171`): 256 squarings + popcount(p-2) = 128 multiplies = **384 fe_mul**, and
  one `fe_mul` measures **140 instructions** in the same aarch64 -O2 listing discipline as the
  step count (`recproc.json:mul_once`), giving **53,760 instructions per inversion**.
* **Back-substitution**: `inv_i = inv * pre[i-1]`, `inv *= d_i` (2(n-1) muls); each walker's
  `lambda_i = num_i * inv_i` is already one of the six multiplies the counted step charges.
* **Batch size is derived**: n = 2^16 walkers per round (`recproc.json:batch_walkers`), 128 B of
  walk state each = 2^23 B, so the marginal inversion charge per step is 53,760 / 65,536 =
  **0.820313 instructions**.
* **Record processing** per report is 40 instructions (`recproc.json:recproc`: record build +
  hash probe + store) at the declared 2^-62 density, i.e. 40 x 2^-62 per step — small, but
  charged rather than assumed away.

Charge = 1,101 + 0.639764 (the v10f rare-event reserve, unchanged and labelled as carried) +
0.820313 + 40 x 2^-62 = **1,102.460077**; W = 127.631517202749 + log2(1102.460077/1000) =
**127.772244**; claim **127.78**, headroom 0.0078 bits.

**The reserve is not an assigned success margin.** v10g carried `0.639140 x (1+2^-10)`; it is
reported here as a carried rare-event term with its provenance, and the success number in the
claim sentence is the *measured sparse CP-95*, with no percentage added.

**Point-walk test — status NOT CLAIMED.** `evidence/p256_point_walk_test.c` includes the shipped
counted step verbatim and drives it on real P-256 points, and `p256_point_walk_test.py` is the
bigint oracle for it; the driver's own table construction still disagrees with the oracle at row 3
of 2,048, so this revision does **not** claim a point-walk pass. The kernel evidence of record
remains the field-primitive test (0 / 10,770 mismatches, control 25,772).

## 3. Flaw 3 — fixture scale and the cap, in exact integers

* The fixture is `(q+1)/2 = 549,756,820,252 = 2^39.000003` classes, **not 2^40.000**. Corrected
  here and in `argument.md`; `ledger.json:fixture.classes_log2` carries the machine value.
* Sparse record layout: 32 B `x` + 32 B packed `(a,b)` + 8 B tag = **72 B**
  (`cap.json:record_bytes`). At S = 2^127.632 steps and 2^-62 density the run emits 2^65.632
  records = **2^71.79 B**, fitting the 2^73 B cap with **0.83 bits** to spare.
* Dense every-step records at the same S would need 72 x 2^127.632 = 2^133.79 B, i.e. **61.17 bits
  over the cap**. That is why production is sparse: the dense configuration is the census
  instrument, the sparse one is the production artefact, and §1 measures the difference instead of
  assuming it away.

## 4. What this package claims and does not claim

Claims: the conditional 2-cycle lemma of v10g (unchanged, with its sign condition), the kernel test
0/10,770 with a 25,772-mismatch mutation control, the dense census CP-95 0.409822, the price
W = 127.772244 at charge 1,102.460077, and the storage arithmetic of §3.

Does **not** claim: the 0.39 floor for the sparse production pipeline (measured 0.2294-0.2422,
below the bar and stated as such); any production-length transfer of that statistic; a point-walk
pass; any authenticated prior submission (`dependencies.json` is `[]`).
