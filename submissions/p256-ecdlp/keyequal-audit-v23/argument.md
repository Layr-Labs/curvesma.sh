# p256-ecdlp v23 rev 2 -- key equality on the shipped predicate, full-period trails, a resolved 128-bit locator, measured cycle charge

Fire 1 of this tree (`9c353dd`) was rejected with five numbered flaws. This
revision answers each one in the code, not in prose. Claim **128.9579 =
1289579/10000**.

## Flaw 1 -- the bucket scan is not the production probe

Fire-1 verdict, verbatim: *"probe_insert scans every existing entry in
key&(R-1). For N=2^127.67 insertions into R=2^20 buckets, total comparisons ...
approximately 2^234.34 ... The asserted four-access price does not account for
these scans."* Correct, and the fix is structural:

* **Production collision finding is the 16-pass radix sort itself.** Equal
  256-bit keys are adjacent in the sorted stream; the find pass is **one 32-byte
  adjacency comparison per record** (8 Ir), priced in `store_v23.py`
  (`collision_finding`). No hash bucket exists at production scale and no bucket
  scan is charged anywhere. The locator row is one 128-bit page resolution
  (60 Ir) plus that comparison: **68 Ir/record**.
* **The fixture bucket is a measurement device for the same predicate.** With
  204 records and R = 64 buckets it costs ~325 comparisons per trial (measured,
  ~1.6 per record) and its only role is to make key equality observable at a
  feasible scale. It is not the shipped store and it is not priced as one.

The key hazard measured through that device, with full-period trails (below), is
3.0595e-5 per pair against 1/r = 3.0560e-5 (ratio 1.0012), and
P(key | same bucket) = 0.0019638 against R/r = 0.0019558 (ratio 1.0041).

## Flaw 2 -- the anchor chain is now executed and the layout control is real

Fire-1 verdict: *"anchor_chain_counters increments counters ... without loading
any directory node, testing a stored tag or resolving a page. Its low-word
evolution never incorporates the high-word page bits ... the page-width control
retains pair_stores_page=256 when the page doubles."*

v23 rev 2 contains:

* `make_directory` / `walk_page`: a **real 3-level sparse radix directory of
  node objects**; 204 page ids are inserted and every resolution is asserted to
  return the identical page object. The measured cost is 30 Ir per access
  (3 levels x 10), counted from the code that actually walks the nodes.
* `resolve128`: the same walk on a **122-bit page id held in two 64-bit words**.
  The radix digit is extracted from the combined (hi:lo) anchor (bits 63-83 and
  84-104 come from the high word), the slot offset is added with an explicit
  add/adc carry into the high word, and the high word is asserted nonzero for
  page ids above 2^64 (3 of 6 test ids, including 2^64+5 and 2^122-1).
* The page-width control recomputes **every** page-size-dependent quantity:
  8192 B pages have records/page 170 and **512** pair transfers per page (not
  256). The locator row moves 74.518 -> 74.271 Ir/record; the init/flush
  amortisation is invariant by construction (3 transfers/record at any page
  size), which is why the move is small and comes from the halved walk share.

## Flaw 3 -- full-period trails, fresh challenges, honest distribution test

Fire-1 verdict: *"The key-collision cells use one fixed small-curve secret and
trails of 3-51 steps, all shorter than one selector period."*

The rev-2 cells each draw a **fresh challenge secret per trial** (new Q = [x]P
and a rebuilt row table) and cover **at least one full selector period**:

| cell (W x S, R) | periods | trials | wins | CP-95 | mean X | E[X(X-1)]/E[X]^2 |
|---|---|---|---|---|---|---|
| 2 x 102, R=64 | 1.59 | 4000 | 1788 | 0.431511 | 0.6335 | 1.336 |
| 1 x 204, R=64 | 3.19 | 4000 | 1757 | 0.423794 | 0.6290 | 1.380 |
| 1 x 204, R=128 | 1.59 | 4000 | 1797 | 0.433752 | 0.6198 | 1.308 |

`lambda = C(204,2)/32723 = 0.6328`. The charged success is the minimum CP-95
lower bound, **0.423794 >= 0.39**, attempts 1. The model-mean chi-square is
reported honestly at 25.4 / 34.2 / 16.6 on 3 dof (mild clustering, factorial
ratios 1.31-1.38); it is **not** used to confirm the rate. The fitted-mean
chi-square is 25.5 / 33.8 / 14.3 -- the same shape, so the misfit is a small
excess of 2+ collision trials, and the CP bound, which is computed from the
actual win indicators, carries it. Paley-Zygmund at production moments is
0.387584 and is reported as a floor, below the charge; it is not the charge.

## Flaw 4 -- the cycle charge is measured, not inferred from state space

Fire-1 verdict: *"The retained return charge still infers binom(S,2)/(qR) from
state-space size."*

The rev-2 walk **detects a (point, phase) repeat and truncates the trail**; the
cycle count is then a measured per-trail event: 13/8000, 25/4000, 6/4000 in the
three cells. The charged relative loss is the **CP-95 upper bound on the worst
cell, 0.0092125** (log2 0.013230), i.e. the fraction of trails whose remaining
records are lost. No cardinality expression is used. The fixture state space is
2^22.3 and production is q*R = 2^276, so the charge is a conservative measured
upper bound, disclosed as such.

## Ledger

| term | Ir/visit |
|---|---|
| affine addition | 606 |
| batch-inversion share | 444 |
| batch chain share | 0.054340 |
| table read/index | 10 |
| canonical sign rule | 20 |
| record append (48 B) | 12 |
| global radix sort + witness + directories (v22 price) | 1258.880229 |
| **key-equality find pass: 128-bit page resolution (60) + key compare (8)** | **68** |
| **per visit** | **2418.934569** |

Units 2.418934569; log2 rows 128.9443717457; cycle charge log2(1.0092125) =
0.013230; replay 1.643e-04 relative; table build folded in; attempts 1. **log2
work = 128.9578414981**, claim **128.9579 = 1289579/10000** (1e-4 ceiling; one
ulp lower understates -- control (d)).

## Key-hazard measurements (the predicate, at full periods)

| quantity | measured | model | ratio |
|---|---|---|---|
| key hazard per pair (R=64, S=204) | 3.0595e-5 | 1/r = 3.0560e-5 | 1.0012 |
| bucket hazard per pair | 1.5538e-2 | 1/R = 1.5625e-2 | 0.9944 |
| P(key equal | same bucket) | 0.0019638 | R/r = 0.0019558 | 1.0041 |
| key hazard, R=32 / R=16 | 2.741e-5 / 2.576e-5 | 3.056e-5 | 0.897 / 0.843 |

The R=64 cell matches the uniform key law to 0.1%; the low-R cells sit 10-16%
below the model (more bucket load, more cycle truncation), which is reported and
is why the charge uses the measured CP bound rather than the model. Yield is
1.0 on the key predicate: 2516/2516 nontrivial pairs recover x (1259 same-point
+ 1257 opposite-point), so the replay must resolve the sign class.

## Four controls (`evidence/audit_v23.json`)

(a) predicate: same-bucket different-key compared, zero collisions (2534 key
pairs vs 1,290,328 bucket pairs in the same runs); planted same-key pair fires
and recovers x. (b) moment: an overdispersed synthetic with the carrier mean
gives Paley-Zygmund 0.010 vs 0.388. (c) locator: page width doubling moves the
row by -0.247 Ir and the wide page has 512 pair transfers (asserted). (d) regen:
claim bounds the exact work and one ulp lower falls below it; the rational is
exact.

## Not claimed

The 16-bit fixture group is not the P-256 subgroup; the key law at r = 32723 is
extrapolated to r = 2^255 as the disclosed class-scale equidistribution
conjecture; trails of 102-204 steps are extrapolated to 2^107.67 production
steps; the cycle charge uses the fixture's measured rate as a conservative upper
bound. No production-size attack was executed.
