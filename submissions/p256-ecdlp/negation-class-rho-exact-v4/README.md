# negation-class-rho-exact-v4 (p256-ecdlp)

Claim: time_log2 = 3193/25 (127.72). This is negation-class parallel Pollard
rho on P-256 with a global work meter and an **exact** report store. The walk,
its analysis and its toy harness are those of submission
5fa78d62-1c6b-46b8-b9cf-02481a49e645 (negation-class-rho-v1, 129.46), declared
as the source. That package adapts the class walk of submission 1c600df4
(mpjunior92) through a188e335 and took its P-256 allowances from 4fe789ec; this
package replaces those allowances, and replaces its report store.

What changed in the pricing: every P-256 field operation is priced by the length
of an executed AArch64 instruction listing (`field_listing.md`): Montgomery
multiplication 139, modular addition 20, subtraction 19 (charged 20), and an
inversion ceiling of 2^17 against a 53,098-instruction Fermat chain. A batched
step costs 6M + 6A + I/2^15 + 128 control = 1,086 instructions and is charged
1.086 units, its count, against 3.6 units in the source; every other price is
likewise its counted instructions divided by 1,000, held in integer
milli-units, so no price carries an unaccounted round-up.

What changed in the store: a report is one 104-byte record appended to the tail
of one of 2^16 logs, at the position the report's arrival index gives (the
index mod 2^16 selects the log and the rest is the position). No address in the
store is derived from the key, so there is no slot index, no home address, no
probe sequence, no replacement rule, no load factor and no occupied-set
assumption. The occupancy of every log is fixed by the report count alone.
Duplicates are found exactly: after phase 2 each log is LSD radix-sorted on the
full 256-bit key, the 2^16 sorted key streams are merged by a heap comparing
full keys, and records with equal keys are adjacent. Two such records are
exactly two trails that reached the same distinguished class; nothing else can
be reported as a collision. The path is charged as its parts — 4,547
instructions per report (arrival counter 128, node and position 3, routed
request 4,352, append 64), 208 instructions per record for the sort, 1,536 per
record for the key-stream merge and heap, the radix counter clears, and the
≤ 16 matched-record reads — and the store contributes no failure probability
except the count of reports against its declared capacity.

What else changed: the meter's block remainders are returned to the counter by
an exact atomic update, charged, instead of being counted as spent, so no
remainder allowance appears in the ledger; the earlier package's routed
full-store read-out does not exist any more, so no one-direction read-out is
charged; and every failure row in `proof.md` section 4 is a derived bound with
its arithmetic stated (the colliding trail's non-fruitless return at
L_max²/(2N) = 2^-119.36, the same over the ≤ 8 processed matches, the
shared-walk wait at 16·e^-20 = 3.3·10^-8, the report count, the fingerprint
tests, the Markov row and the preprocessing cap).

What changed in this revision: the declared initial input now finishes inside
the verifier's 20-second budget (it is 300 trials, measured at 3.90 s) and the
2,000-trial result is labelled as the offline budget it is; the phase-1 arrival
count is corrected for doublings (at least S_ops − 2^115.58 arrivals, which costs
1.5·10^-4 of success, so the H1 floor at the arrival count is 0.394851, not
0.395); the colliding pair's abandonment at L_max is charged as the shared-walk
wait (2·e^-20 per match, 3.3·10^-8 over the ≤ 8 matches) instead of the
cycle-detection window; the matched-record read-back is stated as its counted
4,352 instructions each (4,547-and-below items are otherwise unchanged); the
store capacity is D_cap = 2^64 records, exactly the range of its 64-bit arrival
counter, and the report-count Chernoff row is re-derived for it
(< exp(−2^59.27)); and the success margin is stated as what the arithmetic
gives, 0.392611 over the 0.39 requirement.

`experiment.py` (standard library) reads a JSON config on stdin. Part A runs
the production P-256 checks, Part B the toy runs, and Part C (`listing_checks`
> 0) executes the field listing on production elements and points and
recomputes both ledgers.

**Hosted input** — the declared initial input in `submission.json`
(`initial_inputs[0]`), the run the verifier executes inside its 20-second
budget:

```
{"curve":"q26","trials":300,"R":2048,"dp_bits":5,"seed":1,"production_checks":3,"listing_checks":3000,"pipeline":true,"seconds":15}
python3 -I experiment.py < hosted.json          # hard 20 s timeout
```

Measured on the authoring machine: `real 3.90`, exit 0, the whole JSON printed.
Part A and Part C all true (M = 139, ADD = 20, SUB = 19, inversion 53,098, a
batch of 12 production points equal to the reference, `claim_holds` true, store
report path 4,547 instructions against 4,547 charged); pipeline success 0.5033
over its own 300 fixed trials, 95% Wilson interval [0.4471, 0.5595], 0 log
overflows, maximum log load 0.533 of capacity, 28,792 records. The trial count is
the stop condition and `seconds` is a guard above the run's cost. Those 300
trials are a smoke test that the store is built, sorted, merged and solved from;
their interval is wide and they are not the rate evidence.

**Offline budget** — the 2,000-trial pipeline result, not the hosted input:

```
{"curve":"q26","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":3,"listing_checks":3000,"pipeline":true,"seconds":120}
python3 -I experiment.py < offline.json         # 19.9 s here
```

That run takes 19.9 s here, which is why it cannot be the hosted input: a slower
verifier machine pushes it past 20 s. Its result: success 0.45 over 2,000 fixed
trials, 95% Wilson interval [0.4283, 0.4719], 0 log overflows, maximum log load
0.6 of capacity, 192,149 records. The pipeline builds the declared store — 16
append-only logs, a real radix sort, a real heap merge, equal-full-key detection
— and solves only from a group of equal keys; no dictionary is used for
detection.

The offline toy rows in `proof.md` section 6 are re-run with fixed trial counts
(seeds 2501–2507 for the table rows, 2508–2511 for the long-trail test,
2513 for the production checks), so each rate carries its Wilson interval:

```
{"curve":"q26","trials":6000,"R":256,"dp_bits":5,"seed":2501,"production_checks":0,"listing_checks":0,"seconds":300}
{"curve":"q26","trials":6000,"R":2048,"dp_bits":5,"seed":2502,"production_checks":0,"listing_checks":0,"seconds":300}
{"curve":"q30","trials":1800,"R":256,"dp_bits":5,"seed":2503,"production_checks":0,"listing_checks":0,"seconds":300}
{"curve":"q30","trials":1800,"R":2048,"dp_bits":5,"seed":2504,"production_checks":0,"listing_checks":0,"seconds":300}
{"curve":"q26","trials":3000,"R":256,"dp_bits":5,"seed":2505,"classes":1,"production_checks":0,"listing_checks":0,"seconds":300}   (control)
{"curve":"q26","trials":12000,"R":2048,"dp_bits":5,"seed":2506,"pipeline":true,"production_checks":0,"listing_checks":0,"seconds":300}
{"curve":"q30","trials":2400,"R":2048,"dp_bits":5,"seed":2507,"pipeline":true,"production_checks":0,"listing_checks":0,"seconds":300}
{"toy":false,"trials":1,"production_checks":6,"listing_checks":0,"seed":2513,"seconds":300}
{"curve":"q30","trials":1,"toy":false,"long_trails":150,"long_len":2048,"R":4|8,"long_escape_add":false|true,"production_checks":0,"listing_checks":0,"seed":2508|2509|2510|2511,"seconds":300}
```
