# negation-class-rho-exact-v5 (p256-ecdlp)

Claim: time_log2 = 3193/25 (127.72). This is negation-class parallel Pollard
rho on P-256 with a global work meter and an **exact** report store. The walk,
its analysis and its toy harness are those of submission
5fa78d62-1c6b-46b8-b9cf-02481a49e645 (negation-class-rho-v1, 129.46), declared
as the source. That package adapts the class walk of submission 1c600df4
(mpjunior92) through a188e335 and took its P-256 allowances from 4fe789ec; this
package replaces those allowances, and replaces its report store.

Pricing: every P-256 field operation is priced by the length of an executed
AArch64 instruction listing (`field_listing.md`): Montgomery multiplication 139,
modular addition 20, subtraction 19 (charged 20), and an inversion ceiling of
2^17 against a 53,098-instruction Fermat chain. A batched step costs
6M + 6A + I/2^15 + 128 control = 1,086 instructions and is charged 1.086 units;
the control budget of 128 itemises 104 instructions, including the private
meter update, the two-word trail length and the phase-1 operation down-counter.

Store: a report is one 104-byte record appended to one of 2^16 logs at the
position its arrival index gives. Duplicates are found exactly: each log is
radix-sorted on the full 256-bit key and the sorted logs are merged by a heap
comparing full keys. The path is charged as its parts: 4,547 instructions per
report, 208 per record for the sort, 1,536 per record for the merge, the radix
counter clears and the ≤ 16 matched-record reads.

What changed in this revision (`proof.md` sections 3 to 5):

1. **Every counter fits its word.** The work meter is one shared 64-bit
   counter of super-blocks of 2^70 units, cap n_cap = ⌊2^57.72⌋ − 1, with
   private 64-bit counters inside a block (a chunk index and a milli-unit
   counter over 2^50-unit chunks). An over-cap fetch-and-add halts the run;
   the counter is monotone, so each processor makes at most one. Remainders
   stay reserved, and the stranded capacity (at most 2^90 units) comes out of
   the Markov headroom. Each chunk carry and each phase-1 sub-batch open
   charges its own 8 instructions inside the block. W ≤ 2^127.72 is checked as
   an exact integer inequality over one list of four terms: the blocks, the
   fetch-and-adds on the block counter, a halted processor's straddling charge, and
   preprocessing. Phase 1 ends on a shared 64-bit counter of super-batches of
   2^70 operations, so it performs at least S_ops and fewer than S_ops + 2^79
   operations. Preprocessing has its own meter of the same structure. The
   store's arrival counter halts at D_cap = 2^64 − 2^21, so the largest index
   it ever returns is below 2^64 − 2^20. The counter table in section 3 lists
   every counter with its width, its largest value and its charge.
2. **Abandonment is length-biased.** Both trails of a match are selected
   through the collision, so the fresh-trail tail does not apply. Under H1's
   geometric model each is capped with probability at most
   nθ(1 − θ)^(n−1)/(1 − (1 − θ)^n) with θ = 2^-64 and n ≥ L_max − L_max/4
   (H1's repeat clause: a trail that reaches L_max makes at most L_max/4
   repeated arrivals): < 4.589·10^-6 per trail, < 9.178·10^-6 per match,
   < 7.342·10^-5 over the ≤ 8 processed matches.
3. **Every quoted number equals the ledger, rounded in its safe direction.**
   `check_numbers.py` prints each number quoted in `proof.md`, this file,
   `submission.json` and the submission note next to its ledger value and fails
   on any mismatch.

Result: H1 floor 0.394850 at the phase-1 arrival count, Markov early-halt row
0.002239, other failure rows < 7.35·10^-5, success ≥ 0.392537 against the 0.39
requirement (margin 0.002537).

`claim_holds` in the ledger output is an internal-consistency check of the
ledger arithmetic under the stated assumptions and prices: it confirms that
the stated terms add up, that W ≤ 2^127.72 and every counter bound hold as
exact integer inequalities, and that the success sum exceeds 0.39. It does not
validate H1, the analytic itemisations or the prices themselves.

`experiment.py` (standard library) reads a JSON config on stdin. Part A runs
the production P-256 checks, Part B the toy runs, and Part C (`listing_checks`
> 0) executes the field listing on production elements and points and
recomputes both ledgers.

**Hosted input**, the declared initial input in `submission.json`
(`initial_inputs[0]`), run inside the verifier's 20-second budget:

```
{"curve":"q26","trials":300,"R":2048,"dp_bits":5,"seed":1,"production_checks":3,"listing_checks":3000,"pipeline":true,"seconds":15}
python3 -I experiment.py < hosted.json          # hard 20 s timeout
```

Part A and Part C all true (M = 139, ADD = 20, SUB = 19, inversion 53,098, a
batch of 12 production points equal to the reference, store report path 4,547
instructions); pipeline over its own 300 fixed trials. Those 300 trials are a
smoke test that the store is built, sorted, merged and solved from; they are
not the rate evidence.

**Number check:**

```
python3 -I check_numbers.py proof.md README.md submission.json ../../../NOTE.md
```

**Offline budget**, the 2,000-trial pipeline result, not the hosted input:

```
{"curve":"q26","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":3,"listing_checks":3000,"pipeline":true,"seconds":120}
python3 -I experiment.py < offline.json
```

The offline toy rows in `proof.md` section 6 are re-run with fixed trial counts
(seeds 2501–2507 for the table rows, 2508–2511 for the long-trail test, 2513 for
the production checks), so each rate carries its Wilson interval:

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
