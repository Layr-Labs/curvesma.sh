# Signed BSGS over a distributed merge network

Candidate folder for `tracks/secp256k1-ecdlp`.

- `submission.json` — the single scored claim, its rational score `1305909/10000 = 130.5909`, and the declared experiment.
- `claim.json` — a mirror of that claim with its quantities in plain fields; it adds no assertion.
- `proof.md` — the complete argument: target, construction, exact coverage proof, the itemised machine model and ledger, the memory bound, and the no-heuristic paragraph.
- `accounting.json` — every ledger price, count and integer term, plus the exact bound check.
- `experiment.py` — standard library only, bounded evidence at the scope described below.
- `README.md` — this file.

## The claim in one paragraph

For the protected secp256k1 target, `m = h = 2^127`, `d = 2m-1`, `D = -[d]P`. Two lists are built: `(x([j]P), j)` and `(x(Q+[i]D), i)`. Keys are x-coordinates with an infinity tag. Each list is stored in a complete 64-ary lazy merge tree of depth 22 with one record per leaf, and a single charged engine intersects the two root streams. An equal key gives `x = i*d ± j (mod r)`; the candidate scalar is verified by `[a]P = Q`. The covered scalar set is exactly the union of the `h` adjacent intervals `[id-(m-1), id+(m-1)]`, of `COV = h*(2m-1) = 57896044618658097711785492504343953926464851149359812787997104700240680714240` residues, and `2*COV - r = 92138019644721193389046258963199934143 > 0`, so success is at least `1/2`. There is no heuristic and no expected value used as a bound.

Total instructions `T = 2050001412324430834062893202469342222401024`, so `log2(T/1000) = 130.5908206043` and the submitted bound is `1305909/10000`, checked by the exact inequality `T^10000 < 1000^10000 * 2^1305909` (fails at `1305908/10000`). Peak storage is `115804030901347946931027802910112881136640` bytes, below `2^137`.

## Running the experiment

Standard library only, one JSON request on stdin, one JSON result on stdout:

```sh
echo '{"mode":"all","m":16,"h":16,"block_size":4,"fanout":4}' | python3 -I experiment.py
```

Initial input: `{"mode": "all", "m": 16, "h": 16, "block_size": 4, "fanout": 4}`.

Modes: `all`, `toy`, `merge`, `accounting`. Integer parameters have declared minimum and maximum in `submission.json`.

### Output at the initial input

```text
batch_inversion: {"block": 4, "matched": 4}
coverage: {"covered": "57896044618658097711785492504343953926464851149359812787997104700240680714240",
           "covered_below_order": true, "success_at_least_half": true,
           "twice_minus_order": "92138019644721193389046258963199934143"}
ledger:   {"bound": "1305909/10000", "bits": 141,
           "peak_bytes": "115804030901347946931027802910112881136640",
           "total_instructions": "2050001412324430834062893202469342222401024", plus every term}
merge:    {"comparisons": 6, "duplicate_keys": 2, "emitted": 8, "ends": 4, "input": 8}
toy:      {"b": 2, "covered": 496, "covered_fraction": 0.630241, "fanout": 4, "h": 16,
           "m": 16, "order": 787, "p": 757, "sample": 26, "sample_hits": 22}
ok: true
```

### What the experiment checks

1. **Toy group, end to end.** `find_toy_curve` locates a curve `y^2 = x^3 + b` over a small prime field whose group order is prime and lies in `(3*h*d/2, 2*h*d]`, so the toy's intervals do not wrap and its coverage is between `1/2` and `2/3`, mirroring the production regime. Baby and giant records are built from curve points only; both are merged through the F-ary network; the intersection recovers the scalar and it is verified on the curve. For every scalar in `[0, order)` an independent interval computation predicts coverage, and the exact covered set is checked against the predicted contiguous union and against the count `h*d`. A sample of 12–27 boundary and interior scalars is then run end to end, and recovery must agree with the prediction in both directions.
2. **Merge network.** Sortedness, multiset preservation, duplicate-key retention and left-first tie order are checked on a record list with three duplicate keys; a node with an END child must not starve its siblings.
3. **Batch inversion.** Montgomery prefix/suffix inversion over `block_size` field elements is compared against direct inversion element by element.
4. **Production coverage.** `COV = h*(2m-1)` on the protected integers, `COV < r`, and `2*COV - r` equal to the exact integer above.
5. **Exact ledger and score.** Every term is recomputed from the itemised prices, the total must equal the declared `T`, the exact inequality `T^10000 < 1000^10000*2^1305909` must hold and `1305908/10000` must fail, and the peak-byte formula must reproduce the declared bound.

### Scope

This is bounded evidence at the executed scope. It is not a production run, not a complexity proof, not a timing measurement, and none of its counters is an AArch64 instruction count. The production claim rests on the analytic argument in `proof.md`. Candidate files are untrusted evidence in the shared verifier.
