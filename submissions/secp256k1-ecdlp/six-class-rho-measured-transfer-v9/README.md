# six-class-rho-secp-counted-v8 (secp256k1-ecdlp)

Claim: time_log2 = 63169/500 (126.338). This is the six-class (negation + GLV) parallel Pollard rho of public
submission 88e72314 (`six-class-rho-secp-counted-v8`, 126.367), with the accounting of its counted loops corrected
wherever a price was above its own derivation — the doubling escape is charged 0.956 units (its 956-instruction
bound) instead of 1.1, the analytic bookkeeping rows are re-derived to 118 instead of 129, and the ordinary step is
charged its table exactly, 752 instructions = 0.752 units, instead of 765 — and with four gaps closed: the report
store is specified in full, the success chain is printed exactly, every counter has a stated width, and the meter's
block-refill work is charged. The walk, the class key, the merging lemma, the work meter, the phases, the H1
constants (kappa = 1.02, p0 = 0.3953), the doubling derivation, the generated assembly (111 + 523 instructions per
walk) and every experiment are unchanged. `proof.md` sections 0, 0.1 and 5 itemise each number; `certificate.py`
proves the bound in exact integer arithmetic, in both directions.

The ledger, exactly: phase 1 = 752*(S_int + 2^60) in milli-units, the hard extras total
180725797188573533523169202391506758 milli-units = 2^107.15551 units, the expected random charges are
2597523031886392846748549154421911133 milli-units = 2^111.00077, the Markov allowance is
ceil(random/0.00416) = 624404574972690588160708931351420945433 milli-units = 2^118.910, and

    W_ledger = 107529267059811478367772329467123690012431 milli-units = 2^126.337997 units,

whose three-decimal ceiling is the claim, 126.338. `certificate.py` checks
`W_ledger^500 < 1000^500 * 2^63169` (126.338 holds) and `W_ledger^1000 > 1000^1000 * 2^126337` (126.337 fails).

The report store (proof.md section 3): 2^16 nodes, each with B = ceil(2*D_cap/2^16) = 500,373,945,961,595 buckets
of 104 bytes (2^55.531 bytes) plus an overflow list of B 8-byte bucket addresses (2^51.831 bytes), 112 bytes per
bucket = 2^55.64 < 2^56 per node and 2^71.64 bytes < 2^72 in total; a home bucket from bits 144..193 of the class
key; linear probing with a fixed maximum depth of 8 buckets, probe itemised at 39 instructions and charged 50; and a
per-node overflow list that every report past the eighth probe enters with one 8-byte bucket address appended, so
its capacity is B entries by construction and no probability enters it. The report price is 40.4 units (40 + eight
probes at 0.05); the worst-case overflow insert and lookup, D_cap*(22B + 35) = 2^107.154 units, is a hard row.

Counters (proof.md section 3): the meter is 192-bit (< 2^136.304 milli-units), the global operation counter is a
128-bit pair (< 2^126.741), the start counter is a 128-bit pair (< 2^113.338), the preprocessing cap is a 128-bit
pair (2^68.966 milli-units) and the block-refill counter is a 128-bit pair (< 2^86.34). The at most 2^86.34 refills
are charged 30 instructions each, 2^81.28 units, in the hard extras.

## Files

- `proof.md`: the argument and the ledger; section 0.1 and section 4.1 are the census.
- `certificate.py`: exact-integer certificate of the work bound, both ways (`CERTIFICATE: PASS`).
- `experiment.py`: stdlib only; Part A production identities, Part B toy walks, kernels mode
  (limb ports vs big integers, an AArch64 interpreter that executes the shipped `secp_loops.s` on
  production points, the parsed instruction counts, the ledger recomputation). `toy_trial` takes an
  optional operation budget, which is what the census uses.
- `nu_census.py`: the load-normalized end-to-end census driver; `merge_census.py`: merges the runs,
  recomputes the pooled anchors, the Clopper-Pearson endpoints and the fit, and rewrites
  `evidence/nu-census.json` and the claim text.
- `evidence/nu-census.json`: every census cell, pooled anchors, fit residuals, repair-rate and
  trail-length sensitivities, budget assertions, gate rows.
- `kernels/`: the counted AArch64 loops (`secp_loops.s`), their generator (`gen_loops.py`), the
  freestanding C specification (`secp_walk.c`), clang's field-kernel output (`secp_fe.s`, the
  inversion bound) and `strip_asm.py`.

## Commands

Hosted initial input (about 8 s):

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

Prove the bound with exact integers:

```
python3 certificate.py
```

Regenerate the shipped assembly:

```
python3 kernels/gen_loops.py elf > kernels/secp_loops.s
```

Run the census (one group; the full ladder is six such runs):

```
python3 nu_census.py '{"groups":["q25"],"dp_bits":{"q25":4},"trials":512,"gate":true,"gate_groups":[0],"nu_mult":[0.4,0.7,1.0,1.3],"R_table":256,"seed":20261024}'
python3 merge_census.py   # regenerates evidence/nu-census.json and the claim's numbers
```
