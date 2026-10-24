# argument (v34) -- repairs verdict v33 (rejection of gt-instantiable-40pct-v33, ticket ec9f5ded)

Claim **127.7230** (reduced 127723/10000). W = 5616608965314120933304243253989908740895226, summed in one function from the thirteen rows of
`evidence/ledger34.json`. The rejection had six premises; each is answered by code, not prose.

## 1. The import and the literals

`scripts/toy28.py` is present in this tree and defines `Tower.keys12`; `scripts/protocol28.py` imports it and
`evidence/protocol34.json` is the run of that import plus the key-construction smoke. The v33 tree dropped
`toy28.py` while `protocol28.py` still imported it, which is why the reviewer saw an absent module. Every
2^k token in the claim and ledger now maps to one owning definition in `evidence/ledger34.json`
(`store_spec_to_code.production_lambda` owns 2^128 and nothing else; no 2^120 token remains).

## 2. The store is the shipped store

`dict28.c` is a bucket array of `u64 nb` pointers, `ws_bucket` returns `u64`, chains are capped at `PROBES`
= 4 by FIFO eviction of the tail, and capacity is `4*nb` by construction. It has no 2^128-entry directory
and no lossless overflow segments. v34 therefore **deletes** the rows that priced them -- the 2^128-entry
directory, the 2^135.75 B footprint and the 1.490e39 Ir initialisation -- instead of restating them at a
smaller exponent. What remains is the executed per-access price of the shipped symbols (959 Ir lookup, 1711
Ir insert, 681 Ir slot export, all censused) and a retention pair at the production lambda:

- derived: 1 - P(Poisson(lambda) >= 5) = 0.9999883322852278 at lambda = S/2^128 = 0.2815416777;
- measured: 0.99995911704498 over 1485 trials of the shipped four-probe FIFO policy at that same lambda (n = 1153 arrivals,
  nb = 4096, se 4.840872807401776e-06);
- delta -2.9215240247792984e-05, explained by within-bucket evictions during fill that a terminal-size limit cannot see;
- control: the same arm at lambda 0.2417 gives 0.9999795939189878 and does not match, so v33's lower-occupancy arm fails as
  a transfer control for the right reason.

Retention is reported, not relied on: detection is per-worker trail comparison, so an evicted record cannot
lose the walk. The v33 figure 0.998435 is withdrawn as a store-retention number -- it was the arm's
1276/1278 success-resample fraction -- and the v33 prose value 0.998653 matched neither arm.

## 3. The modular multiplication

The v33 179-Ir disassembly priced undefined behaviour: its object declared `t[6]` and wrote `t[5..8]`, with
the reduction reading `t[6]` and `t[7]`. That row is **withdrawn as void** and no row here cites it. The
coefficient update is charged from the tree's measured Fp span instead: `MUL_CHARGE` 47924 over 54 Fp multiplies
per Fp12 = 888 Ir per Fp multiply, so two 255-bit mulmods plus the 12-lane load cost 1788 Ir/step. That is a
ceiling of a measured quantity, not a disassembly of a broken one.

## 4. The trap bounds the event it charges

The stall counter increments on `br OR guard`, and `Group.step` squares only on guard, so an abort need not
involve 32 squarings. v34 takes the census predicate as the event: `P(br OR guard) <= P(br) + P(guard)`,
with `P(guard) <= C(12,2)/2^384 = 66/2^384` from class-membership birthday structure and `P(br) <= p_s^32` over the
32-step stall window. The measured union rate (4497 admissions in 5040000 steps) is reported as corroboration
*below* the derivation, not as its source. Squaring-only-on-guard is priced as its own sub-event -- the
closed-orbit window -- with charge 131300296263244672882090438236340667904 Ir.

## 5. Success transfer

The transfer arm is run at the production lambda, not a lower one. The reported p5/median keep the retention
factor out of the success resampling because no success premise is imported from occupancy: 0.578125 / 0.625, floor
0.39.
