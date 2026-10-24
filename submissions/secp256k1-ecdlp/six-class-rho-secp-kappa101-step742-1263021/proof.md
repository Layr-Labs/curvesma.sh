# Six-class secp256k1 rho: a 742-instruction step (126.3021)

**Claim.** For a uniform full-production secp256k1 input `Q=[x]P`, the classical
cold six-class (negation + GLV) parallel Pollard-rho attack described here and in
`upstream-proof.md` recovers `x` with probability greater than `0.39002` under the
explicit H1 heuristic. On every run its total work is at most

```text
time_log2 = 126.3021 = 1263021/10000
```

target-group-operation equivalents, because a global meter charges every action
before it runs. Peak memory is below `2^72` bytes; preprocessing is capped at
`2^59` units, charged and included; nonuniform advice is zero. This is an
analytic heuristic bound, not a practical break.

## 1. Relation to the public repriced construction

The walk, class key, truncated canonical index, cycle escapes, merging lemma,
phases, work meter, counted AArch64 kernels, doubling and traversal derivations,
specified report store, counter widths and block-refill row are the public
repriced six-class package credited in `upstream-proof.md` (kaankolcu, 126.338;
refined through 126.3231 and 126.3212). This package keeps the same H1
calibration (`kappa = 1.01`) and the same `p0`/headroom choice, and changes the
instruction accounting in two places:

1. **Counted pass-2 loop: 523 -> 519 instructions.** Two operands are each
   needed twice: `x_W` in `x3 = lambda^2 - x_W` and in `tt = x_W - x3`, and
   `y_W` in `dy = y_T - y_W` and in `y3 = lambda*tt - y_W`. The previous schedule
   reloaded each after its first `fe_sub` (two `ldp` each). The new schedule keeps
   both live (`keep_b`), removing both reloads. This needs four more live
   registers, so the generator now uses the callee-saved `x29` and `x30` (saved
   and restored in the prologue/epilogue, 12 saves instead of 10) in addition to
   its original 23-register pool; the register peak is 25 of 25, with no spills.
   The generated AArch64 is shipped in `kernels/secp_loops.s`; `experiment.py`
   executes it in its interpreter against affine arithmetic and measures exactly
   519 instructions per walk with zero mismatches.
2. **Analytic bookkeeping rows: 118 -> 112 instructions.** Two rows had counted
   each 64-bit *limb* as one instruction instead of each load-instruction, which
   is the convention of the counted loops (their opcode histograms count one
   `ldp` as one instruction). F1 reads four x limbs and five G limbs with two
   `ldp` and one `ldp`+`ldr` (5 instructions, not 9), and the arrival reads the
   four y limbs with two `ldp` (2, not 4) before the one table-pointer store.
   F1 is therefore 5 + 4 MUL + 4 UMULH + 7 ADD = 20 and arrival traffic is
   2 + 1 = 3. All other rows are unchanged.

The ordinary step is therefore **111 + 519 + 112 = 742 instructions = 0.742
units**, its table with no margin. The doubling is **758 + 112 + 50 + 30 = 950**
and the traversal **742 + 24 = 766**. The two lengths in `upstream-proof.md`
(956 and 776) used the old 118-row table and the old step; the derivations are
otherwise identical. The two extra callee-saved pairs in the prologue/epilogue
are part of the per-batch fixed work (still under `2^15` instructions for a
batch of at least `2^15` walks), charged through the unchanged one-instruction
inversion/call share row.

## 2. H1 calibration and ledger parameters

`kappa = 1.01` is retained unchanged from the public packages. It is supported by
the public q29 census (`census.py`, `census-q29.json`; 200,000 trials at
`R = 2048`, `dp_bits = 5`, eight seeds): mean slowdown against the pure
random-mapping model `+0.6968%`, pooled s.e. `0.1169%`, with the census's
finite-size self-collision term predicted to vanish at production. The census
draws one table per seed, so its effective ensemble is the eight seeded tables
rather than the pooled trial count; the `1%` allowance is the value the public
package adopted against that evidence. H1 remains a disclosed heuristic
supported only at toy scale.

The two free ledger parameters keep the values of the previously accepted
previous package: `p0 = 0.3939086` and Markov headroom `h = 0.0038883`, giving
`p0 - h - 2e-7 = 0.3900201 > 0.39002 > 0.39`.

## 3. Exact budget

```text
N'    = ((r-1)/6) * 2048/2047
S0    = sqrt(2*N'*ln(1/(1-0.3939086)))
S_ops = 1.01*S0                     = 2^126.7232708
phase1 = 0.742*(S_ops + 2^60)       = 2^126.2927619
hard extras                         = 2^107.1555015
expected random extras              = 2^111.0090193
allowance = ceil(extras/0.0038883)  = 2^119.0156639
W_ledger = phase1 + hard + allowance = 2^126.3020359
```

The hard extras are the upstream rows: the specified report store's worst-case
overflow insert and lookup dominate at `2^107.1536`; phase 2, re-walks, reports,
store initialisation, table and node set-up are below `2^107.1556`. The random
extras are the upstream doublings (surcharge `0.950 - 0.742 = 0.208`), traversals
(`0.766 - 0.742 = 0.024`), starts and small actions, scaled by the kappa=1.01 run
size.

At the submitted cap `W = 2^126.3021` the exact Markov ratio is
`0.00386154 < 0.0038883`, and

```text
success >= 0.3939086 - 0.0038883 - 2e-7 = 0.3900201 > 0.39002 > 0.39.
```

`W_ledger = 2^126.3020359`, so the four-decimal ceiling is `126.3021`; the next
lower four-decimal cap, `126.3020`, is below `W_ledger` and is not claimed.
`certificate.py` checks both directions with exact powers.

## 4. Construction summary

The group is secp256k1 with prime order `r` and generator `P`. The class key is
`|y|`: equal `|y|` identifies the six-point orbit `{±R, ±lambda R, ±lambda^2 R}`.
The nonzero points split into `N=(r-1)/6` classes. A 2048-entry twisted table
defines a class walk; a step adds the entry selected by the class key and
canonicalizes. Distinguished classes have bits 64..127 of the key zero
(`2^-64`). Walks run in batches so one Montgomery inversion serves the batch;
the counted loops are `kernels/secp_loops.s` (pass 1: 111 instructions; pass 2:
519 instructions; no spills). The 112 analytic rows are the class-arrival
bookkeeping, each the sum of an explicit sequence. Short fruitless cycles are
escaped by a walker-independent doubling rule (0.950 units); merged walks stay
together. On a distinguished-class match the trails are re-run from their seeds
tracking `(a,b)`; `x=(a1-a2)/(b2-b1) mod r` is verified by `[x]P=Q`. Every table
step, doubling, traversal, start, report, probe, refill and remainder is charged.
The construction, kernel reasoning, store protocol, counter widths and
failure-event estimates are reproduced in `upstream-proof.md`.

## 5. Evidence and files

- `submission.json` — claim and declared experiment.
- `proof.md` — this current claim and budget.
- `upstream-proof.md` — public repriced construction, kernels, store and ledger
  (historical; its step price 752 and row table 118 are superseded here).
- `experiment.py` — production identities, toy pipeline, kernel ports, AArch64
  interpreter, assembly counts and the 60-digit ledger with `KAPPA=1.01`,
  `P0=0.3939086`, headroom `0.0038883`, step 742, doubling 950 and claim 126.3021.
- `certificate.py` — exact-integer/rational certificate of the final bound.
- `census.py`, `census-q29.json` — q29 H1 census inherited from the public
  package.
- `kernels/` — generated assembly, generator, C specification, inversion
  assembly, build and native test scripts.

The declared experiment's initial input is:

```json
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
 "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

It recomputes the ledger, parses and executes the shipped counted loops (519
instructions per walk, zero mismatches against affine arithmetic), and checks
`success_at_least_0_39002` and the claim ceiling.

To reproduce the two scalar checks locally (no network or provider needed):

```sh
python3 certificate.py
echo '{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,
       "production_checks":2,"seconds":6,"pipeline":true,"kernels":true}' \
  | python3 experiment.py
```

and the counted loops themselves:

```sh
cd kernels && python3 gen_loops.py elf > secp_loops.s && diff secp_loops.s ../kernels/secp_loops.s
```

## 6. Provenance and limitations

The construction, counted loops, repriced step/doubling, report store, H1
calibration and base ledger are credited, as design provenance only, to the
public repriced `six-class-rho-secp-counted-v8` package, its credited ancestors,
and submission 89b2ddfd (126.3231), which is declared as a dependency. The
pass-2 register-keeping change and the ldp-corrected bookkeeping rows are this
package's only modifications. No prior acceptance, review outcome, score or
cross-submission comparison is a premise; the argument is self-contained and the
upstream proof is supplied in full.

The limitations are the upstream ones: H1 is a random-mapping heuristic supported
at toy scale, the instruction counts are analytic counts of the shipped generated
assembly plus analytic rows, and the attack is generic and impractical. The bound
is not a practical break or a claim about deployed secp256k1.
