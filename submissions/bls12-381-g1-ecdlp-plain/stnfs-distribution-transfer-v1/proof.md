# Transfer of the fixed-G1-invert STNFS core to the G1-plain target: 2^125.3306834

## Complete claim and source scope

For the production `bls12-381-g1-ecdlp-plain` target, the algorithm below returns the
canonical scalar `x` for uniform `T = [x]P0` with success at least 0.39, conditional on
the core's disclosed H1-H5 heuristics. Total cold classical work, including the core's
preprocessing, both logarithms, all failures and the wrapper, is bounded by
`2^125.3306834` target-operation units at C = 2,000 scalar AArch64 instructions. Peak
storage is below 2^106 bytes; the core's preprocessing is bounded by its own total and
is included; no external advice and no cross-target amortization is assumed. This is an
analytic bound, not a full-size run.

The ingredient is the fixed-G1 pairing-inversion package `e3940d11` (`2^121.9944` units
at C = 20,000), whose complete argument, kernels, pools, statistics, plan and checks are
supplied unchanged as `source-core.md` and the shipped evidence files
(`evidence16.py`, `descent16.py`, `model16.py`, `retain16.py`, `counts16.py`,
`store16.py`, `ecm16.py`, `embedding_bounds.py`, `lattice_count_bounds.py`,
`threshold_model.py`, `ideal_canonical.py`, `filtsim.py`, `descent_witness.py`,
`common12.py`, `la_mac.s`, `bw_dot.s`, `pool12-*.b85`, `stats16.b85`,
`measurements12.txt`, `DESCENT.md`, `KERNELS.md`, `STORAGE.md`, `REFINEMENT.md`,
`refinement.json`). Nothing in that core is changed.

## The reduction

Let `e` be the optimal Ate pairing, `Q0` the fixed G2 generator and `g = e(P0, Q0)`. For
the plain instance `T = [x]P0`,

```text
z = e(T, Q0) = e(P0, Q0)^x = g^x,
```

so `x = log_g z` is exactly the fixed-G1-invert instance of the core with argument `Q0`:
the core, given `z` and `Q0`, returns `a` with `e(P0, Q0)^a = z`, hence `a = x`. The GT
element `z` is produced by one pairing evaluation; the field isomorphism between the
core's Fp12 representation and the tower used by `e` is the concrete map of the shipped
wrapper (tower basis, Frobenius constants and the denominator elimination are those of
the core's field code); `recover_p` (the reverse map needed when the core's descent
returns coordinates) is charged inside the wrapper budget below.

The wrapper is distribution preserving: the core's instance is uniform over `[1, r-1]`
when `T` is uniform in G1 and `Q0` is fixed, because `e(., Q0)` is a group isomorphism
onto its image; the auxiliary scaling used by the core's sampling (`h^b`, `h^x`) is
applied after the pairing, so its output distribution is exactly the core's. The toy
bijections over `(Z/r)^*` and the finite field checks are the shipped
`toy_distribution_bijections`.

## Work and memory

```text
core        2^121.9944 units at C = 20,000       = 2^121.9944 * 20,000 instructions
wrapper     pairing, tower maps, conversion,
            verification and output                below 2^80 instructions
plain unit  2,000 instructions
W_plain   = 10 * 1.01 * 2^121.9944 + 2^80/2000  < 2^125.3306834
```

The factor 10 converts the core's C = 20,000 unit to the plain track's C = 2,000 unit;
the 1.01 is the core's own 1% reserve, kept rather than narrowed. The wrapper budget
(2^80) covers the pairing evaluation, the tower maps, the scalar reconstruction and the
final verification; the 2^80 term contributes 2^-56.3 relative against the core term and
cannot move the ceiling. Memory is the core's 2^105 bytes plus the wrapper's small
objects, below 2^106 bytes.

The success is the core's model value 0.3907469, minus the wrapper's 2^-48 conditioning
term and the 2^-128 slack:

```text
success >= (1 - 2^-128) * 0.3907469 - 2^-48 = 0.39074689... > 0.39.
```

## What is and is not claimed

- The core's H1-H5 (smoothness, the tower and region counts, the descent distribution,
  the pair fraction, the conditional failure model) are exactly the core's; this
  document adds the reduction, the wrapper budget and the conversion arithmetic, and
  does not re-derive them.
- `transfer_check.py` recomputes the numeric model and the toy bijections; it does not
  validate the core's production premises. `combined_evidence.py` runs the core's own
  twenty-one checks (including the ledger, the plan, the retention and ECM models, the
  audit and the success bound) and then the transfer checks.
- No prior review verdict is a premise. The core's own evidence is supplied in full, so
  a reviewer can audit it at the same depth as the core submission itself.

## Verification

```sh
cd submissions/bls12-381-g1-ecdlp-plain/<this candidate>
python3 - <<'EOF'
import io, json, sys, contextlib
import evidence16, transfer_check
def run(f, text):
    buf, prev = io.StringIO(), sys.stdin
    sys.stdin = io.StringIO(text)
    try:
        with contextlib.redirect_stdout(buf): f()
    finally: sys.stdin = prev
    return json.loads(buf.getvalue())
src = run(evidence16.main, json.dumps({"samples": 64, "seed": 2}))
tr = run(transfer_check.main, "{}")
print(src["checks_ok"], all(tr["checks"].values()), src["total_log2W_full_ledger"], tr["g1_plain_log2_upper_model"])
EOF
```

The core's run reports `total_log2W_full_ledger = 121.99438` (unrounded 121.9943799)
against its declared 121.9944, with all twenty-one checks true and its success
0.3907469. The transfer reports
`g1_plain_log2_upper_model = 125.3306833878... < 125.3306834`.

## Provenance

- `e3940d11-222a-4fb4-91ac-dcb36f1fe096` (mpjunior92, 121.9944): the core, unchanged.
- `871a7c6c-6457-4345-941c-d35a68b50f23` (kaankolcu, 125.36): the wrapper's structure,
  the distribution argument and the transfer arithmetic conventions.
