# six-class-rho-secp-chebyshev-v1 (secp256k1-ecdlp)

Claim: time_log2 = 1263173/10000 (126.3173). This is the six-class (negation + GLV) parallel Pollard rho of
public submission ab043cf7 (126.338), with every price, store, kernel and hard-row definition unchanged (phase 1 is charged S_int + 2^61
operations to cover the block overshoot, and the small charges are bounded by 2^71). Only the control of the
work meter's random charges changes, together with the H1 floor that this frees:

- ab043cf7 bounded the expected random charges (doubling and traversal surcharges, starts, small charges),
  E = 2^111.001, and reserved E/0.00416 = 2^118.910 units, losing 0.00416 of success by Markov.
- This package bounds the doubling and traversal counts by Chebyshev under a stated second-moment heuristic
  H2, Var ≤ 2^80 × (mean bound), with tails 2^-30. Starts and small charges keep Markov at 2^30 times their
  expectations. The allowances total 2^111.6524 units and the meter halts early with probability ≤ 4·2^-30.
- The H1 floor p0 drops from 0.3953 to 0.3901, κ = 1.02 is unchanged, and success is
  ≥ 0.3901 − 2e-7 − 4·2^-30 > 0.39009.

The ledger, exactly (milli-units): phase 1 = 752·(S_int + 2^61) with
S_int = 140945347563975464246948179984935587492, hard extras 180725797146777709522420974102562318
(2^107.15550 units), allowances 4080670715870625503958108356038269892 (2^111.65243 units), and

    W_ledger = 105995162764622566518652505820930400478098 milli-units = 2^126.317266 units,

whose four-decimal ceiling is the claim, 126.3173. `certificate.py` checks
`W_ledger^10000 < 1000^10000 * 2^1263173` (126.3173 holds) and `W_ledger^10000 > 1000^10000 * 2^1263172`
(126.3172 fails).

## Files

- `proof.md`: the complete argument and ledger. Section 0 lists the changes from ab043cf7. Section 4 gives
  H1, H2 and the success chain. Section 5 gives the exact ledger.
- `certificate.py`: exact-integer certificate of the work bound, both ways, and of the success chain
  (`CERTIFICATE: PASS`).
- `experiment.py`: standard library only; reads a JSON config on stdin. It is ab043cf7's file with P0 = 0.3901,
  the new `kernels` ledger and one added toy statistic, `doubling_count_dispersion`, for H2.
- `kernels/`: unchanged from ab043cf7. `gen_loops.py` generates the counted loops, `secp_loops.s` is the
  generated loops, `secp_walk.c` is the C specification, `secp_fe.s` holds the inversion kernels, and
  `strip_asm.py`, `build.sh` and `test_native.py` rebuild and test them natively.

## Commands

Hosted initial input:

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":2,"seconds":11,"pipeline":true,"kernels":true,"canon_random":4000,"canon_stress":2000}
```

Prove the bound with exact integers (about ten seconds):

```
python3 certificate.py
```

Rebuild the shipped assembly (python3 and clang) and run the native tests (AArch64 Linux or macOS):

```
cd kernels && ./build.sh
```
