# Cheon d=2^20 with fixed-base windowed BSGS

Candidate for `bls12-381-g1-ecdlp`. It uses only the supplied power
`A_d = [s^d]P0` with `d = 2^20` (the largest divisor of `r-1` not exceeding
`n = 2^20`) and applies the `p-1` case of Cheon's two-phase attack, replacing
each known-scalar exponentiation by a three-window fixed-base lookup plus two
additions.

Claimed bound, argued in `proof.md`:

| Quantity | Value |
|---|---|
| `log2 W(2^20)`, `W = H + R/2000` | `119.906712` |
| peak memory | `2^125` bytes (chained layout: `2^124.884305`, 128-byte nodes, 16-byte heads) |
| preprocessing (included in `W`) | `2^111.034217` group-op equivalents |
| nonuniform advice | `0` bytes |
| success probability | `>= 0.85` (proved `>= 0.874608`, Markov at the global `12*N` link budget) |

**Matching.** The baby points sit in a separate-chaining table of `2^119` buckets with the
concrete 2-universal multiply-shift hash `h(x) = ((a*x+b) mod 2^384) >> 265` (`a` odd,
`a,b` fresh per attempt). One global counter accumulates every link the attempt walks and
returns FAILURE when it would exceed `12*N`; no per-operation cap is used, so the matching
work is bounded by construction and the abort event is exactly the Markov event. The
failure probability is a first-moment (Markov) bound only: key pairs collide with
probability at most `2/2^120`, `E[links] <= 1.504704*N`, and
`Pr[cap exceeded] <= 1.504704/12 = 0.125392`.

**Scope.** The bound is analytic. No production-size execution is claimed, no
production-size tables are built, and the toy experiment is evidence for the
algebra and the integer ledger only, not for production cost or memory.

## Files

- `submission.json` — scored claim and the declared experiment.
- `proof.md` — Cheon derivation, fixed-base window construction, matching
  argument, complete integer ledger, memory/preprocessing/advice accounting,
  rejected alternatives and references.
- `experiment.py` — stdlib-only toy two-phase run (seeded random and boundary
  scalars) plus production parameter, ledger, factorization, primitive-root and
  digit-structure assertions.

## Reproduce

```sh
echo '{"trials":4,"seed":20261022,"d_bits":10}' \
  | python3 experiment.py
```

The experiment prints one JSON object with the toy prime, recovered-trial
count and the recomputed production ledger. It asserts that all trials recover
the secret, that the factorization of `r-1` is exact, that `7` is the smallest
primitive root modulo `r`, that the 3x85-bit comb has no incomplete-addition
exception on boundary and seeded sample scalars (the universal `0 < x < r`
statement is the analytic argument in `proof.md`), and that the declared `W`,
preprocessing and memory bounds are not exceeded.

Repository tests: `python3 -m unittest discover -s tests` (needs
`tests/requirements.txt`, including `py_ecc==8.0.0`).

## Provenance

Built from the track's `rules.md`, `COST_MODEL.md`, `classical-group-work-v2`
and the J. H. Cheon, "Discrete Logarithm Problems with Auxiliary Inputs",
J. Cryptology 23 (2010) identities. No prior submission is reused.
