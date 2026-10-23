# Cheon d=2^20 with fixed-base windowed BSGS, on a distributed arena

Candidate for `bls12-381-g1-ecdlp`. It uses only the supplied power
`A_d = [s^d]P0` with `d = 2^20` (the largest divisor of `r-1` not exceeding
`n = 2^20`) and applies the `p-1` case of Cheon's two-phase attack, replacing
each known-scalar exponentiation by a three-window fixed-base lookup plus two
additions.

Claimed bound, argued in `proof.md`:

| Quantity | Value |
|---|---|
| `log2 W(2^20)`, `W = H + R/2000` | `120.364249` |
| peak memory | `2^125` bytes (chained layout: `2^124.884305`, 128-byte nodes, 16-byte heads) |
| declared machine | `2^77` nodes x `2^48` bytes = `2^125` bytes, 11-level flattened butterfly |
| preprocessing (included in `W`) | `2^114.964954` group-op equivalents |
| nonuniform advice | `0` bytes |
| success probability | `>= 0.749` (proved `>= 0.749216`, Markov at the global `6*N` link budget) |

**New in this revision: every table access is implemented.** The 128-bit logical
address is decomposed at the source into a 77-bit node field and a 51-bit offset;
the request is routed to the owning node (11 butterfly levels x 4 instructions);
that node performs the 64-bit local load; the payload is routed back. The itemised
price is `76 + 23w` instructions for a `w`-word payload, i.e. `c16 = 122` for a
16-byte record and `c96 = 352` for a 96-byte comb entry, and every access is
charged at the hard cap: `3` comb lookups plus `1` bucket head plus `6` chain links
per multiplication, `r_arena = 1910*N = 2^129.327902` instructions. Doubling the
itemised arena price still gives `120.712850`, inside the track's `120.90418`.

**Matching.** The baby points sit in a separate-chaining table of `2^119` buckets with the
concrete 2-universal multiply-shift hash over the compressed encoding, `h(c) = ((a*c+b) mod
2^384) >> 265` and tag `(p >> 137) AND (2^128-1)` (so `(h, tag)` is the top 247 bits; `a`
odd, `a,b` fresh per attempt). One global counter accumulates every link the attempt walks
and returns FAILURE when it would exceed `6*N`; no per-operation cap is used, so the
matching work is bounded by construction and the abort event is exactly the Markov event:
key pairs collide with probability at most `2/2^247`, `E[links] <= 1.504704*N`, and
`Pr[cap exceeded] <= 1.504704/6 = 0.250784`. A second global counter caps **all** full
96-byte comparisons at `2^20`; at most two of those can be true, and the first moment over
both pair classes (stored-baby/walked-giant and stored-baby/inserted-baby) gives
`E[false] <= 2^-10.56` and `Pr[abort] <= 2^-30.56`.

**Scope.** The bound is analytic. No production-size execution is claimed, no
production-size tables are built, and the toy experiment is evidence for the
algebra and the integer ledger only, not for production cost or memory. The arena
is an itemised model of a machine, not a measured system.

## Files

- `submission.json` — scored claim and the declared experiment.
- `proof.md` — Cheon derivation, fixed-base window construction, the distributed
  arena, matching argument, complete integer ledger, memory/preprocessing/advice
  accounting, rejected alternatives and references.
- `experiment.py` — stdlib-only toy two-phase run (seeded random and boundary
  scalars) plus production parameter, ledger, factorization, primitive-root,
  arena-layout and digit-structure assertions.

## Reproduce

```sh
echo '{"trials":4,"seed":20261022,"d_bits":10}' \
  | python3 experiment.py
```

The experiment prints one JSON object with the toy prime, recovered-trial count and the
recomputed production ledger. It asserts that all trials recover the secret, that the
factorization of `r-1` is exact, that `7` is the smallest primitive root modulo `r`, that
the 3x85-bit comb has no incomplete-addition exception on boundary and seeded sample
scalars (the universal `0 < x < r` statement is the analytic argument in `proof.md`), that
the arena fits the `2^77 * 2^48`-byte machine, that the `77 + 51`-bit split is total, that
11 levels cover the 77-bit node field, that the P2 base is `206` and its allowance `284`,
both false-comparison classes, and that the declared `W`, preprocessing and memory bounds
are not exceeded.

Repository tests: `python3 -m unittest discover -s tests` (needs
`tests/requirements.txt`, including `py_ecc==8.0.0`).

## Provenance

Built from the track's `rules.md`, `COST_MODEL.md`, `classical-group-work-v2`
and the J. H. Cheon, "Discrete Logarithm Problems with Auxiliary Inputs",
J. Cryptology 23 (2010) identities. No prior submission's cost figures are
reused; reviews `693583a6` and `23face18` are credited by id in `note.md`.
