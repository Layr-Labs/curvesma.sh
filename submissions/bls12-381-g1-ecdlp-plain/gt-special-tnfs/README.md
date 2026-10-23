# bls12-381-g1-ecdlp-plain: GT transfer with special tower NFS

Candidate for the plain BLS12-381 G1 discrete-log track. It maps
`T = [x]P0` to `z = e(T,Q0) = e(P0,Q0)^x` with the pinned public pairing,
solves two discrete logarithms in the order-`r` subgroup of `F_{p^12}` with the
published SexTNFS/STNFS-JP estimates for this prime, and verifies `[x]P0 = T`.

Status: ready for independent review. The submitted bound is analytic and
literature-grounded; there is no production NFS execution. See `proof.md` for
the construction, the exact ledger, the success argument, the memory bound and
the limitations.

Files:

- `submission.json` - scored claim and the declared experiment.
- `proof.md` - construction, correctness, accounting and references.
- `experiment.py` - BLS12-381 constant checks, exact ledger/score arithmetic,
  and a toy factor-base index-calculus pipeline (20-32 bits).

Run the experiment with standard input as JSON, for example:

```sh
echo '{"bits": 28, "seed": 7}' | python3 experiment.py
```

The experiment is toy-scale evidence only. It does not execute a production
NFS, the pinned pairing, or the submitted production bound.
