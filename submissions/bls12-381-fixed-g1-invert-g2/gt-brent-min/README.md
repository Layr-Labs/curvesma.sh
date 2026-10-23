# gt-brent-min

`bls12-381-fixed-g1-invert-g2`: O(1)-memory **Brent** rho in GT for fixed-G1
pairing inversion. Claim `log2 <= 12815/100 = 128.15` in dense-Fp12-multiplication
equivalents (`C = 20,000`), success `>= 39/100`.

One change from the structure the reviewer accepted on `e444842b` (plain Floyd
rho in GT, 129.20): Floyd's tortoise-and-hare is replaced by Brent's cycle
detection with a **single saved point**. No fingerprint, checkpoint grid,
comparison window, distinguished-point set, dictionary, cycle-length bound or
positivity lemma.

| file | content |
| --- | --- |
| `submission.json` | the claim, the score `2563/20`, the evidence list |
| `proof.md` | reduction, the loop, the exact detection law (with proof and its verification against a direct simulation), the cap derivation, the ledger, the success accounting, the toy, memory, disclosures |
| `note.md` | the same argument in note form, with the numbers and the credits |
| `experiment.py` | the runnable evidence (standard library only, about 8 s) |
| `results.json` | the output of `experiment.py` on the declared inputs |

Key numbers: cap `B = ceil(1.60 sqrt(r)) = 2^128.106617`; one full element
comparison per step at 288 instructions, 592 extra instructions per step in
total; main term `2^128.148701`; table build `2^30.0004`; total `2^128.148701`;
claim `128.15`. Model success `0.4534` nominal, `0.4097` at the worst dyadic
alignment; zero-denominator charge `2^-126.75`. Memory `671,156,304` B.

Reproduce:

    python3 experiment.py < inputs.json

with `inputs.json` = `{"production_checks": 2, "table_checks": 2,
"gt_walk_checks": 1, "trials": 600, "toy_bits": 25, "R": 1024, "seconds": 60,
"seed": 1, "model_draws": 2000}`.
