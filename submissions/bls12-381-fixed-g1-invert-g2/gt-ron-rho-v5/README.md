# gt-ron-rho-v5

Fixed-argument pairing inversion on BLS12-381 G1 into G2
(track `bls12-381-fixed-g1-invert-g2`).

* Claim: **127.69** (`12769/100`) log2 reference dense-`Fp12` multiplications.
* Success: **>= 0.39** over uniform `P`, uniform `z`, fresh coins; model success
  at the cap is `0.4501` at the worst dyadic alignment, minus the charged filter
  abort `2^-12.7` and the `12*B/r = 2^-123.6` zero-denominator term.
* Peak memory: **671,164,504 B** (`2^29.32`); 13 stored checkpoint records, 12
  compared per step, `R = 2^20` jump table.
* Detector: staggered checkpoints `(k,W) = (4,8)`, hard cap
  `B = ceil(1.16*sqrt(r)) = 265627019777073510555230149274057194865`.

`v5` supersedes `gt-stagger-rho-v3` and answers the reviews of `a7ac78ea` and
`8b096f2c` at the root: the walk's coefficients are now drawn **uniformly modulo
the group order** (512 coin bits per draw, no rejection sampling, bias
`< 2^-256`), the same regime as the accepted entry `e444842b`, instead of the
`[1, 2^64)` window that made `r/M^2 = sqrt(r)` and left the random-map model
unjustified at production coefficients. The positivity lemma that window bought
is **withdrawn** and replaced by an explicit charge on the success bound. The
jump table is `R = 2^20`, also as in `e444842b`.

## Files

| file | what it is |
| --- | --- |
| `submission.json` | the claim, the score and the declared experiment inputs |
| `proof.md` | reduction, walk, detection law, success model, exact ledger, memory |
| `note.md` | the public short version of the same |
| `experiment.py` | standard-library reimplementation of every number, ~6 s |
| `results.json` | captured output of `experiment.py` at the declared inputs |

## Run it

    echo '{"trials":300,"model_draws":20000}' | python3 experiment.py

Sections: A production facts on BLS12-381; B1 the COMPARED (12) and STORED (13)
checkpoint enumeration at every hare position to `2^20`, B2 detection law versus
direct simulation (1830 pairs, 0 mismatches), B3 the complete toy search with
coefficients uniform mod `q` (the production regime), B4_R_sensitivity the same
at `R = 1024`, B4_degenerate_branch a constructed zero-denominator collision run
end to end; C the exact cost and memory ledger including the hard
`2^80`-comparison filter allowance and the `12*B/r` charge; D1 quadrature of the
random-map model, D2 Monte-Carlo cross-check, D3 dyadic-alignment sweep, D4
success at the declared cap.
