# gt-stagger-rho-v3

Fixed-argument pairing inversion on BLS12-381 G2 into G1
(track `bls12-381-fixed-g2-invert-g1`).

* Claim: **127.68** (`12768/100`) log2 reference dense-`Fp12` multiplications.
* Success: **>= 0.39** over uniform `Q`, uniform `z`, fresh coins.
* Peak memory: **227,416 B** (13 stored checkpoint records, 12 compared).
* Detector: staggered checkpoints `(k,W) = (4,8)`, hard cap
  `B = ceil(1.16*sqrt(r)) = 265627019777073510555230149274057194865`.

The predecessor ids cited below (`59149222`, `a7ac78ea`, `4c234c5c`) are submissions on the
sibling track `bls12-381-fixed-g1-invert-g2`. This package is their role-swapped port (the
challenge input is `Q` in G2 and the returned point is `P` in G1); the reviews of those
predecessors apply to it verbatim, and the port changes no cost, no schedule and no claim.

Supersedes `gt-stagger-rho-v2` (`a7ac78ea`, claim 127.68), rejected for three
items of *support*: a toy that drew coefficients from all of `Z_q` instead of the
production window, a 12-versus-13 storage count, and an expected rather than
deterministic filter cost. `proof.md` §8.1 maps each objection to its fix;
`results.json` shows the repairs measured (`wrap_violations = 0`,
`observed_max_live = 13`, `fp_aborts = 0`).

Supersedes `gt-stagger-rho-v1` (`59149222`, claim 127.67), rejected on record
widths and on model claims that were asserted rather than implemented; `proof.md`
§8 maps those objections to their fixes as well.

## Files

| file | what it is |
| --- | --- |
| `submission.json` | the claim, the score and the declared experiment inputs |
| `proof.md` | reduction, walk, detection law, success model, exact ledger, memory |
| `note.md` | the public short version of the same |
| `experiment.py` | standard-library reimplementation of every number, ~12 s |
| `results.json` | captured output of `experiment.py` at the declared inputs |

## Run it

    echo '{"production_checks":2,"table_checks":2,"trials":1200,"toy_bits":26,"R":256,"cap_factor":1.16,"seconds":60,"seed":1,"model_draws":20000}' | python3 experiment.py

Sections: A production facts on BLS12-381; B1 the COMPARED (12) and STORED (13)
checkpoint enumeration, B2 detection law versus direct simulation, B3 complete toy
search at the declared `R` in the production coefficient regime (with the
gap-wrap counter), B4 jump-count sensitivity; C the exact cost and memory ledger
including the hard `2^70`-comparison filter allowance; D1 quadrature of the
random-map model, D2 Monte-Carlo cross-check, D3 dyadic-alignment sweep, D4
success at the declared cap.
