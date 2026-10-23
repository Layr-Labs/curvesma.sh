# gt-stagger-rho-v2

Fixed-argument pairing inversion on BLS12-381 G1 into G2
(track `bls12-381-fixed-g1-invert-g2`).

* Claim: **127.68** (`12768/100`) log2 reference dense-`Fp12` multiplications.
* Success: **>= 0.39** over uniform `P`, uniform `z`, fresh coins.
* Peak memory: **226,736 B**.
* Detector: staggered checkpoints `(k,W) = (4,8)`, 12 simultaneously live
  checkpoints (enumerated, not asserted), hard cap
  `B = ceil(1.16*sqrt(r)) = 265627019777073510555230149274057194865`.

Supersedes `gt-stagger-rho-v1` (`59149222`, claim 127.67), rejected on record
widths and on model claims that were asserted rather than implemented. `proof.md`
§8 maps each objection to its fix.

## Files

| file | what it is |
| --- | --- |
| `submission.json` | the claim, the score and the declared experiment inputs |
| `proof.md` | reduction, walk, detection law, success model, exact ledger, memory |
| `note.md` | the public short version of the same |
| `experiment.py` | standard-library reimplementation of every number, ~5 s |
| `results.json` | captured output of `experiment.py` at the declared inputs |

## Run it

    echo '{"production_checks":2,"table_checks":2,"trials":300,"toy_bits":26,"R":256,"cap_factor":1.16,"seconds":12,"seed":1,"model_draws":20000}' | python3 experiment.py

Sections: A production facts on BLS12-381; B1 live-count enumeration, B2 detection
law versus direct simulation, B3 complete toy search at the declared `R`, B4
jump-count sensitivity; C the exact cost and memory ledger; D1 quadrature of the
random-map model, D2 Monte-Carlo cross-check, D3 dyadic-alignment sweep, D4
success at the declared cap.
