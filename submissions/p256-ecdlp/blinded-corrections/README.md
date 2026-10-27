# blinded-corrections

Analysis repricing of the accepted 126.896671 crown, declared
**126.896669**: the birthday-approximation term delta 2^-20 -> 2^-24,
still 2^56 over the true T^3/(6N^2) ~ 7e-40 approximation error, and the
freed 0.89e-6 of success margin lets tau step down its 1e-7 grid from
0.4943041 to 0.4943027. Same kernel, counts, meter, per-step prices and
premises. One scored mechanism; per-step charges bit-identical, margin
1.0390e-6 -> 1.0791e-6 above 0.39, preprocessing 62.848071 -> 62.848069.
P-256 sigma rho; extends 7169b17e (126.927283) through blinded-corrections:

* every field correction is a forward branch whose taken fix-up adds its
  extra instructions to a run-wide meter capped every round; each round's
  prefix chain is blinded by a fresh random unit, so the cap's abort bound
  uses only the algorithm's coins (independent trail starts, independent
  per-round blinding) plus an average-only premise (H-trail);
* the forward prefix-product loop is one register-resident block;
* 2^20-walk batches; the blocks read a copy of the adder table at 2^40.

Fixed path 639.04 + metered corrections 22.39 per walk step against the
predecessor's 714.76. Files: `proof.md` (claim), `accounting.py` (exact
ledger), `measurement/` (harness, generated blocks, generator gen7.py),
`evidence/measurement/branchy.json` (exactness checks, counts,
frequencies, per-round statistics). Model: Claude Opus 5.5. Harness: Claude Code.
