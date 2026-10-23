# p256-ecdlp: parallel Pollard rho on negation classes with a sigma walk

Claim: for uniform x in [0,r) and Q=[x]P on NIST P-256, the algorithm below returns x
with probability at least 0.39. It uses at most **2^127.370918** target-group-operation equivalents in total (score 63685459/500000), including all preprocessing, failures,
lookup and verification. Peak memory is at most 2^76.542441 bytes. Preprocessing is 2^47.949581 units, and no nonuniform advice is
used. This is a classical cold attack on one target.

The success bound rests on the disclosed heuristic H1, the standard random-mapping model
for r-adding walks on negation classes, and on M1, the class-data premise (section 2.0 of proof.md).
M2 is a premise on the table: the table and the walk coins come from the stream random.Random(1), seed 1 being a constant of the specification chosen without reference to the target. The first draw randrange(1, r) is discarded, then each point takes two draws randrange(r), and the coefficients are treated as uniform and independent of the target. From M1, M2 and exact algebra
follow the escape, re-entry, abandonment and metered-work bounds (section 2). The cost rests on H3,
average instruction counts measured with the organizer's toolchain, and every run is
charged by a work meter with a Markov reserve. **There is no reference
operation on this track:** every instruction, including all P-256 field arithmetic and
point additions, is charged at 1/C with the borrowed normalisation C = 1000. This is a
carefully accounted generic attack. It is not a practical break.

**Distinguished points and resolution.** Walks run in trails. Each trail starts from fresh
coins and ends at its first distinguished point, which it records once together with its
coins. A key match is resolved by replaying the two trails (section 1.3, Lemmas R and D).
A fruitless cycle therefore cannot store the same point twice, and a degenerate match has
probability exactly 1/r.

**Escape and the order of 2.** An escape goes to 2·canon(M), where M is the cycle's least-key member.
The escape point's coefficients are new within its trail once ord_r(2) > 2·TRAIL_MAX (D1). For the
P-256 order this is verified in `accounting.py` (`ord_2_lower_bound`): ord_r(2) is at least S/2 for a
69-bit smooth part S of r − 1, and the 188-bit cofactor is not factored.

Files:
- `accounting.py` recomputes every number with exact integers (output in
  `evidence/accounting.json`), including the order-of-2 check.
- `measurement/` holds the per-step instruction-count harness and the field code
  (`p256_field.h`, `p256_sub.h`), plus the Python replica `walk_ref.py`.
- `evidence/measurement/` holds the production totals, the 16-step replica check and the checks of the
  restart and abandonment paths.
- `experiment.py` is the declared toy-curve experiment, running the whole procedure;
  `evidence/trails/` holds larger offline runs of the same file.

## 1. Algorithm

The algorithm is in section 1 of proof.md: classes and the sigma walk (1.2), the procedure (1.3), with
the escape rule, Lemma R (replays reach the recorded point) and Lemma D (degenerate matches).

## 2. Success probability

Section 2 of proof.md: the premises (2.0), the formal structure and D1/D2 (2.1), the degenerate-pair and
within-trail bound (2.2), the closure rates (2.3), abandonment (2.4), the metered rare work (2.5), and the
success bound (2.6). The success bound is 0.395898 after the 1/64 meter loss and the failure terms.

## 3. Cost

- The common batched path is measured with whole-run Callgrind totals on production inputs, 847 instructions per walk step.
- Rare-path operations are charged forced-path maxima (`measurement/forced/summary.json`).
- The Markov reserve (40.90 per walk step) covers the expected rare work, which is 0.639140 per walk step.

Preprocessing is the adder table and the initial walk starts, each computed with two scalar multiplications per point;
both are included in the total.

## Reproduction

```sh
python3 -I submissions/p256-ecdlp/rho-sigma2/accounting.py
echo '{"curve":"a3","bits":28,"trials":80,"seed":2801,"adder_bits":11,"dp_bits":5,"walks":16,"tau_e4":5303,"filter_bits":64,"check_all_replays":1,"escape":"double"}' | python3 -I submissions/p256-ecdlp/rho-sigma2/experiment.py
```

The forced-path summary is `measurement/forced/summary.json`; the per-step measurement and the replica are in
`measurement/`.

## Environment and provenance

- **Base.** `Layr-Labs/curvesma.sh` master @ `467130d`. Measurements ran on Apple silicon through Docker linux/arm64
  in the organizer's reference-count image (gcc 14.2.0, valgrind 3.24.0).
- **Setup.** `yukon setup` was not run locally. Experiments were run directly with `python3 -I`. Local results are
  evidence, not a score.
- **Dependencies.** None declared. No submission by anyone else was used.
- **Checks in this checkout.** `accounting.py` reproduces `evidence/accounting.json` exactly, the declared experiment
  runs locally, and the replica agrees with the kernel on every build. The pinned `llm-verifier preflight` has not
  been run on this checkout; it is run before submission.

## What did not work, and limitations

- The forced-path maxima are upper bounds for the measured build; the bounding build is checked for equality and dominance
  on the forced paths, not proved for every input.
- M1 is the premise that carries the success bound. It is not proved for this walk, and the toy runs test it at theta = 2^-5 only.
  H1 is the standard random-mapping heuristic, used for collisions.
- Practicality: about 2^68.5 distinguished points must be stored with their coins. This is not practical.
