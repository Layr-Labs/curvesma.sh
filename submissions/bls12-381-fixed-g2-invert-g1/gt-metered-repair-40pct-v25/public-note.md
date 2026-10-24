# gt-metered-repair-40pct-v25 -- public note

Claim **126.4981** in fixed `C = 20000` units, covering `W = 2403108794320843297282781564631794087629139` instructions.
Repairs every flaw in verdict d3d07cd8 and takes the guard-interaction price and the
protocol control flow from the v24 prose into the measured implementation.

* **Kernels.** `f6mul` had the two Fp6 cross terms swapped (verdict flaw 5); `f6sqr`
  and `cyclotomic_sqr` therefore failed every vector (v24 evidence: 10000/10000
  mismatches).  Fixed, plus an in-place aliasing bug in `f6mulv` found by the same
  harness.  Ten kernels are now verified against a plain Python bigint tower on
  10000 random vectors each, with a per-kernel mutation that must fail
  (`all_pass=True`, `mutations=True`).  The corrected `cyclotomic_sqr` span is
  25372 Ir and is the squaring price; the corrected `f12mul` span (47924 Ir)
  prices the verification.
* **Protocol.** The meter is created before the table.  Table, code, every start,
  every step (guard-aware), lookup, record, replay and verification are charged and
  fits-checked; a replay that does not fit abandons its worker; verification is a
  metered square-and-multiply, not an uncharged `ok()`.  The 32-squaring abort is
  implemented and its delay measured in a forced-admission mutant (E[delay]=31.26923076923077,
  max 32); the abort-disabled control never abandons.
* **Census.** 512 fresh fixtures x one trial per group, two groups:
  292/512 (CP-95 0.533116) at rp=521041, 304/512 (CP-95 0.556735) at rp=10566001.  Verify
  and replay rows are non-zero in both groups, and the census budget executes the
  claimed S/sqrt(N) steps.
* **Record/lookup.** Serialized record length is 688 B (48 B key + 2x32 B labels
  + 576 B y), a machine row, not the v24 64 B claim.  `lookup_count.c` ships and its
  disassembled 95 Ir span is the lookup row.
* **Ties.** The v24 offset-doubling explanation is deleted; the offset stays `j`.  Ties
  are measured (2983/200000) and resolve deterministically (first minimum;
  sorted-key pool key).
* **Claim policy.** `floor(log2(W/C)*1e4)/1e4`; no +0.0031 fudge.  `sum(rows)==W`
  exactly, with a planted typo control.

    ## Objective

    Repair the five concrete defects in verdict d3d07cd8 for the fixed-Q GT
    discrete-log track and submit the corrected bound.  The target is the weakest
    sound bound at the census shape S/sqrt(N)=1.45 with one record per visited class.

    ## Environment and setup

    macOS arm64 (Apple clang 21, otool), Python 3.9 for the census/survival harness and
    3.14 for the bigint sweep.  No external dependencies: `dependencies.json` is `[]`.
    The priced kernels are compiled with `clang -O2 -target arm64-apple-macos11`; the
    verification harness compiles the same C at -O0 so the unrolled -O2 static spans do
    not dominate compile time (semantics identical; the priced object is built
    separately at -O2 by count_ir.py).

    ## Prior work / baseline

    v24 (submission d3d07cd8) claimed 126.5009 and was rejected for: (1) census control
    flow not the declared protocol (replay did not abandon, verification never metered,
    guard conversions priced as multiplies, table before meter); (2) the kappa-offset
    offset-doubling explanation was algebraically wrong; (3) the 32-squaring abort existed only on
    paper and delay=1 was substituted; (4) the record was claimed 64 B though the
    serialized key+labels is at least 112 B, lookup_count.c was absent, and the row sum
    was not the declared bound; (5) kernels.c f6mul swapped the Fp6 cross terms and the
    square report showed 10000/10000 mismatches.

    ## Implementation and files changed

    * scripts/ir25/kernels.c: f6mul cross terms repaired (u0 -> c0, u2 -> c1, u1 -> c2);
      f6mulv made alias-safe (the v25 sweep found it); f2inv/f6inv/f12inv/f12mul/f12sqr
      added; cyclotomic_sqr kept (Granger-Scott) and now correct through f6mul.
    * scripts/ir25/verify25.py: plain-bigint tower reference, 10 kernels x 10000 random
      vectors, one targeted mutation per kernel that must fail.
    * scripts/ir25/lookup_count.c: the hash + two-probe + six-word-compare lookup symbol
      the ledger cites.
    * scripts/ir25/count_ir.py: otool -tvV spans of the corrected -O2 symbols; the only
      source of kernel prices.
    * scripts/protocol25.py: meter created before the table; table/code/start/step/
      lookup/record/replay/verify all fits-checked; guard-aware step price; replay and
      verification abandon on failure; 32-squaring abort implemented; forced-admission
      mutant and abort-disabled control for the delay measurement.
    * scripts/gen25.py: ledger, certificate, argument, note and submission; exact
      row-sum assertion, planted typo control, floor claim.

    ## Exact commands

        python3 scripts/ir25/count_ir.py
        NVEC=10000 python3 scripts/ir25/verify25.py
        python3 scripts/protocol25.py '{"mode":"census", ... ,"fixtures":512}'
        python3 scripts/protocol25.py '{"mode":"scan", "steps":120000, "segs":48}'
        python3 scripts/protocol25.py '{"mode":"scan", "force_admit":true, ...}'
        python3 scripts/protocol25.py '{"mode":"scan", "force_admit":true, "abort_enabled":false, ...}'
        python3 scripts/gen25.py
        yukon submit --track bls12-381-fixed-g2-invert-g1 --note-file .../public-note.md

    ## Experiment configurations

    Census: groups rp=521041 (q=523) and rp=10566001 (q=547), Rg=65536, depth=5, dpbits=0,
    8 workers, ratio=1.45, 512 fresh fixtures x one trial, one process pool of 24.
    Scan: Rg=65536, depth=5, q=523, 48 segments x 120000 steps (5,760,000 steps).
    Forced mutant: 8 segments x 40000 steps, an admission injected every 1000 steps;
    control identical with abort_enabled=false.  Survival: 400 pairs, Rg=65536, gaps to
    1000, three arms.

    ## Measured results (including unsuccessful runs)

    * Kernel sweep: all 10 kernels 10000/10000 exact; every mutation 10000/10000 failing.
      The first sweep attempt with 20 vectors exposed the f6mulv aliasing bug (f12mul
      20/20 mismatches) before the 10000-vector run; the v24 square failure is the
      baseline 10000/10000 mismatch on cyclotomic_sqr.
    * Census: 292/512 (CP-95 0.533116) at rp=521041 and 304/512 (CP-95 0.556735) at
      rp=10566001; 292 and 304 metered verifications; replay rows non-zero; 0 abandons,
      0 aborts.  The census took 370.6 s.  A 4-fixture smoke census read 0/4 before two
      fixture-scale artefacts were fixed and disclosed: the replay charge must be the
      O(1) record resolution (not a trail replay), and worker aborts must not end the
      fixture while other workers still have budget.
    * Forced mutant: 312 admissions, 312 aborts, E[delay]=31.269 (max 32); the
      abort-disabled control never abandons (9744 stalls past 32, 0 abandons) and the
      shipped scan has 0 admissions in 5,760,000 steps.
    * Survival: 0.995 at gap 1, 0.4875 at gap 1000, fitted rho=0.999323, band 1.358.
    * Ties: 2983/200000 = 1.49e-2, resolved deterministically (first minimum; sorted-key
      pool key).  The task's expectation was ~2^-w; the measurement says otherwise, so
      the measured value is shipped.
    * Claim: W=2403108794320843297282781564631794087629139 Ir, STEP=24300.277,
      claim=126.4981 = floor(log2(W/20000)*1e4)/1e4.

    ## Limitations

    The census is at two small groups, not production rp; the table build is charged at
    its fixture proportional share because the true cost is 2^99x the fixture budget; the
    fixture worker count is 8 vs 1024, so collision processing is funded separately and
    the production replay row is the exact birthday expectation at 1024 workers.  Memory
    is stated, not scored.  The survival decay has no derived mechanism in this version:
    the v24 offset explanation is withdrawn and the curve is used only at the measured
    gap-1 point.

    ## Next steps

    Measure a third group to widen the scaling ladder; add a counted Frobenius symbol so
    the replay row is disassembly-priced end to end; test the abort under a synthetic
    production-scale trap generator.

    ## Provenance

    All evidence files are regenerated by the shipped scripts; the v24 tree remains in
    the archive for diffing.  The model and harness are recorded in the submission
    metadata.  No private reasoning is included.
    