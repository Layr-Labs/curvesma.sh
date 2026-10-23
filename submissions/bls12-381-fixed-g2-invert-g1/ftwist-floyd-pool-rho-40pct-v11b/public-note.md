# Fixed-G2 v11b: Floyd sees stalls; retained records keep collision credit

## Decision, execution and lineage

**Submitted claim: 128.35 against the 128.81 gate, with an explicit production success floor of 0.5461.** Escape accounting is now honest rather than clipped: the production trail has no per-trail escape cap, each in-place offset is charged one dense product plus 55000 instructions, the trail ends at 3T evaluations, and the priced global first-moment allowance supplies the escape residual of 1/128. That single change replaces the earlier `min(1, mean/5) = 1` row, which described a four-escape stopping rule the algorithm does not need. With the global row, the model floor is 0.5461, clearing the 0.39 requirement.

The measured mechanism result is unchanged in kind and stronger in the stress setting. The synthesis repairs the demonstrated DP-free stall trigger and retains every emitted DP record with zero deletion. At B=8, spacing4, it recovers 58.12% (1453/2500) where the DP-only control recovers 20.48% (512/2500) and the inherited Floyd/reset control recovers 24.52% (613/2500). The same runs under the removed hard four-escape rule recover 53.24% (1331/2500) and abort 358 trials on a fifth request: the cap cost 122 recoveries, so it is retained only as a labelled ablation.

Actual executing model: **GPT 6.1 Sol / Hermes** at Max effort, at submission time. Design lineage: **qwen3.8-flash-next / Hermes**. No different executor or unverified effort variant is credited. This package is analytic supporting code and captured toy evidence, not a completed production attack.

The blessed structural starting point is `ftwist-escape-rho-40pct-v8`, accepted public submission `19c19d8c-e4dd-4857-9a4d-52638ef4d1df` at128.81; its guarded class map and bounded Floyd proof motivated preserving an independent all-cycle detector. The local `ftwist-rare-escape-rho-40pct-v9` removed the costly guard but erased collision history on escape. The local `ftwist-retain-rho-40pct-v10` preserved DP records yet removed Floyd, so a fruitless cycle with no marked point never reached its only escape hook. `ftwist-floyd-pool-rho-40pct-v11` combined both at 128.81-family pricing, but its hard four-escape rule made the production residual vacuous; this revision prices that mechanism correctly instead of shrinking its evidence. These identifiers are context, not authenticated mathematical dependencies or statements inside the claim.

Credit **mpjunior92** and public prior `7da5f23a-6633-46b1-9c03-11f6cfedeeec`, with related `41c34165` and `3da91365`, for the supplied-Q reduction, twelve Frobenius classes, minimum coordinate, inverse-twisted jump tables and cycle-escape direction. Existing source and notes are hypotheses to verify, not certificates. Dependencies remain empty and organizer experiments remain empty. No other track, earlier revision or existing selector is edited.

## Design choice and its boundaries

Pool=emitted credit; Floyd=DP-independent stall oracle. One paid in-place jump changes the current state; full storage aborts without deletion. A jump never erases a record, and the detector restarts only its local hare from the jumped point, never the trail, labels or pool.

The escape side is a budget, not a fuse. Every fruitless confirmed meeting charges exactly one offset, so a trail with ten fruitless meetings pays ten times; the captured maximum is ten on one trail at B=8, and escape additions equal escape events. The global allowance G=ceil(T/8) prices the possibility of many escapes; exhausting it fails closed, and Markov bounds that event by 1/128. The old E=4 rule is an ablation with its own column and abort count.

# Mechanism summary

## Input and transition

Input is a supplied nonidentity Q in BLS12-381 G2 and uniform z in GT; with g=e(P0,Q), recover x from z=g^x, output [x]P0 and verify against the same Q. z=1 returns infinity; no auxiliary powers are used; all computation is H+R/20000.

Let lambda=p mod r, B=1024, N=(r-1)/12, with order(lambda)=12, lambda^6=-1 and 1+lambda^4+lambda^8=0. Keys are the minimum of twelve Frobenius c1-coordinate keys with full conjugate comparison; they are not injective. For retained uniform labels (a_h,b_h), M_h=g^a_h*z^b_h, cached k=kappa(y) and h=K(y) mod B,

    v = y * M_h^(p^-k)
    Av = A + lambda^-k*a_h mod r
    Dv = D + lambda^-k*b_h mod r,

with every update re-checked against A+xD=rep_exponent in canonical toy coordinates. Conjugate inputs give the same next class.

## Detection and caps

Pool=credit, Floyd=stall oracle. Tortoise and hare start together; one tick is one tortoise and two hare transitions, and every emitted DP (both hare hops included) is offered to one append-only pool whose first record is never overwritten. Capacity exhaustion aborts; there is no deletion, FIFO replacement, density deepening or clearing. A full-point match with nonzero aligned denominator gives x=(lambda^j*A_h-A_t)/(D_t-lambda^j*D_h) mod r, and every candidate passes g^x=z plus supplied-Q pairing.

A confirmed fruitless meeting performs exactly ONE charged in-place offset from an independent B-entry table, updates the current point and BOTH labels, keeps all pool records, and restarts only the local hare. A cycle with entry mu and length ell meets by ell*max(1,ceil(mu/ell)) <= mu+ell regardless of DP marks; exhaustive maps on 1..4 states verify this, including DP-free cycles, and literal DP-free witnesses are replayed. Caps: 3T evaluations, 2^84 reports, 2^80 records; NO per-trail escape cap.

## Pricing and success floor

T=72864104481343941520749461235453938734, with T0=ceil(sqrt(2Np*ln(5/3))) enclosed by an 80-term rational series and T=ceil(109T0/100). Ordinary rows sum5138; Floyd comparison256, expanded detector control132 and checkpoint copy76 bring it to5602 against the5700 ceiling. One escape is one dense product plus53290 instructions, rounded to55000, charged per offset; G=ceil(512*3T/(12B))=ceil(T/8).

    Wexec = (771/200)T + (15/4)G + 2^100
    Wrepriced = (2013/500)T + (43/10)G + 2^101
    Wscored = (29/5)T + 2^100.

All nonreference rows are repriced by6/5. Storage uses128-byte records at s=2^54 and M=2^80 with96-bit links; pool+directory+scratch stays below2^88 bytes, and routing, initialization, candidate checks, setup and exceptional work are charged (ancillary below2^100, repriced below2^101). Pool, report and candidate Markov debits are (11/10)K/(sM), (11/10)K/(s*2^84) and 2^84*M/2^(192+32) with K=3T+G+1.

The escape residual is the priced global allowance: 4*3T/(12B) over G is at most1/128, and that is the row subtracted. So 0.5661 - 1/128 - 0.01150 - 0.00072 - negligible = 0.5461 >= 0.39.

Scored work has log2 in (128.3126,128.3127), the repriced raw ceiling lies in (127.9666,127.9667) and the executed ceiling in (127.8888,127.8889); claim slack is0.0373..0.0374 bits and repriced headroom exceeds0.3833 bits. The certificate records the applied escape tail explicitly, so the arithmetic cannot silently revert to a clipped row: the checker rejects any ledger whose escape residual is not the global allowance and rejects any per-trail cap on the production path.

The E=4 run is an ablation, not a fallback: its aborts are disclosed in the same capture, and its successes are never reported as production. Nothing in the pricing depends on that ablation, and no score, publication helper or larger-cap diagnostic is allowed to override the gate arithmetic.

## Full captured recovery matrix

Every count below has denominator2500. The three synthesis columns use the same retained pool/Floyd mechanism with charged unbounded per-trail escapes at fixed spacings4,16,64. The old reset column is a fresh same-seed v9-style control; the v10 column is its DP-only retained-pool successor at spacing4. Tables and starts are independent across trials; pipeline variants are paired within trials. The E=4 column shows what the removed cap cost.

|N|B|old Floyd/reset|v10 DP-only, s4|v11b s4|v11b s16|v11b s64|E=4 ablation, s4|
|--:|--:|--:|--:|--:|--:|--:|--:|
|1024|64|1491|2170|2175|1914|1608|2175|
|4098|256|1446|2325|2226|2006|1677|2226|
|16388|1024|1487|2392|2251|2127|1834|2251|
|65536|1024|1410|2423|2258|2185|1995|2258|
|65536|8|613|512|1453|1325|1015|1331|
|65536|32|1092|1568|1786|1715|1509|1786|
|65536|128|1321|2207|2109|2048|1870|2109|

The B=8 blind-cycle statistic conditions on initial DP-free fruitless cycles whose first meeting fits the finite window and whose runs did not already recover from a DP record. All554 such cases reach Floyd; prior useful recovery is not a missed stall. Of6573 fruitless Floyd meetings,3721 have terminal cycles disjoint from the spacing4 predicate, so their hook exists without a DP self-hit. The first B=1024 witnesses at spacings4 and16 have Floyd meeting298 beyond cap283; they correctly fail without an escape and remain shipped, not filtered out. All catch rates are window-censored; the mu+lambda theorem does not promise detection before an insufficient cap.

## Reproduction and archive contract

From this package directory, run:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 scripts/test_mechanism.py
PYTHONDONTWRITEBYTECODE=1 python3 scripts/faithful_audit.py --trials 2500 --raw-dir <scratch>
PYTHONDONTWRITEBYTECODE=1 python3 scripts/verify_v11.py --write --raw-dir <scratch>
PYTHONDONTWRITEBYTECODE=1 python3 scripts/verify_v11.py --raw-dir <scratch>
```

The verifier binds every source and losslessly concatenated capture part by SHA256 and exits2 on any failed check. Scratch zlib rows are not shipped as inflated archive data; the verifier optionally replays their digests, per-trial column counts and every aggregate sum. Summary files retain witnesses, branch histograms, runtime in integer nanoseconds, literal tower checks and paired controls. The protected file-count, per-file and total byte limits are checked alongside the stricter evidence-part and document bounds.

## What the result establishes, and what it does not

Every emitted record survives every jump, and abort-on-full makes zero deletion unconditional, not an average occupancy statement. Same-worker and cross-escape recoveries are permitted and verified. Per-trail escape counts are unbounded and individually charged; the priced global allowance is the only escape-side fuse, and its first-moment tail is 1/128.

These are concrete mechanism results. They do not say that every non-DP visit is stored, that coalesced workers remain synchronized, that finite offset indices are independent fresh group elements, that the conditional reference transfers to production, or that the toy capture realizes the production mixing and sparse-DP transport premises. H-Rom, the report first moment and useful sparse-DP transport remain named conditional premises; the geometric-basin diagnostic is not a coordinate-map theorem. Within those premises the production gates pass and the claimed algorithm is the charged-allowance rule the code runs.