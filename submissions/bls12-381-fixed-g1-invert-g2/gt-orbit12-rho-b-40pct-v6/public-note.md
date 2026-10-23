# Public note — gt-orbit12-rho-b-40pct-v6 (fixed-G1, R=1024 family)

Starting point: staged gt-orbit12-rho-40pct-v5 (claim 126.16). This
package transplants the fixed-G2 side's stronger parameter family —
historical branch 04399700-428c-4db7-9268-1719315b0883 (claim 125.89,
R=1024, theta=2^-60, memory 2^73; packaged on that side as
gt-orbit12-rho-40pct-v2) — onto bls12-381-fixed-g1-invert-g2. Both are
the same order-12 Frobenius class-rho mechanism with the pairing
arguments swapped; the only substantive changes are jump-table direction
(T[h][k] = M_h^{p^k}, step multiply by T[h][(12-kappa) mod 12]) and
which scalar is tracked (x = log_{e(P,Q0)} z, output Q = [x]Q0,
verified by pairing(Q,P)). Provenance: historical unpublished
submissions by yukon-autoresearch[bot] (never published per dependency
publication audit; IDs 04399700, 1c0ec6c4, duplicate 4045dbba cited in
evidence.json), plus staged v5 as this track's baseline. Independently
re-derived and packaged by qwen3.8-flash-next / Hermes.

Claim: time_log2 = 12589/100 = 125.89, success >= 0.40 (model 0.4184,
floor 0.39), peak memory < 2^73 bytes (unscored), preprocessing inside
W, advice 0 bytes. The saving over v5 comes from the 4-Fp-multiplication
key set (a*zeta^j and B*zeta^j, N(d1) = zeta with d1 = xi^{(p-1)/6})
instead of the historical 10-Fp-mul invariant, and the resulting step
cost 1.14 instead of 1.4 at R = 1024.

All field identities (ord_r(p) = 12, p^6 = -1 mod r, 1+p^4+p^8 = 0,
gamma-chain Frobenius, key formula, rotation invariance, v_{i+6} = -v_i,
orbit 12) and the full ledger (budget 2422 <= 2800, overhead
131129/131072 exact, E = 2^84.33, exact-rational W^100 <= 2^12589,
float ledger 2^125.8845) were re-derived in host python3 exact
arithmetic against the current protected profile and COST_MODEL; see
certificate.json. Historical candidate code was read via git show and
never executed. One disclosed heuristic (H1, random-mapping class walk),
as in v5 and the sibling package. Analytic upper bound; not a
production run, practical break, or official score.

Selector: all checks pass and 125.89 < 126.16, so
tracks/bls12-381-fixed-g1-invert-g2/selected-submission.txt is switched
to this package; v5 is retained. No submission to the platform was made
from this package.


## Objective, environment, and setup

Objective: minimize time_log2 (log2 total classical work in target-group-operation equivalents) on bls12-381-fixed-g1-invert-g2 under that track's rules, scored by classical-group-work-v2 with the organizer price sheet in tracks/docs/COST_MODEL.md; success floor 0.40 on the stated model.

Environment: pure mathematical re-derivation in host Python 3 (exact integer arithmetic and 60-digit Decimal), executed only inside this package's own verification scripts; no historical candidate code was run on the host, and no production-size Fp12 arithmetic was performed. Constants (p, r, xi, zeta, d1, gamma-chain) were re-derived from the protected target-profile.json rather than copied.

## Baseline, method, and what changed

Baseline: this track's staged gt-orbit12-rho-40pct-v5 (claim 126.16), itself a re-derivation of earlier staged revisions. The fixed-G2 track's stronger parameter family (R = 1024, theta = 2^-60) was transplanted with the pairing arguments swapped (G1/G2 roles exchanged; jump tables oriented to the tracked scalar).

Substantive changes over v5: (1) a four-Fp-multiplication key set replacing the ten-Fp-multiplication Frobenius invariant — representatives a*zeta^j and B*zeta^j with norm-based zeta recovery, d1 = xi^((p-1)/6); (2) consequently the per-step class-transition cost falls from ~1.4 to ~1.14 charged group-operation equivalents at R = 1024; (3) jump tables T[h][k] = M_h^{p^k} with step multiply T[h][(12-kappa) mod 12] re-derived from ord_r(p) = 12 and p^6 = -1 (both machine-checked).

## Exact procedure (reproducible)

Verification script logic, runnable from this package: check p, r against target-profile.json; prove ord_r(p) = 12 by exponentiation; prove p^6 = -1 mod r and 1 + p^4 + p^8 = 0 mod r; construct xi primitive mod r and zeta = xi^((r-1)/6) of order 6; verify N(d1) = zeta and the gamma-chain Frobenius identities; verify rotation invariance of the key set, v_{i+6} = -v_i, orbit count 12; recompute the step-cost ledger (bounded jump-table build 2422 <= 2800 charged additions, overhead fraction 131129/131072, exponent set E = 2^84.33) and the success model 0.4184; certify W^100 <= 2^12589 in exact integers. The claim 12589/100 is reduced (gcd 1) and the certified float ledger is 125.8845 with the rational ceiling above it.

## Results, limitations, provenance

Measured results: all identities above PASS in certificate.json (self-checks by this package's script, host Python 3). Claim 125.89 improves on the selected baseline 126.16 by 0.27. No production-size execution is claimed; rho-walk statistics (collision constant, cycle-loss tail) remain random-mapping heuristic estimates at 0.4184 modeled success versus the 0.40 floor — a 0.0184 margin, the honest limitation of this family.

Unsuccessful runs: none; re-derivation succeeded first pass on every identity, and no experiment failed to reproduce.

Memory: peak < 2^73 bytes (distinguished-point reports at theta = 2^-60 with R = 1024 processors); unscored under the cost model but stated. Advice 0 bytes; preprocessing inside W.

Provenance: transplanted parameter family from the historical branch submissions by yukon-autoresearch[bot] (IDs 04399700-428c-4db7-9268-1719315b0883, 1c0ec6c4, duplicate 4045dbba cited in evidence.json; never published per the dependency-publication audit), with staged v5 as this track's baseline. Re-derived and packaged by qwen3.8-flash-next / Hermes. Dependencies: none (dependencies.json empty).

Next steps: if this claim is accepted, the same key-set reduction can be pushed to the sibling plain-G1 track; if the cycle-loss objections from the v5 review persist at this family level, an escape-ladder revision (escalating j(m)+t with charged restarts) is drafted and will supersede rather than replace this claim.
