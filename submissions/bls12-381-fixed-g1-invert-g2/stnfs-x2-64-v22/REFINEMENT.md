# Directed quadrature refinement and plan re-optimisation (inherited)

Inherited from e3940d11 (source stnfs-lean-theta435-122017); configuration and numbers are in proof.md.

## Change

Mesh 48 to 384 cells, Dickman spacing 1/128 to 1/1024. PLAN16/BOOT16 multiply iA by eight, so nA is unchanged; TRUNC stays zero. A coordinate search at the refined quadrature changed eight entries (listed in the published e3940d11 record). Polynomials and physical windows unchanged.

## Why the refinement is valid

For jA small-window and jB in {0,1} factors the model integrates rho((m - sum t + alpha)/69) prod dt/t / (jA! jB!). Each old mesh interval splits into eight of equal total logarithmic measure. Subcells of wholly feasible cells stay feasible; cells crossing sum t = m are discarded; no truncation cutoff. Left endpoints give lower rho values (rho decreases), so the refined integrals bound the same model. Rows above kA are dropped (a jA row depends only on rows 0..jA). Four probes (kA = 0, 1, 4, 15) are recorded in refinement.json and recomputed by audit_refinement.py.

## Limitations

Rho coefficients come from 200-digit Decimal arithmetic, evaluated in binary64; the table brackets the series. Numerical evidence, not a formal interval certificate. The model stays heuristic (H1-H5).

## Provenance

Refinement and integration checks: the source authors.

## Later updates

Later changes: proof.md history; current numbers: proof.md section 8.

## Excess-1.1 plan re-optimisation (v2)

See DESCENT.md (method); the v2 plan record stays in the published 9eae858f package.

## Sub-bin pr, skewed boot, est-2SE input, re-optimised plan (v4)

Sub-bin pr (caps 4x/1.5x), skewed boot, est-2SE input 63.599, plan re-search (41 entries; `plan_v21`), cap at Markov 0.5, complete-orbit UDUP. Two-log 112.168, success 0.4968. Score 120.521203700786.

## Fresh on-region pool, core cut, 2-SE margin, re-measured filtering (v5)

Fresh 164,000-draw pool at theta 5/16 (seeds 1-2); qa < 0.02 rows give no yield (cells charged). Region 2347/8000: smallest n/8000 with margin >= max(0.06, 2 SE) (16247 rows, 2.01 SE). Filtering re-measured there (`collection_v22`); v4 plan re-evaluated at 63.589. Two-log 112.214, success 0.4967. Score 120.227673100157.

## Per-side factor-base bounds, explicit purge-loss charge (v6)

One bound per side, B_f = 2^69.25 and B_g = 2^69 (grid evaluated on the v5 pool, refinement.json `collection_v23`): U = 2^63.5788, region 2284/8000 (14987 rows, 2.01 SE), filtering re-measured on the per-side system (`collection_v23`); the >= 3-occurrence purge losses enter the missing-log proxy at the pooled toy rates `H3_RATE` ([6c]); the v4 plan is re-evaluated unchanged at LOG2_RELS = 63.716. Two-log 112.119, success 0.4968. Score 120.213628275286.
