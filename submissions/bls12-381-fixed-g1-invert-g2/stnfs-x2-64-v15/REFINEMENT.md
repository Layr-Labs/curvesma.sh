# Directed quadrature refinement and plan re-optimisation (inherited)

Inherited from e3940d11 and its source stnfs-lean-theta435-122017 (commit 852f57de5e67dd4a066e37e3f3d926c1e4bec592, submission 3a2e69ca); this package's configuration and numbers are in proof.md. Analytic exploratory claim, not a full-size run.

## Change

The descent integration mesh goes from 48 to 384 cells and the Dickman lookup spacing from 1/128 to 1/1024. PLAN16 and BOOT16 multiply iA by eight, so nA = 69 + iA (n - 69)/mesh is unchanged; TRUNC stays zero; retention uses the same global mesh. A coordinate search over every plan entry and the boot at the refined quadrature (mu >= 16, phi <= 0.05, psi <= 0.8) changed eight entries; refinement.json lists them and the recurrence before and after. The polynomial pair and physical prime windows are unchanged.

## Why the refinement is valid

For jA small-window and jB in {0, 1} large-window factors the model integrates rho((m - sum t + alpha)/69) prod dt/t / (jA! jB!). Each old mesh interval splits into eight; its logarithmic measure is exactly the sum of the new ones. Every subcell of an old wholly feasible cell stays feasible; cells crossing sum t = m are discarded; there is no truncation cutoff. On each subcell the left endpoint gives a lower value of rho (rho decreases), and moving left endpoints right can only raise it; the finer rho table moves the upper-rounded argument toward the actual one. So the refined integrals bound the same model, with no change to a smoothness heuristic. Only convolution rows through the requested kA are computed; a jA row depends only on rows 0..jA, so dropping rows above kA leaves every requested sum unchanged (the upper retention companion likewise). Four nonzero probes (kA = 0, 1, 4, 15) at the 1/1024 grid are recorded in refinement.json and recomputed by audit_refinement.py; they matched the full-row calculation exactly in binary64.

## Limitations

The rho coefficients come from 200-digit Decimal arithmetic and are evaluated in binary64; the table brackets the computed series. This is numerical evidence, not a formal interval certificate. The probability model and its production extrapolation stay heuristic (H1-H5). No favorable finite-density correction and no Samuels tail is used; Markov is kept.

## Provenance

Refinement and integration checks: the source authors (Codex, GPT-6 Astra/Luna). Candidate code here ran only in an isolated no-network container.

## Slack, charge, UDUP, thin and plan updates (superseded)

Slacks (295,224)->(292,223); LA pad 1.10->1.096, 20.7/entry, 2710/row;
UDUP 0.25->0.5; THIN_MAX 90->85, MIN_RELS 16->18; joint-move re-search
(35 entries). Two-log 116.3024->115.8431, meter 2^116.5664.

## Both-bars alpha (superseded inputs; charges kept)

Alpha (5.15, 0.85); descent alpha_f 5.0->5.15. Current inputs in proof.md.

## UDUP 0.6, descent alpha_g 0.85, slacks (294, 225) (superseded)

UDUP 0.5->0.6 (support superseded below), alpha_g 1.0->0.85,
slacks (295,224)->(294,225). Total 2^121.32592.

## 100k pool, theta 2901/8000, slacks (295, 226) (superseded)

Pool 20k->100k (11908 rows). Theta 19/50->2901/8000, allowance
1.0674->1.0318, slacks (294,225)->(295,226): collection 2^120.9506.
Total 2^121.15796.

## UDUP 0.65 (superseded support; model kept)

UDUP 0.6->0.65 (support superseded below). Total 2^121.15081.

## 150k pool, theta 2881/8000 (shipped config)

Pool 100k->150k (seed 2 to 148000; 17874 rows; re-derivation exact).
Theta 2901/8000->2881/8000 (6628 rows, rmin 869), floor 0.12->0.11
(margin 0.1104 at 2.09 SE), allowance 1.0318->1.0263, slacks kept:
collection 2^120.9254. Inputs 63.940, 746.1/499.9; plan unchanged.
Two-log 115.6787->115.697, cap 116.41->116.43. Total 2^121.12966.

## UDUP relation-weighted evidence (superseded sampler; kept design)

v12c rejected the UDUP support (unconditioned E[1/m]; need E[s/m]).
Model/plan/numbers UNCHANGED (121.12966). yr = E[shat/m]/pr per regime;
unit_mates triples; live [11] check; stats table. Anchors 128->64.

## Uniform-cell UDUP re-measurement (v13)

v12e rejected the study's sampler (round-and-reject not uniform over
cells). Model/plan/numbers UNCHANGED (121.12966). New: `int_ball_point`
(box rejection, exact uniformity by construction); 52-cell re-run (nodes
5.18-150, min CI 4.70; boot 1.10-1.25, CI 1.06, binding; all sextile pr); truncation
x5 (<=11.5%)/bound-8 closed; live [11] relw + shell-chi2 smoke check
(< 60); stats table replaced. Maxima studies keep the old sampler
(maxima+margins argument unchallenged). Analytic.

## Norm-sub-bin pr integration (v14)

pr over norm sub-bins (24+10 cells, uniform grids, margins kept, edges
pinned; boost capped 4x/1.5x); cap 115.635 (Markov held). Two-log
114.902. Total 2^121.10595. Plan/UDUP unchanged. Implied
integrated-pr floors: 4.70/4 = 1.175 nodes, 1.06/1.5 = 0.7067 boot.

## Descent plan re-optimisation (this revision, v15)

Coordinate search with joint moves at the integrated pr (128 improving
steps, mu>=18/phi<=0.05/psi<=0.8 kept, boot frozen); 45 entries differ
from base (compact diffs in refinement.json). Two-log 114.9016->113.6097,
cap 114.343 (Markov held). Total 2^121.08658. Claim 0 states the
integrated-pr floors 1.175/0.7067 (cap-division of the v13 floors).
