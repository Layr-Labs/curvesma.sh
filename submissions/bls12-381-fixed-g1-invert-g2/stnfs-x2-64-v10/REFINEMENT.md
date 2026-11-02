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

## Slack, charge, UDUP, thin and plan updates (this lineage)

Slacks (295,224)->(292,223) at theta 19/50, X2 2^64 (rels 2^64.13021, margin
0.3003, allowance 1.0671, collection 2^121.1444 with route 10.01). LA pad
1.10->1.096 (N/U 0.05163), C_ENTRY 21->20.7, C_ROW 2750->2710: LA 2^117.7103,
generator 2^107.1997. UDUP 0.25->0.5 (36% below 0.779; [11] 12 points at
0.944). THIN_MAX 90->85 (bin 90 keeps searched flag; cap 2^90 covers Qq'),
MIN_RELS 16->18 (boot failure 2.39e-04). Joint-move re-search (mu>=18,
phi<=0.05, psi<=0.8; 35 entries differ from base, boot unchanged;
refinement.json): two-log 116.3024->115.8431, meter 2^116.5664. Same H1-H5;
analytic, not a full-size run.

## Both-bars alpha (5.15, 0.85), theta 19/50, slacks (295, 224) (this revision)

Both-bars alpha (5.15, 0.85): f-side full 50k sweep 5.085035 (313 s; only
p = 19 has repeated f roots below 20k), charge 5.15 with 0.065 margin (clears
the 50k bar that rejected 6802beb); g-side 0.731543 with primitive-correct
p = 2 term, charge 0.85 with 0.118 margin ([14] checks both 50k values).
Descent alpha_f 5.0->5.15 (two-log 115.8581, cap 116.58 refit). Theta
759/2000->19/50 (1040 rows, margin 0.3014, allowance 1.0674), slacks (295,224)
kept: collection 2^121.1553. Inputs 64.131, 746.0/498.9; plan unchanged.
Analytic, not a full-size run. Supersedes v7 (3a00ffc8) under the new f-side bar.

## UDUP 0.6, descent alpha_g 0.85, slacks (294, 225) (this revision)

UDUP 0.5->0.6 in the descent model (mu, boot trials, retry chain), supported
by a 200-point unit-mate convergence study (f-side 260-bit special-q lattice
seed 13300, R = 96, seed 13307): distinct/raw 0.8514/0.8394/0.8353 (SE
0.023) at 500/840/1190 units, flattening; 0.6 is 23% below the 480-point
0.779 and 10 SE below the 1190-unit value. Descent alpha_g 1.0->0.85,
matching the collection charge with the same 50k support ([14] checks both
50k values on both sides). Slacks (295,224)->(294,225) (margin 0.3012,
collection 2^121.1547, allowance 1.0674 kept); inputs unchanged (64.131,
746.0/498.9), plan unchanged. Two-log 115.8581->115.5527, cap 116.58->116.28
refit for success 0.3943, boot failure 2.63e-04->4.80e-05, E[N]
2^17.63->2^17.26. Total 2^121.32592. Every manifest/proof bound realigned
with the audited model (meter, Markov, nodes, retention, ECM schedules).
Analytic, not a full-size run.

## 100k pool, theta 2901/8000, slacks (295, 226) (this revision)

Pool 20k->100k draws: seed 2 extended from 18000 to 98000 indices with the
same sampler (`v11-production-ideal-key`, theta = 1/2), every canonical row
kept (11908 rows, no cherry-picking); all 20000 old indices regenerate
bit-identical. The pool file stores delta-varint canonical index lists plus a
200-row field-anchor table with deterministic field regeneration (`load_pool`);
[2] checks membership both ways plus every anchor field, and all 100,000
indices were re-derived offline with an independent code path (exact match).
Theta 19/50->2901/8000 (4518 rows, a-radius 872, certified products
0.971219877/1.027872797), margin floor 0.30->0.12 (margin 0.1201 at 2.09
SE, pool SE 0.057), allowance 1.0674->1.0318 (work SE 0.01590),
slacks (294,225)->(295,226): collection 2^120.9506. Inputs 63.949,
746.9/499.8; plan unchanged (min mu 33.8, max psi 0.644 boot, max phi
0.034). Two-log 115.5527->115.8579, cap 116.28->116.59 refit for success
0.3953, boot failure 4.80e-05->3.60e-04, E[N] 2^17.26->2^17.74.
Total 2^121.15796. Every manifest/proof bound realigned with the audited
model (meter, Markov, nodes, retention, ECM schedules, lattice products,
pool header). Analytic, not a full-size run.
