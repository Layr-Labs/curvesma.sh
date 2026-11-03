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
generator 2^107.1997. UDUP 0.25->0.5 (unconditioned mate study, now
superseded). THIN_MAX 90->85 (bin 90 keeps searched flag; cap 2^90 covers Qq'),
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

## UDUP 0.6, descent alpha_g 0.85, slacks (294, 225) (superseded)

UDUP 0.5->0.6 (model; mate-study support superseded below), descent
alpha_g 1.0->0.85 ([14] checks both 50k values), slacks (295,224)->(294,225).
Two-log 115.8581->115.5527, cap 116.58->116.28. Total 2^121.32592.

## 100k pool, theta 2901/8000, slacks (295, 226) (superseded)

Pool 20k->100k draws (seed 2 to 98000; 11908 rows; 200 anchors;
re-derivation exact). Theta 19/50->2901/8000 (4518 rows, floor
0.30->0.12, margin 0.1201 at 2.09 SE), allowance 1.0674->1.0318,
slacks (294,225)->(295,226): collection 2^120.9506. Inputs 63.949,
746.9/499.8; plan unchanged. Two-log 115.5527->115.8579, cap
116.28->116.59. Total 2^121.15796.

## UDUP 0.65 (superseded support; model kept)

UDUP 0.6->0.65 (model; mate-study support superseded below).
Collection/LA/plan/inputs unchanged. Two-log 115.8579->115.6787, cap
116.59->116.41. Total 2^121.15081.

## 150k pool, theta 2881/8000 (this revision)

Pool 100k->150k draws (seed 2 to 148000; 17874 rows; 128 anchors; full
re-derivation exact). Theta 2901/8000->2881/8000 (6628 rows, rmin 869,
products 0.971121721/1.027970242), floor 0.12->0.11 (margin 0.1104 at 2.09
SE, pool SE 0.053), allowance 1.0318->1.0263 (work SE 0.01314), slacks
kept: collection 2^120.9254. Inputs 63.940, 746.1/499.9; plan unchanged
(min mu 36.2, max psi 0.614 boot). Two-log 115.6787->115.697, cap
116.41->116.43 refit for success 0.3957. Total 2^121.12966. Bounds
realigned (meter, Markov, nodes, products, pool header). Analytic.

## UDUP relation-weighted evidence (this revision)

v12c rejected the UDUP support (unconditioned E[1/m]; need E[s/m]).
Model/plan/numbers UNCHANGED (121.12966). New: yr = E[shat/m]/pr per
regime (nodes 5.25-145, min CI 4.64; boot 1.06-1.25, CI 1.02, binding).
Truncation x5/bound-8 closed. unit_mates triples; live [11] check; stats
table replaced. Anchors 128->64 (bytes). Analytic.
