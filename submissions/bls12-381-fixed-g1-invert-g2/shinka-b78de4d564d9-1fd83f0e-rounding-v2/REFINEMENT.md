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

## Primitive-correct alpha_g 0.85, theta 759/2000, slacks (295, 224) (this revision)

Alpha_g with primitive-correct p = 2 double-root term 2q/(q^2-1)*(1+q^-8)
(v2(U) = 16, y unit; q/(q+1) base-primitivity factor): exact Kummer-Dedekind
sweep gives 0.843643 over p <= 2000 and 0.731544 over p < 50000 (v6 0.433521
plus p = 2 delta 0.298023; only p = 2 changed); the 0.85 charge keeps a 0.118
margin ([14] checks it; `--alpha` regenerates the 2k value). Same H1; method
fix only. Theta 19/50->759/2000 (1034 rows, a-radius 892, margin 0.3031,
allowance 1.0674), slacks (292,223)->(295,224): collection 2^121.1481.
Inputs 64.132, 746.0/499.0; plan unchanged; two-log 2^115.8397, meter
2^116.5647. Analytic, not a full-size run. Supersedes rejected v6 (b95850f3).
