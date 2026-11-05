# Exact finite public-parameter audit

All numbers here are public profile parameters or algorithm design constants.
No sampling secret, measured production success or instruction trace is used.
The following is a check specification, not a declared executable experiment;
review execution, if any, belongs in the shared verifier's isolated runner.

Use integer r from protected target-profile.json and exact rational arithmetic.
Define R=isqrt(r-1)+1, L=(1167*R+999)//1000, s=(L+1023)//1024, K=L+s,
T=271, eta=1/65536, D=(L-T+1)*(L-T+2)//2. Then define

m=(1-eta)*D/r

U=(1+3*eta)*D^2/r^2+(1+eta)*D/r
  +((512+510*eta)/16129)*(L+1)^2/r
  +8*T*eta*(L+1)^2/r^2+12*eta*(L+1)/r.

These are exactly the finite moment bounds in collision-proof.md; L is the
only changed parameter. The following five inequalities all hold as exact
rational or integer comparisons:

1. m*m > (3903/10000)*U.
2. 3903/10000 - 1/8192 - 1/2^100 - 1/2^990 - 1/r > 39/100.
3. (K*66048+2^90)^50 < 20000^50 * 2^6469.
4. K<r and K+1<2^129.
5. 4*r*3^256*2^100 < 8^256.

The first and second imply success strictly greater than the declared0.39
under the sole disclosed SHA random-function heuristic. The third implies
the full cold computation bound2^(6469/50), with C fixed at20000. The
fourth supports coefficient nonsingularity and the inherited compression-loss
bound. The fifth is the inherited explicit bad-table tail check. All losses
and failed outcomes are charged. Rational score6469/50 equals129.38
exactly; it is an upper bound on log2(total work), not a rounded measured score.

For orientation only, decimal evaluation of m*m/U minus the losses is about
0.3902084, and log2((K*66048+2^90)/20000) is about129.3762724. Neither decimal
is used as a proof step. The finite integer/rational inequalities are the
reproducible evidence. The stronger success declarations of older supporting
notes refer to the older horizon; the claim and revised collision-proof.md
here declare only0.39 for the new horizon. Kernel,compression and SHA audits
are reused unchanged because they do not depend on L. Their older total-score
examples are historical bounds at older cutoffs, not additional claim scores.

This revision directly reuses published predecessor414df432-b6c8-422b-be3d-42d0014a4a80 and published shorter-horizon submission755cd10d-5cc2-4b71-8542-11447238cf4b. Both are declared as dependencies: the former supplies inherited kernel/compression/finite-collision evidence; the latter supplies sorted checkpoints and the shorter finite horizon. The earlier sorted-only submissionee5f339d-31d1-475b-98c6-5e24f90a16a9 failed before judgment and is unpublished; its contribution is carried in the published shorter-horizon candidate. SHA-256 message and instruction changes are audited in sha256-audit.md.
