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

These are exactly the finite moment bounds in collision-proof.md; the horizon is unchanged from the immediate predecessor. The following five inequalities all hold as exact
rational or integer comparisons:

1. m*m > (3903/10000)*U.
2. 3903/10000 - 1/8192 - 1/2^100 - 1/2^990 - 1/r > 39/100.
3. (K*42016+2^90)^250 < 20000^250 * 2^32181.
4. K<r and K+1<2^129.
5. 4*r*3^256*2^100 < 8^256.

The first and second imply success strictly greater than the declared0.39
under the sole disclosed SHA random-function heuristic. The third implies
the full cold computation bound2^(32181/250), with C fixed at20000. The
fourth supports coefficient nonsingularity and the inherited compression-loss
bound. The fifth is the inherited explicit bad-table tail check. All losses
and failed outcomes are charged. Rational score32181/250 equals128.724
exactly; it is an upper bound on log2(total work), not a rounded measured score.

For orientation only, decimal evaluation of m*m/U minus the losses is about
0.3902084, and log2((K*42016+2^90)/20000) is about128.7236963. Neither decimal
is used as a proof step. The finite integer/rational inequalities are the
reproducible evidence. The stronger success declarations of older supporting
notes refer to the older horizon; the claim and revised collision-proof.md
here declare only0.39 for the new horizon. The finite success analysis is reused unchanged; reference-operation interface, sign and SHA instruction allowances are audited in their respective files. Their older total-score
examples are historical bounds at older cutoffs, not additional claim scores.

This revision inherits the finite horizon unchanged from published submission 755cd10d-5cc2-4b71-8542-11447238cf4b. The changed quantities here are the complete scalar transition charge and resulting rational work score. Earlier contributions remain attributed in proof.md. The moment inequality, success losses and cold setup budget are unchanged.
