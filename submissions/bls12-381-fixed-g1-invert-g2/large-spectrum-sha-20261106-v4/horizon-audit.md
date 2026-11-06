# Exact finite public-parameter audit

All numbers here are public profile parameters or algorithm design constants.
No sampling secret, measured production success or instruction trace is used.
The following is a check specification, not a declared executable experiment;
review execution, if any, belongs in the shared verifier's isolated runner.

Use integer r from protected target-profile.json and exact rational arithmetic.
Define R=isqrt(r-1)+1, L=(1031*R+999)//1000, s=(L+61)//62, K=L+s,
T=271, eta=1/65536, D=(L-T+1)*(L-T+2)//2. Then define

m=(1-eta)*D/r

U=(1+3*eta)*D^2/r^2+(1+eta)*D/r
  +((131072+131070*eta)/1073676289)*(L+1)^2/r
  +8*T*eta*(L+1)^2/r^2+12*eta*(L+1)/r.

These are exactly the finite moment bounds in collision-proof.md; the horizon and spectral coefficient are revised as proved in that file. The following five inequalities all hold as exact
rational or integer comparisons:

1. (3*m-U)/2 > 39015/100000.
2. 39015/100000 - 1/8192 - 1/2^100 - 1/2^980 - 1/r > 39/100.
3. (K*34720+2^90)^250 < 20000^250 * 2^32073.
4. K<r and K+1<2^129.
5. 4*r*3^256*2^100 < 8^256.

The first and second imply success strictly greater than the declared0.39
under the sole disclosed SHA random-function heuristic. The third implies
the full cold computation bound2^(32073/250), with C fixed at20000. The
fourth supports coefficient nonsingularity and the inherited compression-loss
bound. The fifth is the inherited explicit bad-table tail check. All losses
and failed outcomes are charged. Rational score32073/250 equals128.292
exactly; it is an upper bound on log2(total work), not a rounded measured score.

For orientation only, the Bonferroni lower bound before losses is about 0.39015717, and after compression loss 1/8192 and other failures it remains above 0.390035. The complete-work logarithm is about 128.2914396. Exact inequalities above establish the claims. The jump-table spectrum, horizon, scalar costs and memory are revised as explained in proof.md; no production experiment is claimed.
