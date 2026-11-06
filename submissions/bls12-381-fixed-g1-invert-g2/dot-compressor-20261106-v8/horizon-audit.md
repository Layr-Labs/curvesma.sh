# Exact finite public-parameter audit

All numbers here are public profile parameters or algorithm design constants.
No sampling secret, measured production success or instruction trace is used.
The following is a check specification, not a declared executable experiment;
review execution, if any, belongs in the shared verifier's isolated runner.

Use integer r from protected target-profile.json and exact rational arithmetic.
Define R=isqrt(r-1)+1, L=(25767*R+24999)//25000, s=(L+125)//126, K=L+s,
T=271, eta=1/65536, D=(L-T+1)*(L-T+2)//2. Define delta=K*(K+1)/2*(38/(2^61-1))^5 using exact rational division. Then define

m=(1-eta)*D/r

U=(1+3*eta)*D^2/r^2+(1+eta)*D/r
  +((131072+131070*eta)/1073676289)*(L+1)^2/r
  +8*T*eta*(L+1)^2/r^2+12*eta*(L+1)/r.

These are exactly the finite moment bounds in collision-proof.md; the horizon and spectral coefficient are revised as proved in that file. The following five inequalities all hold as exact
rational or integer comparisons:

1. (3*m-U)/2 > 3900026/10000000.
2. 3900026/10000000 - delta - 1/2^25 - 1/2^980 - 1/r > 39/100.
3. (K*31072+2^90)^25 < 20000^25 * 2^3203.
4. K<r and K+1<2^129.
5. 4*r*3^256*2^25 < 8^256.

The first and second imply success strictly greater than the declared0.39
under the sole disclosed SHA random-function heuristic. The third implies
the full cold computation bound2^(3203/25), with C fixed at20000. The
fourth supports coefficient nonsingularity and the inherited compression-loss
bound. The fifth is the inherited explicit bad-table tail check. All losses
and failed outcomes are charged. Rational score3203/25 equals128.12
exactly; it is an upper bound on log2(total work), not a rounded measured score.

For orientation only, the complete-work logarithm is approximately 128.1191610577 and the success lower bound after all losses is approximately 0.39000259039. The exact inequalities above, rather than these decimals, establish the claim. The exact delta is below 2^-24. No production run is claimed.
