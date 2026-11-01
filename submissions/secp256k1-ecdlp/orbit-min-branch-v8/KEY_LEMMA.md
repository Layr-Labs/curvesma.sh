# Alternative key classifier: replacement lemma

This is the class-key component of the complete argument in proof.md and LEDGER_EVIDENCE.md. Let W=2^64, R=W^4, C=2^32+977, p=R-C, L=W^3. Constants B,D=A-B satisfy B beta=D and B beta^2=-A modulo p, beta of order3. The alternative class key is min(Bx,Dx,-Ax) modulo p.

## Global weak-input approximation

For c=B or D, write c=c0+c1W<W^2 and x=sum(x_i W^i)<R. KEY_APPROX computes F=lo(x3c0)+lo(x2c1)+hi(x2c0)+hi(x1c1) modulo W. The omitted carry into word3 is the floor of

frac(x2c0/W)+frac(x1c1/W)+(x1c0+x0c1)/W^2+x0c0/W^3.

It is nonnegative and bounded above by4-5/W+1/W^3<4, hence is an integer e0 in0..3. Set z=cx,H=floor(z/R),t=(z modR)+HC. Since H<W^2 and W^2 C<W^3, the fold adds at most one word3 carry e1. Thus t/W^3 is congruent to F+e0+e1+theta modulo W,0<=theta<1. Also t<2p. For v=cx modp and V=Wv/p, V is congruent to t/W^3+epsilon modulo W, with epsilon=tC/(pW^3)<2C/W^3<1. Consequently V=F+b modulo W for some0<=b<6. This covers every weak representative x<R, including x>=p and fold-overflow cases.

## Classifier correctness

Let F0,F1 approximate normalized Bx,Dx; sigma=(F0+F1) modW, m=min(F0,F1). Exclude windows F_i>=W-16; sigma>=W-32; (m+sigma) modW>=W-64; and cyclic F0-F1 distance<=32 (the implementation's asymmetric endpoint inclusion is contained in this conservative closed window).

Outside the individual windows, V_i=F_i+b_i as ordinary numbers. Outside the sigma window, adding b0+b1<12 does not cross the relevant modular boundary. The third normalized residue is V2=W-sigma-b0-b1. If m+sigma carries, V2 is below both V0,V1. Otherwise the excluded e-window gives m+sigma<W-64 and V2-m>52; whichever F_i is smaller corresponds to a V_i within6 of m, hence below V2. The difference window prevents the two approximation errors from reversing their order. Therefore the branch classifier gives the exact minimum outside these windows.

## Fixed bad set and endpoints

All windows are contained in fixed circular arcs of canonical field-linear values. In normalized W units, use forms and arc lengths: Bx22, Dx22, Ax44, (2B+D)x82, (B+2D)x82, (B-D)x76. The sum and minimum-plus-sum errors expand by at most12 and18 respectively; the difference expands by less6 at each endpoint. A circular arc of length a contains at most ap/W+2 normalized field residues, charging both endpoints after splitting at zero. Since p/W<L, (a+2)L is conservative. The six contributions are24L,24L,46L,84L,84L,78L, totaling340L<352L. Each coefficient is nonzero, so multiplication permutes the field.

Define Bbad as that fixed union of canonical x-values, independently of weak representatives. Under a uniform nonidentity group point, its probability is at most2*352L/(r-1). Since352L=11*2^197, this is below the retained old allowance2*(13*2^197+2^64 C)/(r-1).

## Meter domination, without conditioning on clean nodes

For an ideal nonidentity point P define J(P)=1 if the exact new minimum uses index2 OR x(P) belongs to Bbad; set J(O)=0. Exact endomorphism orbits allocate index2 to one-third of nonidentity points. Thus E[J]<=1/3+d_wj and Var(J)<=2/9+d_wj. On a computationally clean node every counted index2 path is dominated by J of its ideal leaf; fold drops only reduce its count.

For the unconditional ideal leaf sum in one node: disjoint incidences are independent; a shared Y and independent X are independent; a shared X with distinct nonzero Y gives covariance at most1/(4(r-1)^2); different leaves of one incidence are an invertible linear image of(X,Y), with Y excluding identity, giving covariance at most1/(r-1). Bounding all ordered pairs by the latter yields aggregate less than n0^2/(r-1)<1 at the fixed parameters. Therefore Var(sum J)<=(2/9+d_wj)n0+1. This supports the existing conservative Chebyshev slack2^42-6q2-1-n0*d_wj. Computational defects are union-charged separately by q1*d_batch; do not condition this variance proof on cleanliness or on all four collision-witness generators.

Wrong-index and meter failures are node-local. Conditional on the defective node's coins, the other witness leaf remains uniform, so its class-collision probability is at most6/r. This permits the existing cross-node loss coupling, without a union bound over all leaves. Fold-drop/noncanonical-output exception charges remain necessary and unchanged.


# Transfer of capacity and intruder bounds

This argument supplements the classifier and meter lemmas above. Staying inside an endomorphism class alone does not establish the old covariance assumptions.

Retain the full global arithmetic exceptional-set charge, including LT_BAD=(2^66+2)C0 and KEY_BAD=3C0. Define the bad event by membership in these conservative sets, not merely observed execution failures. A congruent weak leaf >=p requires true x<C0, covered by LT_BAD. Noncanonical key-output exceptions are inverse images of the relevant near-zero/near-p intervals under B,D,-A; their union is bounded by3C0 regardless of the selected branch. The retained global bound is g_def=NODES*q1*d_batch.

Unconditionally define the ideal point-only map F(P): use the unique canonical x(P), run the deterministic new classifier on it, and return the canonical residue c_j*x(P). Identity produces no record. Wrong-index points are included in this definition. Analyze ideal records under the original generator distribution, never the distribution conditioned on cleanliness.

Outside the global bad event every retained actual key equals F(P). Fold, meter and capacity drops only remove records. For each specified key k, F(P)=k implies x(P) lies among k/B,k/D,-k/A: at most six points. Thus a bit cell of at most2^(256-f) integers has at most6*2^(256-f) preimages, exactly the retained KB(f). The factor3 is already included, not an additional factor.

Conditioning on the other ideal point fixes its ideal key even if its classifier index is wrong. The free point must enter one fixed cell, so KB(f) applies once. Existing free-generator home factors and invertible same-incidence linear forms are unchanged. Fixed-point indicator covariance bounds therefore have the same premises as before. Without the global canonical-leaf exceptional-set charge this conclusion would not follow, because a weak representative could depend on the entire batch.

Actual slot/bucket loads are bounded by ideal records and actual candidate/request/relay loads by ideal ALL-UNORDERED-PAIR envelopes. Do not use adjacency-list subset claims: deleting records can create adjacent matches. The inherited h/g/v bounds use the appropriate all-pair envelope. Wrong-index competitors need no extra global charge; the witness's own wrong-index event is still separately charged.

