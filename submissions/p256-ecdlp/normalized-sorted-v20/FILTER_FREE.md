# Filter-free labels and complete load domination

Write R0=2^256, delta=R0-p and k for the canonical Montgomery x-key of an ideal nonidentity leaf. Set U(O)=empty. The hot arithmetic and its denominator/identity exception charge are unchanged.

## 1. Exact output set

The raw Montgomery square satisfies 0<=T<R0+p. The common-sum representation c lies in[0,R0), hence z=T-c lies in(-R0,R0+p), and h=floor(z/R0) is -1,0 or1. The assembly stores u=(z-h*p) modulo R0. Its canonical field answer is k=z modulo p.

- For h=-1 and z>=-p, u=z+p=k. If z<-p, then u=z+p+R0=k+delta, with k=z+2p>p-delta.
- For h=0, u=z. If z<p this is k; otherwise u=k+p and k<delta.
- For h=1, u=z-p lies in[delta,R0). It is k if below p; otherwise it is k+p with k<delta.

Thus every output belongs to U(k)={k}, with k+p added when k<delta and k+delta added when k>p-delta. Since delta<p/2, the exceptional bands are disjoint and |U(k)|<=2. Every extra label is in[p,R0), not arbitrary machine garbage. This covers the noncongruent negative-wrap case as well as redundant congruent residues. The emitter now outputs u without a representation test, subject only to its fully charged slot cap.

## 2. Witness exception loss

Let Bbad=[0,delta) union(p-delta,p). A true negation collision shares its canonical k. This set has at most2delta keys and4delta curvepoints, each with at most two negation-collision partners. For two independent cross-node uniform points, a collision with k in Bbad therefore has probability at most8delta/r^2<2^227/r^2. Summing over Gamma gives the retained Gamma*2^227/r^2 allowance. It is applied unconditionally, before any capacity/survival conditioning. No global S*2^-32 term is omitted: exceptional competitors are handled by the set-valued load bounds below, not excluded globally.

## 3. Single-cell bounds

For a fixed f-bit pattern let h be the number of its fixed positions in[224,256). Because p>R0-2^224, every extra label has all its top32bits set. The containing interval[R0-2^224,R0) has at most2^(224-f+h) integers satisfying the pattern. Each extra raw label has at most two inverse translated canonical keys, and each canonical key corresponds to at most two curvepoints. Extra point preimages are therefore at most4*2^(224-f+h).

For the35-bit slot masks and96-bit home/bucket mask, h=0. FOURIER_CAPACITY.md bounds canonical point counts. In the96-bit case its error is below4721808*2^128, while the new extra fibers add at most4*2^128. Their sum4721812*2^128 is below2^151, the spare1/512 fraction of2^160. For35bits, canonical error plus4*2^189 is also below2^221/512. Consequently both refined coefficients remain513/512.

For home70, the crude canonical count2*2^186 plus4*2^154 gives coefficient2+2^-30. The group144 and matched182 masks include sixteen fixed positions in[224,256), so their coefficient is2+2^-14. Divide these counts by r-1 where a free generator excludes one group point. These are mask-specific bounds, not claims about arbitrary f-bit subsets.

Use the fixed indicator A_C(P)=1[U(k(P)) intersects C]. It is a deterministic function of the ideal point, even when the actual representation depends on a whole batch. Every actual record in C is dominated by this indicator. Define all ideal variables under the original independent generator distribution; do not condition on clean execution or chosen raw representatives.

## 4. Pair loads and dependence

For a fixed opposite ideal point P', its set U(k(P')) has at most two labels. A free point P whose possible labels can match P' on182bits lies in a union of at most two fixed cells. Hence the matched-degree bound is M=2*KB182. Home and relay-cell factors use the single-cell KB bounds. With rho=(r/(r-1))^4, retain the source's conservative shared-generator cases with

h=rho*(S^2*KB70*M+(2S+20)*M),
g=2*rho*n0*S*M,
v=18*rho*((S/SIDE)*S*KB35*M+(2S+20)*M).

The witness is charged separately if bad; otherwise its actual key is its unique canonical label and its home/group are fixed by the witness generators. Thus its intruder rate is rho*(S+12)*KB144, without another factor2. Slot and bucket indicators remain fixed point predicates, so the source covariance argument applies: independent incidences/free Y cases, the one-excluded-point bound for a common X, invertible two-sign forms for a shared incidence, and the conservative60-partner allowance for witness-generator incidences. The phase slot uses qphi, while the persistent bucket uses all q1.

Actual partial/full matches can be false, including between two exceptional labels. They are nevertheless members of the ideal all-unordered-pair envelopes above. Caps and verification account for them. Removing records can create adjacent pairs, so no subset claim about an earlier adjacency list is used. Source labels, signs, generator dependencies and exact final verification are unchanged.

## 5. Work and evidence

Deleting the two CMP/B.HS pairs and one threshold initialization changes the emitter upper bound from33n+12 to29n+11. The arithmetic is unchanged. Its full-record envelope already charges up to two records for every incidence; formerly filtered outputs are not free new work. Exact numerical success and cost rows, including the enlarged h/g/v bounds, are recomputed in success_probe.py. Boundary and differential tests explicitly assert membership in U(k), not equality of every raw output with k. These finite tests supplement this argument; they do not independently prove production success.
