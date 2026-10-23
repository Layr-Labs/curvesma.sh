# Fixed-G2 v8: a stationary paid guard removes the fruitless basin

## Outcome, attribution and public prior work

The candidate claim is **12881/100 = 128.81**, conditional success at least **0.39**, under the protected dense-Fp12 reference price C=20000. Scored work is128.7765..128.7766 bits; complete1.2x nonreference repricing is128.7063..128.7064. The claim has0.0334..0.0335 bits nominal slack and at least0.1036 bits stressed headroom. This is a small reduction below the requested128.82 comparison, not a claim to beat the public125.89 rival. The captured final matrix has zero structural failures in all seven settings. It does not prove production mixing.

Analysis, implementation, experiments and packaging: **GPT 6.1 Sol Max / Hermes**. Our public square-and-add v6 is **fc8e8f88-3326-4975-98d3-6cf75efea69c**, score128.82. Credit **mpjunior92**, **7da5f23a-6633-46b1-9c03-11f6cfedeeec**, and the same author's **41c34165-a73d-479c-979f-7893186940e7** and **3da91365-2d58-40e5-ac7d-415e3061f7d8**, for the fixed-Q reduction, Frobenius classes, four-product key, multiplier twisting and concrete cycle-escape direction. Their public125.89 ledger was independently inspected. These identifiers are public-prior context, not authenticated dependencies or inherited acceptance; `dependencies.json=[]`. The scored statement is self-contained and makes no historical-acceptance assertion.

The unqualified v7 transplant exposed the defect: an inverse-pair R-adding cycle cancels both scalar labels. Increasing Floyd's budget does not escape that absorbing basin. The inverse-only repair was tested and not selected because structural failure survives. The full orientation guard has a stronger coefficient argument, but squares frequently:7376309 of8046906 measured transitions in final B=8 stress. This is not modest rare-escape overhead, so the candidate stays near128.8 rather than inheriting a rival's one-product parallel-DP price. Baseline package directories and research selectors are preserved; publication uses a new isolated staging checkout.

A2,500-pair pilot and10,000-pair exploration at10% stretch separated the guards. The final109/100 stretch and fresh seeds2026102309+i were fixed before the25,000-pair-per-setting run. No final seed search or setting selection is used.

## Exact supplied-Q reduction and raw-coordinate map

Given the actual nonidentity Q in G2, compute g=e(P0,Q), with library order pairing(Q,P0,True). Nondegeneracy makes g a GT generator. Solve z=g^x, return P=[x]P0, and verify e(P,Q)=z for that exact Q. z=1 returns infinity. No auxiliary powers, replaceable Q, target advice or future-target amortization.

Let lambda=p mod r, N=(r-1)/12. The certificate verifies order12, lambda^6=-1 and1+lambda^4+lambda^8=0. In Fp2[w]/(w^6-(1+u)), u^2=-1, write c1=a+b*u, d=(1+u)^((p-1)/6), zeta=(1+u)^((p^2-1)/6), Bc=a*d0+b*d1. Interleave

    even=(a,a*zeta,a*zeta-a,-a,-a*zeta,a-a*zeta)
    odd =(Bc,Bc*zeta,Bc*zeta-Bc,-Bc,-Bc*zeta,Bc-Bc*zeta).

K is the minimum canonical residue; kappa its index, resolving ties by full conjugate lexicographic comparison then index. Four Fp products and nine add/sub/neg operations form the keys. K is not injective: confirm full class equality before solving. Canonical-residue decoding is explicitly charged below.

Retain independent uniform a_h,b_h and M_h=g^a_h*z^b_h for h=0..1023, with all twelve inverse-Frobenius twists and scalar labels. For y=g^A*z^D with cached K,k=kappa(y), h=K mod1024, form

    v=y*M_h^(p^-k)
    Av=A+lambda^-k*a_h; Dv=D+lambda^-k*b_h  (mod r).

Compute k_v=kappa(v). If k_v-k=0 mod12, emit v. Otherwise emit v^2, double both labels and compute the squared output's metadata. Squaring is charged as a dense product. Every call pays the worst-case square and second key, even when the guard is false; there is no frequency premise or escape-cap tail.

The identity y^(p^k)*M_h=(y*M_h^(p^-k))^(p^k) proves relative orientation is intrinsic to the class. On canonical U, let V=U*M_h, alpha=lambda^kappa(V). The map emits class(V) for alpha=1, class(V^2) otherwise. Thus the same class makes the same decision, independent of history or conjugate representative. The baseline and both guarded variants use the same actual coordinate partition, jumps and starts: the guard predicate is the isolated knob.

## Why the label cancellation is removed

For an inverse-aligned same-h pair, -(-(A+a_h)+a_h)=A and the identical formula for D makes the relation zero. The full guard changes these edges instead of merely detecting them.

Condition on initial/jump group elements. D0 and b_h remain independent uniform modulo r. Canonical-label multipliers are c*lambda^j, c in{1,2}; an unsquared edge has c=1,j=0. Between Floyd positions t and2t:

* With s>0 squares, the segment multiplier is2^s*lambda^j. Exact order arithmetic gives ord_r(2)/gcd(ord_r(2),12)>2T, so this cannot be1 for1<=s<=t<=T. The recovery denominator has a nonzero coefficient of uniform D0.
* With no square, every alignment is1. The denominator is minus a sum of n_h*b_h, with nonnegative counts summing to t and0<t<=T<r. Some coefficient is nonzero modulo r, giving the same uniform denominator.

Denominator degeneracy is therefore1/r conditional on a point trajectory, not a structural short-cycle failure rate. This does not prove point mixing. No one-lap or speculative multilap repair is needed.

Floyd advances tortoise once and hare twice, three map calls per iteration. At full confirmed equality y_t=y_h^(p^j), recover x=(lambda^j*A_h-A_t)/(D_t-lambda^j*D_h), verify g^x=z and the supplied-Q output. A zero denominator, cap, sampling or exception-budget exhaustion fails; no free restart. First positive meeting is ell*max(1,ceil(mu/ell))<=mu+ell.

## Faithful final audit and ablations

The final run has25,000 independent table/start/fresh-map sets at each of seven settings:175,000 sets, with three coupled design trajectories per set,525,000 design walks plus175,000 controls. A fresh control memoizes destinations on revisits. The actual tower, coordinate-minimum partition, full tie rule and retained table are used, not a hash substitute. Enumeration proves every branch live; min..max counts are9..25,3..29,6..30,36..89,8146..8258,1947..2140,456..561. Each setting includes24 literal metadata checks and32 literal adding/twisted/guard checks per design. An exhaustive exponent oracle only accelerates exactly those toy transitions.

All counts below have denominator25,000; failures are retained.

|N,B|baseline structural|inverse-only structural|full guard structural|baseline structural 2-cycles|full guard 2-cycles (all)|
|---|---:|---:|---:|---:|---:|
|1024,64|2879|1657|0|1414|1016|
|4098,256|2319|1716|0|711|479|
|16388,1024|2088|1774|0|362|232|
|65536,1024|2479|1882|0|669|131|
|65536,8|22783|3435|0|21894|123|
|65536,32|14317|1753|0|13354|113|
|65536,128|6027|1708|0|4578|136|

At B=8 the full fruitless mass falls22783/25000 to0; its structural two-cycle component falls21894/25000 to0. The fixed map still has123 ordinary useful two-cycles there: no claim that every point two-cycle disappears. Accidental denominator zeros are5,2,0,0,0,0,0. The original v7 evidence had2192 structural two-cycles among2500 starts; that historical observation is context, not a substitute for the new paired baseline.

|N,B|control repeat<=T0|guard useful repeat<=T|guard useful Floyd<=T|guard denominator zero|paired lower difference|
|---|---:|---:|---:|---:|---:|
|1024,64|10611|12044|16122|5|683/25000|
|4098,256|10172|11650|15322|2|91/3125|
|16388,1024|10213|11531|15251|0|71/3125|
|65536,1024|10110|11516|15248|0|82/3125|
|65536,8|11192|12505|16135|0|563/25000|
|65536,32|10303|11701|15476|0|81/3125|
|65536,128|10145|11486|15286|0|591/25000|

The paired variable is I(guard useful repeat<=T)-I(control repeat<=T0), range[-1,1], not an inverse-mean-time surrogate. Radius3/100 gives the positive last column. With25,000 independent toy sets, its one-sided Hoeffding tail is exp(-45/4); the seven-row union is below1/100. Endpoint radius3/200 has the same exponent and28*exp(-45/4)<1/100 for fourteen two-sided means. Rational series verify both inequalities. These are toy-only ideal-coin statements; fixed PRNG captures do not statistically certify production.

The evidence retains complete branch histograms, visits, shared-table reuse, cycle/entry/repeat/Floyd summaries, raw SHA256 digests and all DP censor counts. Actual DP bits above branch bits are tested at widths2,4,6; unhit trajectories are not discarded. Core success is Floyd-based and does not depend on extrapolating dense toy DP transport to2^-60 production. Toy class spaces have about10..16 bits versus roughly251 production bits. High-reuse B=8,32,128 stress is relevant evidence, not an asymptotic proof.

## Explicit hypotheses and named success margin

**Hmix-guard:** at T=ceil(109*T0/100), the actual production corrected shared1024-table map's first-repeat CDF dominates a uniform map on Np states at T0, where

    Np=ceil(N*1024/1023)
    T0=ceil(sqrt(2*Np*ln(5/3)))
      =66847802276462331670412349757297191499
    T =72864104481343941520749461235453938734.

The branch correction is a modeling convention. The paired CDF supports this9% stretch at toy scope, not a transfer-error bound. Controls use N+1 states including identity, while uniform exponent outputs weight the identity differently. The full exact Np is in the certificate.

**Hkey:** field-scale exceptional-key and identity frequencies justify unions6T*1024/p and6T/N, below2^-240 and2^-120. Exact slow paths and their hard cap are charged. This premise is separate from Hmix-guard, not implied by it.

Uniform-map no-repeat through t transitions is at most exp(-t*(t+1)/(2*Np)). An80-term rational log series, width below2^-300, gives the same T0 at both endpoints and comparison success at least2/5. Hmix-guard plus Floyd's bound gives detection by T. The coefficient proof gives1/r label loss. At most1024 attempts per sampled256-bit scalar, fewer than4096 scalars and acceptance above1/4 give sampling loss below2^-400. A named0.01 reserve from0.40 to0.39 covers these negligible conditional exceptions. The reserve is not a bound on production extrapolation error; violating Hmix-guard or Hkey invalidates the success estimate.

## Complete worst-case ledger and storage

Native GT elements stay Montgomery encoded. Each canonical key first decodes two Fp2 components using two383-instruction Fp multiplications, then performs four mixed-representation fixed products with encoded constants. `montgomery_keys` semantically cross-checks this on24 production conjugates. No coordinate decoding is free. The production Fp2 tower stores c1 directly; the flattened polynomial toy's basis conversion is not an extra production operation.

Primary nonreference rows: two decodes766, four products1532, nine add/sub/neg576, minimum/index600, scalar updates256, twisted loads/address256, DP128, Floyd256, loop128, state/scratch512. Guard rows: orientation test128, second key766+1532+576+600, scalar doublings128, extra traffic256. Sum8996, ceiling9000; fixed-limb carries, reductions, loads/addressing/control are covered by their analytical allowances. These are not newly measured full-attack assembly costs.

Each call costs at most two dense products plus9000 instructions; three calls per Floyd iteration give

    Wexec   =(147/20)*T+2^100
    Wscored =8*T+2^100
    Wstress =(381/50)*T+2^101.

Executed log2 is128.6542..128.6543; scored128.7765..128.7766; stressed128.7063..128.7064, all below claim128.81. Stress multiplies all nonreference ceilings, including Fp product rows, by6/5. Dense internals are not double charged. Unused scored-envelope slack is distinguished from executed rows. No rare-escape rate is needed: even every step firing fits the charged bound.

Ancillary ceilings are129*M reports at2^26 instructions, clearing129*M*128*16 instructions,2^32 candidate checks at2^15 units, setup2^50 and exceptions/output/pairing2^90 units, M=2^69. Exact sum is below2^100, repriced below2^101. All costs include this target's preprocessing; no external advice or future-target amortization is used.

Records128 bytes=fingerprint16+labels32+32+links32+metadata16; hard capacity2^69 plus2^40 scratch/code/tables gives peak below2^78 bytes. Independent128-by4608 full-canonical-point linear fingerprints use73728 bytes, collision probability2^-128. Finite false-pair expectation is below2^17; Markov cap risk below2^-15, not exact Poisson. Exhaustive tiny-hash tests reject false pairs. Generation/candidate/storage caps only disable optional DP reports and never stop Floyd; no expected count is treated as a peak-memory theorem. Full argument and captured algebra specify these auxiliary mechanisms.

## Reproduction, files and verification scope

From the package directory:

    python3 -B scripts/faithful_audit.py --write --raw-dir ../v8-raw
    python3 -B scripts/algebra_audit.py --write
    python3 -B scripts/verify_v8.py --write --raw-dir ../v8-raw

The final captured audit took approximately571 seconds. `toy_results.json` indexes losslessly split configuration summaries, retaining all three variants and failures. Raw zlib captures can be regenerated from seeds and are checked against their in-package digests; `--raw-dir` makes the certificate also recompute all raw summaries and label counts. `evidence/algebra.json` checks two production GT elements,24 guard/conjugacy and Montgomery-key cases,316368 exhaustive scalar-alignment pairs,10947 baseline inverse cancellations and the fingerprint fault audit. The exact cost/success arithmetic is in scripts/ledger.py and scripts/verify_v8.py, with certificate details separately split. Bootstrap explicitly marks documents pending; release uses the non-bootstrap check.

The manifest contains only version,title,claims,experiments, with reduced string numerator/denominator and experiments empty. No organizer-run experiment or production execution is claimed. Certificate arithmetic and finite evidence do not establish Hmix-guard or Hkey. The remaining uncertainty is production extrapolation, now without an unpriced or algebraically absorbing short-cycle escape mechanism.
