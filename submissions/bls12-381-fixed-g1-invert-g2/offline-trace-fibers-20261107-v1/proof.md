# Offline trace-fiber collision search

Proposed reduced time_log2=31867/250=127.468. Complete full-size fixed-G1 pairing inversion has success>0.39024467 under the disclosed standard-BLAKE3 bridge H and explicit finite banked-storage model. This analytic candidate is not an official score,qualification,promotion or production execution. Full protected p/r,supplied nonidentity P,six tracks,one computation metric,threshold0.39,C20000,zero advice and no supplied auxiliary powers are preserved. All tables are generated and paid per target. This version sacrifices some speed to remove ordinary field arithmetic,online dictionary/head initialization,and every walk/mixing/interval lemma.

## Concrete algorithm

Validate P and GT z. For z=1 return canonical G2 infinity and verify. Otherwise h=pairing(Q0,P) generates GT for EVERY allowed P. Write z=h^a only in the proof,a!=0;the algorithm never receives a. A verified a returns Q=[a]Q0 for the original P,pinned API pairing(Q,P).

Let W=2^92,R=min{integer R:2R^2>=r},L=ceil1031R/1000. Sample2W independent full Fr values v_j,w_k through the actual standard SCAL protocol in hash.md,at most1024 attempts per scalar. Construct T0[j]=h*z^v_j,T1[k]=z^w_k,storing complete640-byte initialized records:rawGT576,canonical scalar32@608,unused zero padding. Binary exponentiation costs<=512 dense references per entry,including h where needed. Every scalar/hash/rejection/selector/record/control cost is in coldD below;these powers are not auxiliary inputs or advice.

For t=0..L-1,full standard PAIR hash selects disjoint independent92-bit indices j,k. Copy T0's raw/scalar data to shared scratch BEFORE switching banks. Multiply by raw T1[k] with ONE dense reference and add beta=v_j+w_k modulo full r before leaving that bank. Thus X=h*z^beta and A=1 for EVERY sample. Samples are not asserted independent conditional on tables;the finite rank proof below handles reused indices.

Extract key S=Tr_(Fp12/Fp4)(X)/3 by copying native a0 and b1,192 bytes total. This identifies exactly the THREE-element orbit{X,X^(p^4),X^(p^8)},as trace.md proves. It does not restrict GT to a smaller subgroup:all exponents and points still have the full protected order r. Abort on identity. Compute the NH160 fingerprint U20 from hash.md. Write an initialized256-byte record(S192,beta32,U20,zero12). Store all L records.

Use sort.md's ten stable16-bit radix passes to sort records by all160 fingerprint bits. This is deterministic counting/scattering on native128-bit counters,not an ideal sorting oracle. Scan equal-U groups and enumerate their unordered record pairs. Confirm all24 exactS words;unequal S is just a fingerprint collision and is ignored. For equal S,compare all four beta words. Equal beta is exactly the same raw point and is ignored cheaply. Different beta gives a useful trace-fiber collision. Explicitly cap pair confirmations at Q=floor(L/16);cap exhaustion fails. All generation,sorting,scan,ignored pairs and failed paths are priced.

For a useful pair,the actual relation is Xt=Xs^(p^ell) for ell4 or8;ell0 is excluded by beta inequality. Try BOTH public choices. For each,compute D=beta_t-xi^ell*beta_s,xi=p_field modr. Skip D=0. Otherwise set a_candidate=(xi^ell-1)/D and verify h^a_candidate=z. The true ell has D!=0:actual equality gives h^(1-xi^ell)*z^D=1,and D=0 would contradict prime-order h and xi^ell!=1. At least one candidate therefore verifies. Return canonical Q=[a]Q0 only after curve/subgroup and exact pairing(Q,P)=z checks. No raw reconstruction,orientation minimum,unchecked scalar-only output or random-denominator premise is used. Both fixed attempts and verification are cold paid.

## Finite useful-collision probability

For fixed input,z=h^a with a!=0. Define A_j=1+a*v_j,B_k=a*w_k;these are independent uniform Fr variables in two disjoint labeled families. A sample raw exponent is x_t=A_J+B_K. Distinct parameter types(J,K) give jointly independent uniform raw pairs,even if one family index is shared. Let C count time pairs with distinct types,nonzero x_s,and x_t=alpha*x_s for alpha in{xi^4,xi^8}. These are two distinct NONIDENTITY phases. Hence
m=E C=binom(L,2)*(1-1/W^2)*2*(r-1)/r^2.

For two useful-collision constraints,expand their two phase choices. A row x_v-alpha*x_u is a nonzero homogeneous linear form in A/B variables. Every variable in its labeled support has nonzero coefficient:1,-alpha or1-alpha. Independent rows have probability1/r^2;dependent rows have probability1/r. Dependence requires identical A-support AND B-support sets. For a three-time shape sharing one endpoint,the third parameter type must equal the other old type,probability1/W^2. For four times,the other ordered type pair has at most four possibilities matching both supports,probability<=4/W^4. These are necessary conditions only,so overcounting is safe. There are exactly three two-edge trees on three endpoints and three matchings on four. Thus
U2=3*2^2*[binom(L,3)*(1/r^2+1/(rW^2))
             +binom(L,4)*(1/r^2+4/(rW^4))]
>=E binom(C,2).
No generic-position,independent-collision,Poisson or walk hypothesis is used. Pointwise1[C>0]>=C-binom(C,2),so useful collision probability>=m-U2. Exact rationals give m=.5314805,U2=.141235760940125. If C>0,a corresponding equal-S/different-beta pair occurs in an equal-U group and is enumerated unless a separately bounded failure intervenes.

## Exact confirmation-cap and failure union

Let D count equal-S time pairs in the full extended sample sequence. Equal parameter types have probability1/W^2;distinct types give independent uniform raw values and trace-key equality probability<3/r. Thus E D<=binom(L,2)*(1/W^2+3/r). Let F count unequal-S pairs with equal U. The NH160 collision bound delta=27/2^160 holds for every fixed unequal key pair,and its seed keys are independent of the complete table/PAIR path. Therefore E F<=delta*binom(L,2),without adaptive-query assumptions. Every confirmed group pair is in D or F. Exact Markov arithmetic gives
Pr[D+F>Q]<=E(D+F)/(Q+1)<2^-24.
Fingerprint collisions are CONFIRMED,not mistaken for key equality or assumed absent. Their work is in the explicit cap. Identity abort is<=L/r<2^-100. Cold scalar-cap failure<2^-970. The union,without failure independence,gives
m-U2-2^-24-2^-100-2^-970>0.39024467945523>0.39.
The identity target is solved directly. Under H,this bound applies for every allowed P and nonidentity target,and hence uniform full-production inputs. H for the actual full standard hashes at the full query scale and finite storage admissibility remain material uncertainties.

## Complete computation

Generation native bounds per sample:key copy/test128,NH2048,PAIR2048,full Fr addition256,all main/two-table/record/control640,total5120. Sorting uses ten passes,each<=64L counting+160L stable scattering=2240L. Group scan<=128L. Each full pair confirmation<=512;Q<=L/16 gives<=32L plus one fixed failing-check allowance. Total ordinary scalar<=7520L+fixed,retained8192L with672 spare. All memory/address/control/serialization/full limbs and every ignored/failing path are included;details in sort.md/hash.md. Exactly one dense reference per generated sample is charged separately. No ordinary field multiply or Frobenius computation is hidden in the trace projection.

ColdD=2^120 instruction-equivalents pays both tables/all rejections/code/constants/seed/buffer/controller setup,radix count-array zeros and prefix sums,every fixed remainder and both verified terminal attempts. Two2^127-record256-byte sort buffers plus tables/code/scratch give peak<2^137 bytes. No online head array or free bulk-memory primitive exists. There is no cross-target amortization.
W_total<=L+(8192L+2^120)/20000=(28192L+2^120)/20000.
Exact(28192L+2^120)^250<20000^250*2^31867 proves proposed127.468,above diagnostic127.46787527606543. It is weaker than unqualified faster proposals but below the published AI-valid unpromoted127.756 comparator,with a much simpler native/key implementation. No previous candidate qualifies this one.

## Material sources and execution status

Direct local material predecessor1d975e82-0f5f-4d10-828d-8a313c9f14f6 supplies the two-random-table index-support collision argument and constant-A/nontrivial-twist extraction. The present variant uses only the raw three-fiber trace,offline radix processing,and a work-cap argument that confirms short fingerprints exactly. a3e0fdf5-5da1-41bf-890f-37795da2a627 supplies the earlier independent/table sampling route;4e784401-a15a-415f-b5ad-f8079fb64135 and013773fa-ec39-496b-b07b-1e0a04de70f0 compact accounting;11a0ed66-fac6-43f9-a7d5-573e7c234e1a exact native trace-fiber mechanism;bf5be6c0-7dd5-4acb-9f20-d2597efc47a6 NH/compression method. The public note credits the full earlier material chain and deferred unpublished links. Published23d/13f supply supplied-P/cold-accounting and full BLAKE3/rejection/Fr-add foundations;no old walk,full-twelve minimum,optimized native body count or qualification is inherited. External XTR/NH and pinned primary definitions are credited. Every candidate file is static untrusted evidence,never host-executed or run with credentials;shared llm-verifier owns scientific execution/review/scoring.
