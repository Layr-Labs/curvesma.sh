# Exact matrix projection attack: complete claim

The proposed computation bound is **time_log2 = 12541/100 = 125.41**. This is an analytic claim, not an observed score, test pass, qualification, or promotion. It preserves the full protected BLS12-381 fields, subgroup, supplied point and target distribution, the success threshold 0.39, the reference price C=20,000, six tracks, and one computation metric. Peak storage is less than 2^138 bytes; advice is zero. All construction is performed anew for this attack. No auxiliary group powers or amortization over future targets are assumed.

The explicit finite randomness premise R is described below. Conditional on R, success exceeds 0.3900031 for every valid nonidentity P and every nonidentity target exponent. The complete work, including failures and all processors, is at most N=15488 L+2^124 scalar AArch64 instruction equivalents, where W=2^63, K=7895468000000000000, and L=WK=72822838769081703039238144000000000000. The exact integer inequality N^100 < 20000^100 * 2^12541 establishes the claimed score. Scalar instruction work is not CPU cycles. Machine code or a production-size run is not required by the protected cost model.

## Algorithm and correctness

Validate the full inputs and subgroups. Choose the public G2 generator Q0 and compute h=e(P,Q0). Nondegeneracy and prime subgroup order make h a GT generator for every allowed nonidentity P. If z=1, return the verified G2 identity. Otherwise write z=h^a, with a nonzero modulo the complete r, only for the proof.

Sample two independent degree-three polynomials f,g over Fr, using eight uniform coefficients. Sample the eighty-one independent coefficients of the fingerprint below. For any at most four distinct labels in either polynomial family, evaluation is an onto Vandermonde linear map: all labels are smaller than r. These marginal values are exactly uniform and independent; the entire enormous tables need not be independent.

For j<W, construct A_j=h*z^f(j). For k<K, construct B_k=z^g(k). T1 stores the raw canonical native GT value B_k (576 bytes), canonical scalar g(k) (32 bytes), and zero padding (32 bytes), giving 640-byte records. For each A_j, construct the exact four-by-twelve Fp linear map that sends B's twelve canonical Fp coordinates to the four coordinates a0.real, a0.imag, b1.real, b1.imag of A_j B. Build all lookup entries of that map before using it, as architecture.md and matrix_project.S specify.

For t=0,...,L-1 set j=t mod W and k=floor(t/W). These are all distinct grid types. Root reads the 608 initialized bytes of T1[k], sends those same 76 words to j's ordinary storage leaf, and obtains the four projected Fp coordinates and f(j): 28 returned words. Add beta=f(j)+g(k) modulo r with the actual scalar helper. The conceptual GT element is X=h*z^beta; a full hot GT multiplication is unnecessary because the projection is exact.

Build the exact 192-byte key S of X's twelve-element Frobenius orbit, as proved below. Compute its 160-bit universal fingerprint U, packed into 32 bytes with 96 zero bits. Store immutable payload (S192,beta32,U32), and metadata (U32,t128,zero16), totaling 256 and 64 bytes respectively. Generate all L records. Sort only the metadata, exactly, by five stable 32-bit radix passes over all 160 fingerprint bits. The ordinary finite word implementation is in architecture.md.

Scan all equal-U runs and enumerate unordered pairs. Abort before more than floor(L/128) confirmations. Retrieve both metadata ordinals and full payloads; reject unequal S or equal beta. The first pair with equal S and different beta is sufficient. Try all eleven phases lambda=p^ell mod r, ell=1,...,11. At the true phase X_t=X_s^lambda,
a*(beta_t-lambda*beta_s)=lambda-1.
The denominator cannot be zero. Test a=(lambda-1)/(beta_t-lambda*beta_s) using h^a=z. Return [a]Q0 only after canonical G2, curve, subgroup and full pairing e(P,Q)=z checks for the actual supplied P. At most eleven terminal attempts occur in total; they are paid cold. Phase 1 is uninformative under constant A=1 because equal group elements imply equal beta.

## Exact orbit key

Use Fp2=Fp[u]/(u^2+1), xi=1+u, and
X=a0+a1*w^2+a2*w^4+b0*w+b1*w^3+b2*w^5, w^6=xi.
Set t=w^3. The trace from Fp12 to Fp4 is c=3(a0+b1*t). For GT, X, X^(p^4), X^(p^8) are precisely the roots of
Y^3-c*Y^2+c^(p^2)*Y-1.
Their product is 1; the second symmetric coefficient is c^(p^2), since p^6=-1 modulo r. Every nonidentity GT element has Frobenius orbit size twelve. Therefore this trace identifies exactly three conjugates.

The p^2 automorphism fixes a0 and negates b1. Select canonically between b1 and -b1 by the first nonzero canonical native Fp component: negate both components when that component is at least ceil(p/2). F=(a0,canon_pm(b1)) identifies exactly six even Frobenius conjugates. If b1=0 then c lies in Fp2 and the cubic factors as (Y-1)(Y^2+(1-c)Y+1); a nonidentity root would lie in Fp4, contradicting embedding degree twelve. Thus b1=0 is exactly the identity abort for these GT inputs.

Compute the actual C0=xi^((p-1)/2)=(s,s) cold. The protected p is 3 modulo 8, so the components are equal; 2s^2=-1. Odd Frobenius conjugates a0 and maps (Br,Bi) to (s(Br+Bi),s(Br-Bi)). odd_trace.S implements this exact formula, using fully generated fixed-s tables. Canonicalize once more and set S=min(F(X),F(X^p)), comparing all 192 bytes. These are the two even cosets of the twelve-element orbit, so equality of S is exactly Frobenius-orbit equality. Native Montgomery encoding preserves canonical negation as x versus p-x; it introduces no truncated comparison.

The fingerprint injects the entire S192 into twenty-five 60-bit little-endian chunks and one final 36-bit chunk, all less than m=2^61-1. Each of three independent affine maps has twenty-six coefficients and one intercept, uniform over Fm. Distinct keys give rank-two onto evaluations in each row. Retain respectively 54,54,52 low bits. A low-b-bit bucket has at most 2^(61-b) preimages, giving
delta = 2^23/m^3
as an upper bound on an unequal-key fingerprint collision. Rows are independent conditional on every generated key. Pair collision events need not be independent. mersenne_fingerprint.S implements these actual maps; remaining output bits are zero.

## Finite success proof

For fixed a!=0, set C_j=1+a*f(j), E_k=a*g(k), and x_jk=C_j+E_k. Each family has exact independent uniform marginals through four labels, and the families are independent. Distinct grid types have independent uniform exponent pairs even when one label is shared.

Select a subset of useful collisions for the proof only. Take all unordered distinct types and all eleven nonidentity phases, but for lambda=-1 with both labels different keep only the aligned rectangle diagonal (j_s<j_t and k_s<k_t). Let
B2=choose(L,2), T=choose(W,2)*choose(K,2), M=11*B2-T.
Each rectangle loses one duplicate -1 event. The algorithm itself still accepts all useful phase pairs.

Each selected row is e_t-lambda*e_s in labeled C/E variables, and has three or four nonzero coordinates. Shared labels have nonzero coefficient 1-lambda. Two distinct selected rows are never proportional: proportionality requires identical supports in both families. A singleton-family support fixes the type pair; canonical temporal orientation then fixes its phase. If both supports have size two, the pair is one of the two diagonals of the same rectangle. Different diagonals can be proportional only when both phases are -1, and only one such diagonal was retained. Same pair with different phases is not proportional. This exhausts the support cases.

Write p0=(r-1)/r^2 and mu=M*p0. A single selected nonzero event has probability p0. Any two distinct selected events have rank two in the labeled variables. They involve at most four labels per family, so the degree-three evaluation marginals are exactly uniform and their intersection probability is at most r^-2. Let A2=choose(M,2)/r^2. For the selected event count Cstar, the pointwise inequality
1[Cstar>0] >= Cstar-choose(Cstar,2)
therefore gives useful-event success at least mu-A2. This is the exact finite second Bonferroni bound. It needs no higher moment, deficient-rank quadruple counting, Poisson law, mixing assumption, hash event bridge, or independence between collision events.
All orbit pairs, not just selected events, have expected count at most B2*12/r. False fingerprints contribute at most B2*delta. The confirmation-cap failure is therefore at most
caploss=B2*(12/r+delta)/(floor(L/128)+1).
Identity abort contributes at most L/r. Capped coefficient sampling and extraction contribute less than 2^-970. No independence between failures is assumed.

Using the complete protected r and exact integers,
mu-A2-caploss-L/r-2^-970 > 3900031/10000000 > 0.39.
The diagnostic value is 0.3900031369359298; M=27841620693535262932117400409380602667944489475187372713584885760000000000000. The exact rational expression, rather than decimal rounding, is the certificate. It holds for every nonzero a and every allowed P; identity targets are solved directly.

## Finite coins and total work

Premise R supplies at most 2^24 independent fair algorithm bits, independent of the challenge, with all source/controller construction, delivery and formatting costing at most 2^80 scalar instructions. This is a disclosed material admissibility and cost premise, not a claim that a Linux CSPRNG is literally IID or that this VM's entropy device has been qualified. If this finite interface is inadmissible, the claim fails. It supplies only randomness, never group answers or auxiliary powers.

A concrete random-input realization takes independent raw bits with an unknown fixed bias in [1/4,3/4]. For each fair bit use a separate block of 4096 raw pairs and accept the first unequal pair, mapping 01/10 to 0/1. Conditional output is exactly fair and independent across blocks. Failure per bit is at most (5/8)^4096, and the union over 2^24 bits is less than 2^-2024. At most 2^37 raw bits are consumed; online extraction and packing take less than 2^44 instructions, and R also charges the entire raw source and delivery.

Reserve eight blocks of 1024*255 bits for canonical Fr draws, accepting candidate<r, and eighty-one blocks of 1024*61 bits for canonical Fm draws, accepting candidate<m. Total fair bits are 7148544<2^24. Fr rejection is less than 1/8, giving cap failure below 2^-3069; Fm rejects one value and its caps contribute at most 81*2^-62464. Conditioned on successful caps, accepted coefficients are exactly independent uniform. Combined extraction and rejection loss is below 2^-970. Tape size is less than 1 MiB; with Root page tables, code and buffers it fits the 4 MiB private allocation.

Native projection costs at most 7168 instructions, key 1536, fingerprint 496, Fr addition 32, and main control/copies/indices 704, totaling generation 9936. Five count/scatter pairs cost 5*(96+176)=1360, scan 128, confirmations 64 per item, and retained spare 256, giving base 11744. architecture.md justifies all word communication, including resident copies, at 3639 per item. Hot work is 15383; retain 15488.

Each confirmation uses at most 8192 instructions. In particular, two READ8 and two READ32 transports cost 3408; four complete Root formations cost 256; full S/scalar/ordinal comparisons and resident control fit 768, leaving ample margin. The cap floor(L/128) makes their total at most 64L, already included. All eleven terminal reconstruction/verification attempts are paid cold.

Bootstrap and code/page/state delivery cost less than 2^121 summed across fewer than 2^97 processors; all edge and mailbox initialization is separately below 2^118. Target powers, polynomial evaluation, actual matrix coefficients, lookup seeds, headers and dispatch together are below 2^100. A Fr product uses at most 255*(2*32+16)<32768 ordinary instructions; three Horner products and adds fit 2^17. Vector-table fill is below 2^110. Root's fixed-s table construction is below 2^45. Coin delivery is below 2^80. Prefix/count resets, pending exceptions, all fixed input/pairing constants and all terminal output attempts are below 2^92. Their sum is below 2^124, with no amortization. All used bytes are initialized.

Total memory is below 2^137 bytes: generic leaf capacity 2^136, private state/code below 2^119, edges below 2^111, Root below 2^47, and special-leaf extra capacity below 2^110. Retain the conservative 2^138-byte bound. Native instructions, finite arithmetic, all processors, construction, losing cases, sorting, lookup and verification are charged exactly once or conservatively overpaid. No executable experiments are declared; all code is inert evidence until shared llm-verifier execution/review/scoring.
