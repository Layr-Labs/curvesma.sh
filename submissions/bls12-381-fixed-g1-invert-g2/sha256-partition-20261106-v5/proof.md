# Finite pairing inversion with a SHA-256 walk partition

This is an analytic classical attack at the protected full BLS12-381 parameters. The sole score claim is 320583/2500 = 128.2332. Total cold work is less than 2^(320583/2500) dense-Fp12 reference equivalents at fixed C=20000. Success is at least 0.39 conditional on the joint standard-SHA random-function idealization below. Peak storage is below 2^50 bytes, nonuniform advice is zero, and all preprocessing is included. No full-size attack is executed or claimed practical.

## Exact algorithm and sole computational idealization

Read p,r,P0,Q0 and pairing semantics from the protected target profile. Validate the exact supplied nonidentity G1 point P and GT input z, including canonical decoding and subgroup checks. Compute h=e(P,Q0) using pairing(Q0,P) in the pinned library. Nondegeneracy and prime r imply that h generates GT for every allowed P. If z=1 return {"Q":null}. Otherwise there is a unique a modulo r with z=h^a. The eventual answer is the actual G2 point [a]Q0 for this supplied P.

Use the 16-byte ASCII tags FAPI-RHO-INIT-v1, FAPI-RHO-SCAL-v1, FAPI-RHO-UKEY-v1 and FAPI-RHO-WALK-v1, abbreviated INIT, SCAL, UKEY, WALK. Canonical P_bytes consists of its two 48-byte coordinates; z_bytes consists of its twelve 48-byte coefficients. Let seed=SHA512(INIT||P_bytes||z_bytes). Initialize one unsigned 64-bit counter to zero. Every sampler attempt hashes its tag||seed||counter_be64 and then increments the counter. There are no external random coins or advice services.

First sample five compression keys by taking the low 61 bits of SHA512(UKEY||seed||counter_be64), rejecting values >=q=2^61-1. Sample all scalar residues by taking the low 255 bits of SHA512(SCAL||seed||counter_be64), rejecting values >=r. Limit each sampler to 1024 attempts. Exhaustion terminates with failure. Use this scalar sampler for independent pairs (u_j,v_j), j=0,...,w-1, with w=2^40, then alpha0. Construct M_j=h^u_j*z^v_j. Each table record stores the full 576-byte canonical Montgomery GT element and both 32-byte residues, 640 bytes total. Initialize (X,A,B)=(h^alpha0,alpha0,0). The invariant is X=h^A*z^B.

For a canonical GT element X, compute the injective 289-byte norm-one encoding E(X) and the keyed 40-byte compressor U(E(X)) specified below. Define j(X) as the low 40 bits of SHA256(WALK||seed||zero48||U(E(X))), where zero48 means exactly 48 zero bytes. Execute one transition by X<-X*M_j, A<-A+u_j mod r, B<-B+v_j mod r. The label depends on X and public instance initialization, never its coefficient history. The walk uses standard SHA-256; seed, keys and scalar samplers use standard SHA-512. All hash computations are actually performed as scalar instructions and paid for. Do not substitute SHA-512 for the walk in this revision.

Let R0=isqrt(r-1)+1, L=ceil(1031*R0/1000), s=ceil(L/62), K=L+s. Save index-zero state. For t=1,...,K, advance once and binary-search a sorted array of pointers to saved complete GT records. Compare keys lexicographically as 72 canonical unsigned 64-bit words in a fixed order. On first equality to an earlier record, set D=B_current-B_saved mod r. If D=0 fail immediately. Otherwise a=(A_saved-A_current)/D mod r. Verify h^a=z, construct Q=[a]Q0, verify its subgroup and e(P,Q)=z, then output {"Q":Q} in the protected encoding. If there is no equality and t is one of s,2s,...,62s, save the full current 640-byte record and insert its pointer at the sorted insertion position. Never search a checkpoint created at the current t. There are at most 63 saved records. At K with no useful relation, fail.

The only computational heuristic H is: for average-case success over the protected uniform independent P,z, these domain-separated SHA-512 initialization computations and SHA-256 walk computations may jointly be analyzed as independent random functions on their distinct domains. This is unproved for the fixed standard functions at this enormous scale. It explicitly includes both algorithms; no unconditional theorem or observed production success is asserted. It is the familiar hash idealization for a fully computed adding-walk partition, not a supplied free random oracle. SHA-256 changes the function family from earlier SHA-512-only candidates, so this joint idealization must be reviewed explicitly. All finite probabilistic bounds below are proved inside that idealization; no independent random-mapping assumption for the adding walk is made.

## Compression, exact sign, and scalar audit

Use Fp2=Fp[u]/(u^2+1), Fp6=Fp2[v]/(v^3-u-1), Fp12=Fp6[w]/(w^2-v). This matches the protected polynomial basis under u=w^6-1 and v=w^2. Represent base coefficients by six canonical Montgomery limbs. Convert protected inputs once in charged setup. Call exactly the pinned reference dense Fp12 multiplication for each walk step. Its output coefficients are canonical. No custom tower-product price is used.

Write X=a+bw with a,b in Fp6. On GT, p^6=-1 modulo r (an exact protected-parameter check), so X^(p^6+1)=1 and a^2-vb^2=1. For fixed b there are at most two choices a,-a. Let the encoding be b's six canonical 48-byte coefficients followed by a sign byte: one iff the coefficient tuple of a is lexicographically greater than that of -a. If a=0 the sign is zero. This distinguishes all elements injectively, including identity; full 576-byte elements remain stored and compared.

Compute this same sign without any field negations. Inspect coefficients in the prescribed tuple order. Zero coefficients equal their own negations. For the first nonzero canonical Montgomery integer x, x>p-x iff x>floor(p/2), equivalently x>=H0=(p+1)/2. At each coefficient load six limbs (six LDR), OR them (five ORR), and CBNZ to its comparison label: at most 12*6=72 instructions. At the selected label load H0's six limbs, subtract low-to-high with SUBS and five SBCS, then CSET CS: 13 instructions. This comparison uses the final no-borrow flag, not a single-limb zero flag. Fixed pointers, optional caller BL, zero-path handling, sign store and return fit 32 more. Total <=117; charge 128. This leaf needs only caller-saved registers. The constant table is prepared in paid setup. Montgomery representation causes no exception: canonical modular negation is still p-x for nonzero x.

The compressor modulus q=2305843009213693951 is prime. A complete Lucas certificate has base 37 and q-1=2*3^2*5^2*7*11*13*31*41*61*151*331*1321. The listed factors are prime by trial division through 37; modular exponentiation checks 37^(q-1)=1 mod q and gcd(37^((q-1)/f)-1,q)=1 for every distinct factor f. Charge certificate construction/checking in setup. Interpret E as a big-endian integer M<2^2312 and form 39 digits c_l=(M>>(60*l)) & (2^60-1), l=0,...,38. These fixed-length digits determine M. Define f_M(t)=sum c_l*t^l over Fq. U(M) is the five evaluations f_M(k_i), each as eight big-endian bytes. Different M give different polynomials, so Pr[U(M)=U(M')] <=(38/q)^5 for fixed unequal encodings.

Do not union-bound an already key-dependent walk. Use the key-independent independent-index surrogate constructed in the finite collision proof below. Conditional on that entire surrogate, union-bound its at most binom(K+1,2) distinct-state pairs. Since K+1<2^129, q>2^60 and 38<64, the probability of any unequal-state compression collision is below 2^257*(2^-54)^5=2^-13. Sequentially couple the ideal-hash walk to that surrogate until its first group repeat or distinct-state compression collision. Before either event every compressed input is fresh and its 40-bit SHA-256 label is uniform independent. Subtract only this bounded compression event. No additional compression heuristic is introduced.

For each of five Horner evaluations, every product-plus-digit update h<-h*k+c mod q uses this exact scalar sequence, with 64-bit registers:

    MUL lo,h,k
    UMULH hi,h,k
    ADDS lo,lo,c
    ADC hi,hi,xzr
    LSR tmp,lo,#61
    ORR tmp,tmp,hi,LSL #3
    AND lo,lo,q
    ADD lo,lo,tmp
    LSR tmp,lo,#61
    AND lo,lo,q
    ADD lo,lo,tmp
    SUBS tmp,lo,q
    CSEL h,tmp,lo,HS

Here h,k<q and c<2^60. The 128-bit sum is folded at bit 61 using 2^61=1 mod q; first folded sum is <2q and the second <=q, so one conditional subtraction canonicalizes. Both halves of multiplication and all carries are charged. Fully unroll 39 descending digits shared across five accumulators: 195*13=2535 instructions. Aligned digit extraction uses at most two loads, two shifts, an OR and mask; charge eight per digit, 312 total. Five accumulators, five keys, q, digit, lo, hi, tmp, two pointers and extraction temporary fit at most 19 registers. A reserved x18 can use one saved pair; no per-update spills or calls are needed. All initialization, key loads, saved pair, constant creation, output reversals/stores and return fit 128. Thus the complete compressor <=2975; charge 3072.

The aligned digit buffer is prepared in the separate hash preparation allowance. Read b as big-endian words B_0,...,B_35; form d_0=(B_35<<8)|sign, d_j=(B_(35-j)<<8)|(B_(36-j)>>56) for 1<=j<=35, and d_36=B_0>>56. Load/reverse each B once, retain its neighbor, and store all 37 aligned words. This takes <512 instructions including pointers and control. Offsets are fixed below 296 bytes. No unaligned read or uninitialized padding occurs. It reconstructs M exactly.

## SHA-256 implementation and its complete price

Use the exact standard algorithm in [FIPS 180-4](https://csrc.nist.gov/files/pubs/fips/180-4/final/docs/fips180-4.pdf), sections 4.1.2 and 6.2.2, for walk hashing. The initialization remains SHA-512. Standard padding is preserved. WALK||seed||zero48 is 128 bytes, exactly two SHA-256 blocks; compute those two fixed blocks once in paid setup. The remaining 40 payload bytes, 0x80 byte, zero padding and eight-byte original-message length (1344 bits) form one final 64-byte block. This is ordinary prefix precomputation, not a custom compression hash or a modified IV.

Use A64 scalar W-register integer instructions, never SIMD or SHA extensions. [Arm's instruction reference](https://documentation-service.arm.com/static/6245c734b059dc5ff9a8bdab) permits logical shifted operands including ROR. Each big sigma is three instructions: one ROR followed by two EORs with rotated operands. The Sigma0 rotations are 2,13,22; Sigma1 rotations 6,11,25. Each small sigma is likewise three instructions: sigma0 uses ROR 7, EOR ROR 18, EOR LSR 3; sigma1 uses ROR 17, EOR ROR 19, EOR LSR 10. W writes implement arithmetic modulo 2^32 without extra masking.

An expanded schedule word needs four loads, two small sigmas (six instructions), three ADDs and one store: 14 instructions. Fully unroll 48 such words with immediate frame offsets below 256 bytes. One round needs two big sigmas (six), choice (XOR(f,g),AND(e,tmp),XOR(g,tmp), three), majority (ORR(a,b),AND(c,tmp),AND(a,b),ORR of these terms, four), two constant/schedule loads, and seven additions for T1,T2,new A,new E. Total 22. The working-state rotation is a compile-time register permutation: overwrite old H and D with new A and E after their last uses. Charge 24 per round, including two extra instructions of margin. There are 64 unrolled rounds. Eight state registers, five reusable temporaries, two word loads and three base/output pointers fit 18 registers, below x18; expansion reuses registers before state initialization. No hidden stack spills or runtime operation loops occur.

Fixed overhead: 16 LDR/REV/STR triples for initial schedule words cost 48; eight state loads cost 8; eight feed-forward LDR/ADD/STR triples cost 24. Total 80. Stack frame, pointer setup, epilogue and return fit the remaining 48 of a 128 fixed allowance. No nested call occurs. Keep the entire old 2048 preparation allowance: <512 for the 37-word compressor input; <512 for chaining-state copy, padded block, label extraction and fixed transfers; remaining >1024 for caller preparation, pointers, spills and formatting. Low 40 digest bits can be assembled from the final two 32-bit words with a mask and shift/OR, all included. Nothing assumes a free arbitrary-length hash.

Complete per-step SHA cost <=64*24+48*14+128+2048=4384 scalar instructions. INIT, UKEY and SCAL hashes, and both prefix blocks, are separately covered by the fully charged cold setup. Their standard algorithms and domains are not altered.

## Checkpoint detection, relation, and remaining scalar work

For a deterministic map with tail mu and cycle lambda, first repeat is tau=mu+lambda. If tau<=L, the checkpoint a0=ceil(mu/s)*s exists among 0,s,...,62s. It lies in [mu,mu+s-1], and a0+lambda<=tau+s-1<K. Thus the algorithm detects a previously saved equality by K whenever the first repeat was by L. The sorted search finds every equality; saving only after search avoids an index-zero-length relation. No hypothesis about tail/cycle proportions is used.

Conditional on jump group values M_j, every v_j is still independent uniform modulo r: each fixes exactly one u_j with M_j=h^u_j*z^v_j. The trajectory, sorted key order, checkpoint schedule and equality selection inspect only group values. At the selected equality, coefficient difference is sum n_j*v_j, where nonnegative n_j count jumps between indices and sum to d0 with 0<d0<=K<r. At least one coefficient n_j is nonzero mod r, so the denominator is uniform even conditional on the selected trajectory. Zero occurs with probability exactly 1/r. Otherwise the invariant proves the unique recovered a and verified [a]Q0. This introduces no separate useful-relation heuristic.

At most six full-key comparisons suffice among 63 saved keys. Unroll 72 words, four instructions per word (two LDR, CMP, B.ne). Allow 96 more per comparator/search iteration for its call/return, bases, pointer lookup, midpoint/lo/hi and control: 384. Search costs at most 2304 even if 71-word prefixes agree. Records are immutable; sorted pointers are updated only after an unsuccessful equality search. At most 63 insertions shift at most 62 pointers; even 16 scalar instructions per moved pointer is far below the 2^30 fixed insertion allowance, paid in setup/fixed work rather than amortized over other targets.

Each modular scalar addition loads two four-limb operands (four LDP), uses ADDS and three ADCS, loads r (two LDP), subtracts with SUBS/three SBCS, selects four limbs with CSEL, and stores with two STP. This is 20 instructions. Since r<2^255 the sum has no unrepresented fifth carry. Allow 128 pointer/control/spill instructions per addition; both updates <296, charged 512. Reserve also 256 for table index/addressing, 512 for counters/checkpoint control, 512 for an occasional complete 640-byte record copy, 256 for equality/dispatch and 1024 for other spills/control. The full outer sum is 2304+512+256+512+512+256+1024=5376. Multiword step counters have at most 129 bits; all their limb work fits these explicit control allowances. Pointer insertion is already paid separately.

One step uses the exact named reference dense Fp12 multiplication, H_step=1. The protected classical-group-work-v2.json and COST_MODEL.md expressly allow H+R/C, with a reference operation costing one unit. Do not replace C with a candidate timing or charge the multiplication twice. Other work is R_step=3072 compressor+4384 hash+128 sign+5376 outer+384 reference-call interface=13344. The 384 interface allowance covers at most two scalar instructions for each of three pointers, one caller BL, buffer selection/swap and residual interface control. The reference writes straight to the next 576-byte state; no hidden full product copy or custom tower routine is used. All transfers are charged as loads/stores.

## Cold construction, memory, exact finite arithmetic and attribution

All setup, bounded rejected samples, table construction, parsing/subgroup checks, basis conversion, constants and prime certificate, initial pairing/power, sorted insertion, final scalar inverse, exponent verification, G2 output and independent final pairing verification are charged within 2^90 scalar instructions per cold input. At most (2w+6)*1024<2^52 sampler attempts, at <2^16 instructions per scalar SHA-512 hash, cost below 2^68. Each of 2^40 jump entries uses <1024 dense multiplications for two binary exponentiations and a product. Even a complete dense implementation/interface bound 2^20 instructions gives total below 2^70. Such a bound is elementary for fixed six-limb Montgomery products and tower formulas and is much looser than the protected measurement. No large-operand arithmetic is hidden there.

A base inverse by p-2 needs <768 products and Fp12 inverse by p^12-2 needs <9216. The pinned 64-bit Miller loop plus <4608-bit final exponent uses <2^20 fixed-size extension operations. Allowing even 2^58 scalar instructions per such operation including inversion bounds each pairing below 2^78. Two pairings, table work, hashes and all other fixed work fit far below 2^90. Complete point formulas branch for infinity, equal/opposite inputs and zero y before division; Fp2 inversion uses (a-bu)/(a^2+b^2). At most 4096 point steps with even 2^52 instructions each remain negligible in this budget. All failed outcomes are charged. Preprocessing alone is below 2^76 reference equivalents, included in total, with no future-target amortization.

The jump array is exactly 640*2^40=703687441776640 bytes (640 TiB). Code, constants, sorted pointers, at most 64 current/checkpoint records, scratch and retained seed/keys fit in the remaining margin below 2^50 bytes. This is analytic storage, not allocated VM memory. Use an aligned low base; every array address stays below 2^50. Forty-bit j gives j*640=(j<<9)+(j<<7) below 2^50, ordinary scalar addressing. There is no astronomical-width unit-cost pointer or free communication oracle. Advice is zero; public profile constants are excluded as permitted, while construction of algorithm-specific certificate/constants is charged.

Total H<=K and R<=13344*K+2^90 give W<=K+(13344*K+2^90)/20000=(33344*K+2^90)/20000. For the protected r, exact integer/rational checks are:

1. With D=(L-270)*(L-269)/2, eta=1/65536, m=(1-eta)D/r, and U as proved in the finite collision appendix, (3m-U)/2>39015/100000.
2. 39015/100000-1/8192-2^-100-2^-980-1/r>39/100.
3. (33344*K+2^90)^2500 < 20000^2500 * 2^320583.
4. K<r, K+1<2^129, 4*r*3^256*2^100<8^256.
5. w/(4*32768^2)=256, (2w+6)*1024<2^52, and 640*w+2^30<2^50.

These establish the declared rational work bound and conditional 0.39 success, not an official benchmark score. Decimal work logarithm is about 128.233099835; the exact inequality, not this approximation, is used. The full finite collision appendix follows. No executable experiment is declared and no candidate code executes on the host. Shared llm-verifier owns review/scoring.

Material sources: published submission 755cd10d-5cc2-4b71-8542-11447238cf4b supplies the walk, finite correlated interval bounds, compression and sorted checkpoints and is declared. The immediate direct source is our published larger-spectrum submission 30a7c6df-a9c5-47fe-895b-a5bdff9433cb, accepted at 128.292 but awaiting manual promotion, declared as a dependency. Its algorithm, finite proof and scalar audits are reused. Earlier still-pending or infrastructure-failed own work includes: 2274e1f2-cab0-442d-9ace-c733636ddd62 (retry 7dc67838-0672-4518-b7f8-49165ee7ec8a) for reference/sign implementation, 7f59eb7c-8345-4654-b502-d8a42e304654 (retry 36064d96-ffcf-488b-adaa-849d8a9e548c) for Bonferroni and 62 checkpoints, with their features also carried in the published immediate source. Their relationships will be added when publication permits; no independence is claimed. The new change here is the explicitly disclosed SHA-256 walk function and its complete scalar schedule, with SHA-512 initialization retained. This is exploratory, conditional mathematical evidence, not a production execution, unconditional break, novelty certificate, human approval, qualification or promotion claim.

# Finite correlated collision appendix

# Finite collision bound in the ideal-hash model

Let N be the protected prime subgroup order, w=1099511627776, and M_1,...,M_w independent uniform elements of GT. The implementation generates these from uniformly idealized initialization scalars; no test for a good table or discrete-log oracle is performed. This note proves a finite collision bound in that idealization. It does not assert that the joint fixed-SHA functions is a mathematically random function.

## Good-table probability

For every nontrivial character chi of the cyclic group, define lambda_chi=(1/w)*sum_j chi(M_j). Its real and imaginary summands lie in[-1,1],have mean zero,and are independent over j. If |lambda_chi|>1/32768,one component has magnitude>1/(32768sqrt(2)). Hoeffding's inequality on each component and a union bound give probability<=4exp(-w/(4*32768^2)). Union over N-1 characters gives bad-table probability<=4Nexp(-256)<2^-100 for the protected N<2^256. For an exact check, e>8/3 implies this is below4N*(3/8)^256; the claimed dyadic bound is the integer inequality4N*3^256*2^100<8^256. This is an analytic public-integer inequality; no experimental result is claimed.

Condition henceforth on a good table, so every nontrivial eigenvalue of its uniform-addition Cayley transition matrix has modulus<=rho=1/32768. Duplication among the sampled M_j is allowed and does not invalidate Hoeffding or the transition weights. The attack pays for all samples and accepts bad tables too, counting their possible failures.

## Markov surrogate and coupling

Consider independent uniform indices J_t in[1,w] and the Markov walk Y_(t+1)=Y_t*M_(J_t). A random-function partition assigns independent uniform indices to previously unvisited elements. Therefore its deterministic iterated walk and this independent-index walk can be coupled exactly until the first repeated group element. They need not agree after that repeat. Only the distribution of the first repeat is transferred below. This revision first compresses encoded states with an independently keyed polynomial family. Coupling also stops on a collision of distinct states under that compression. universal-proof.md bounds the probability of that extra stopping event by less than2^-13 on the independent-index surrogate. Subtract this event after the first-repeat bound. Equally spaced checkpoint detection runs on the actual deterministic map,not on a fresh-index walk.

Write p_d(x) for the d-step increment distribution and epsilon_d(x)=p_d(x)-1/N. Cyclic convolution is diagonalized by the characters,including complex eigenvalues for this possibly nonreversible chain. Fourier inversion and Parseval imply

sup_x |epsilon_d(x)|<=rho^d, and (sum_x |epsilon_d(x)|^2)^(1/2)<=rho^d.

These also hold at d=0. Set R=ceil(sqrt(N)), L=ceil(1031R/1000), T=ceil(log2(N))+16=271 and eta=2^-16. Since rho^T<=eta/N, for d>=T the return probability p_d(1) lies between(1-eta)/N and(1+eta)/N.

## Long-pair second moment

Among states Y_0,...,Y_L,count collisions C=sum I_(i,j) over i<j and j-i>=T. The number of such pairs is D=(L-T+1)(L-T+2)/2,so E C>=(1-eta)D/N. All equality events depend only on independent increments over their intervals; the starting group value does not affect the argument.

We bound ordered pairs of indicators in C^2. Identical intervals contribute at most(1+eta)D/N. For distinct intervals,include a baseline(1+3eta)/N^2 per ordered pair,whose sum is at most(1+3eta)D^2/N^2. The following bounds account for every excess above that baseline.

Disjoint increment intervals,including intervals touching at one endpoint,are independent and have joint probability<=((1+eta)/N)^2,already covered. Intervals sharing their first or last endpoint belong to the nested case below; identical intervals have already been separated.

For nested intervals with consecutive segment lengths a>=0,b>=T,c>=0 and a+c>0,the inner return and the outer residual return depend on disjoint increments. Their joint probability is p_b(1)*p_(a+c)(1),at most(1+eta)/N^2+(1+eta)*rho^(a+c)/N. There are at most(L+1)^2 choices for starting position and b. Sum over a,c with a+c>0 gives sum rho^(a+c)=1/(1-rho)^2-1=65535/1073676289. Both orders of the two indicators give total excess<=(131070/1073676289)(1+eta)(L+1)^2/N. This includes a=0 or c=0,so shared-endpoint nested intervals are counted.

For crossing intervals,the four ordered endpoints have positive successive gaps a,b,c with a+b>=T,b+c>=T. Independence of these three increment blocks gives joint probability sum_x p_a(x)*p_b(x^-1)*p_c(x). Expand each p as uniform+epsilon. Linear terms vanish. Cauchy-Schwarz on quadratic terms and sup/l2 bounds on the cubic term give

joint <=1/N^2+(rho^(a+b)+rho^(a+c)+rho^(b+c))/N+rho^(a+b+c).

Sum the rho^(a+c) term over positive a,c:the sum is(rho/(1-rho))^2=1/1073676289. Free starting position and b,plus the two indicator orders,give excess<=(2/1073676289)(L+1)^2/N. For the other two quadratic terms,sum_(a,b>=1,a+b>=T)rho^(a+b)<=2T*2^-T. The free remaining segment and start,with both orientations and both terms,give excess<=8T*eta*(L+1)^2/N^2.

For the cubic term,its constrained sum over a,b,c is below6*2^-T (dominate termwise by rho=1/2). Indeed,b>=T-1 contributes at most4*2^-T; for1<=b<=T-2,the conditions require a,c>=T-b and the geometric sum contributes below2*2^-T. After both indicator orders and the starting position,the contribution is below12eta*(L+1)/N. No interval pattern is omitted:ordered distinct intervals are disjoint,nested,or crossing.

Combining terms proves

E C^2 <= (1+3eta)D^2/N^2 +(1+eta)D/N +((131072+131070eta)/1073676289)(L+1)^2/N +8T*eta*(L+1)^2/N^2 +12eta*(L+1)/N.

This is a finite bound,not an assumption of independent collision indicators. In particular the short-return correlations of an adding walk are explicitly included.


For every integer c>=0, 1_{c>0}>=(3c-c^2)/2: the right side equals 0 at 0, 1 at 1 and 2, 0 at 3, and is negative thereafter. Therefore P(C>0)>=(3 E C-E C^2)/2>=(3m-U)/2. The exact production comparisons in the main argument give idealized repeat probability above 0.39015. Compression collision, bad table, bounded sampling and singular relation losses are all subtracted there, leaving more than 0.39. The coefficient-conditioning and checkpoint proofs in the main argument transfer this to the required verified G2 output. No independence of interval indicators or ideal random-mapping hypothesis is used.
