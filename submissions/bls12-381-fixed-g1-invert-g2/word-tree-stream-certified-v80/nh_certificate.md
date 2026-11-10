# Composed NH128 and affine fingerprint, full160-bit collision proof

## Self-contained wrapping-NH bound (a conservative2^-63)

The NH construction is credited to Black, Halevi, Krawczyk, Krovetz and Rogaway, UMAC (1999), https://www.cs.ucdavis.edu/~rogaway/papers/umac-full.pdf . The bound USED HERE is proved below without relying on the paper's sharper theorem.

Put N=2^64. For unequal W,W', select a differing pair and exchange its two positions if needed so the second words differ. Fix every key except its first k. Then c=(W_second+K_second) modN and c'=(W'_second+K_second) modN are distinct. Let x=(W_first+k) modN,d=(W'_first-W_first) modN; x is uniform and x'=(x+d) modN. All other pairs contribute a fixed difference y moduloN². A collision forces
(c-c')x = y+c'd-Nc't moduloN²,
where t=0 if x<N-d, otherwise1. Within either region, two solutions would give nonzero (c-c')(x1-x2) divisible byN², impossible since its absolute value is belowN². Thus at most two of N keys collide. Averaging the fixed keys proves2/N=2^-63, including d=0. This uses unsigned wrapping additions and the full128-bit accumulated output.

This supplies every NH premise used by the following affine composition. The two-region bound is intentionally looser than the paper's theorem, avoiding any need for an unproved cross-region exclusion. It is a rigorous finite bound for all distinct full192-byte strings, independent of empirical checks.

Interpret full S192 as24 unsigned little-endian64-bit words W. For each of three INDEPENDENT rows sample24 independent uniform64-bit keys K. Compute H=sum over i0..11 of ((W2i+K2i) mod2^64)*((W2i+1+K2i+1) mod2^64), reduced modulo2^128. Do NOT truncate H. Encode H injectively as three chunks of60,60,8 bits. Each chunk is strictly less than q=2^61-1, so no different128-bit H values acquire equal chunk triples modulo q. Independently sample four uniform Fq values b,a0,a1,a2 for each row and compute Y=b+a0*C0+a1*C1+a2*C2 modulo q. For two unequal chunk triples, the affine map from its four coefficients to (Y,Y') has rank2, so the pair is exactly uniform Fq². Keep low54,54,52 bits of the three Y values.

For b output bits, a residue class has at most ceil(q/2^b)=2^(61-b) field elements. Thus, CONDITIONAL on distinct NH results, the collision probability of low bits is at most2^(61-b)/q. Conditional on equal NH results it is1, whose NH-key probability is at most2^-63. A row's total match probability is at most2^-63+2^(61-b)/q. The three rows' entire key sets are mutually independent, even for the same fixed S,S', so their match events are independent. Therefore delta=(2^-63+128/q)^2*(2^-63+512/q). This bound applies to ANY distinct192-byte keys, conditional on the complete generated sequences; no independence between different pairs of records is assumed. Equal exact keys always get equal fingerprints. The full S192 still decides equality at confirmation.

## Native scalar source and exact layout

nh_fingerprint.S uses ONLY x0..17 and flags, with no SP,x18,callee19..30 or LR changes and no BL. Arguments are x0 output32,x1 S192,x2 keys672. The72 raw keys are interleaved by input pair: each48-byte group contains the three rows'16-byte pairs;12 groups occupy576 bytes. The three affine rows (b,a0,a1,a2) occupy32 bytes each at576,608,640. This is a deterministic permutation of the same independent keys, not shared randomness between rows. All field coefficients are canonical less than q. Keys are disjoint from output. Output and S may be exactly aliased or disjoint, never partially overlapping.

First zero six accumulator registers (three pairs). Twelve pair rounds each use one LDP of S and, for EACH row, LDP key, two wrapping ADDs, MUL, UMULH, ADDS low and ADC high. Thus each round costs1+3*7=22, total264, plus6 zero instructions. Wrapping high ADC implements exactly modulo2^128. Input words are read once, and ALL192 input bytes have been read before any output write.

For each row, extract its three chunks using AND,EXTR,AND,LSR (4). Initialize affine128-bit accumulator from b and zero high (2). Three LDR/MUL/UMULH/ADDS/ADC groups cost15. Its sum is below3*q*2^60+q<2^123, so it fits128 bits with no lost carry. The existing full-width Mersenne reduction logic costs11 INCLUDING output mask: split into61/61/6-bit parts, fold, one correction, then low-bit mask. The folded intermediate after the second fold is less than2q, so one conditional subtraction is sufficient. Total per row32. Retain first result in x1 after S reads finish, then x4, then x5; pending NH accumulators remain live until their row is consumed. The scratch x10..17 never overlaps an unconsumed accumulator.

One global MOV creates q (logical-immediate alias). Counts:1 modulus MOV+1 key-base ADD+6 accumulator zeros+264 rounds+3*32 affine+6 pack/store instructions+RET=375 throughRET. The six final instructions concatenate54/54/52 bits and write the full32-byte output with its upper96 bits zero. Four caller instructions give379, retained384 (was496). No hidden vector instruction, multiplication oracle, extra argument loads or omitted memory writes.

## Complete probability and cost

The exact four-family group event counts m,L,s,M,G, all success/failure caps, full production fields and fixed C20000 remain unchanged. Substitute ONLY the proved delta above in caploss. Exact finite-rational evaluation gives final probability0.3900268131920969, strictly above3900031/10000000 and0.39. Decimal values are diagnostics; the displayed formulas determine the exact rational bound. All false fingerprints, ignored baby pairs, giant identities, capped coins and failures are included.

Generation=2816+832+384+32+704=4768. Sorting1360+scan128+confirmation1+spare256 gives base6513. Network1369 gives hot7882, retained7968. Complete N=7968L+2^124 satisfies the EXACT integer inequality N^1000<20000^1000*2^124753, giving proposed time_log2=124753/1000=124.753. Cold124/peak138/advice0 still pay all actual tables, processors, random source, keys, transports, failures and terminal verifications. The score is a proposed analytic bound, not observed execution, hardware qualification, acceptance or promotion. Shared llm-verifier alone runs/reviews/scores.

The modulus q is prime: the Lucas-Lehmer recurrence v0=4,v_(i+1)=v_i²-2 modulo(2^61-1) has v59=0. This is a finite prime certificate and is separate from any candidate execution claim.

The fingerprint and its complete finite probability bound are unchanged by the fused recurrence transform. The full hot sum is7882, retained7968 with86 instructions above the displayed sum. This changes no fixed reference price C20000, success threshold or group parameter. Exact integer comparison for the full displayed L proves the124753/1000 bound.

The key-base ADD and all immediate encodability constraints are detailed in stream_repair.md. The shared audit enforces the pair-load ranges and expects375, while the retained384 budget is unchanged.
