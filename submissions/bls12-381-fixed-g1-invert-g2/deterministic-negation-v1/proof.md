# Centered deterministic inversion with exact inverse-orbit keys

Proposed time_log2=142 (rational 142/1), success probability 1, peak memory_log2_bytes=145, preprocessing_log2=142 included cold in total work, zero nonuniform advice bytes. This is a full-production analytic accounting control under C=20,000. It is not a practical production attack or a claim to beat the recorded heuristic comparator 126.425. No auxiliary powers are supplied or assumed.

## Material sources

The supplied-P reduction and sorted deterministic BSGS correctness foundation reuse public submission 487d6285-ba02-4b3a-bb79-e87d17a61fd2, declared as a dependency. This revision also builds directly on our submitted f00694ca-67df-4231-b803-8825a49ef896 (routing-v2, proposed 144), whose finite routing model and cold ledger are retained. At preparation it was still validating rather than published; the CLI skill permits a published dependency to be added after both endpoints publish. The direct use is disclosed now and will be formally declared when eligible. No acceptance, official score or promotion is inherited. The current change is an exact inverse-orbit key and centered coverage proof that reduces baby-table storage and sorting while preserving deterministic success.

The historical fixed-Q repository proof was inspected but its old word convention and prices are not used. Public note e750ade was inspected as a stronger comparator; its twelve-orbit trace compression, finite-rank random sampling and hash bridge are not used. The two-element inversion symmetry here is derived directly below, without its trace formulas. The reference tower definition is the pinned blst source used by the organizer: https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fp12_tower.c . All required algebra is included here; no external lookup is necessary to reconstruct the claim.

## Exact key and local cost

Use the native reference tower Fp2=Fp[u]/(u^2+1), Fp6=Fp2[v]/(v^3-u-1), Fp12=Fp6[w]/(w^2-v). Write X=A+B*w with A,B in Fp6. The nontrivial Fp6 automorphism sends X to bar(X)=A-B*w. It equals p^6-Frobenius, since this is a quadratic field extension. The protected integers satisfy p^6=-1 modulo r (checkable by six modular multiplications), so every order-r element satisfies bar(X)=X^(-1). This includes identity. Neither a subgroup reduction nor a replacement challenge point is being made.

All tower coordinates use canonical six-limb Montgomery residues. The scale and field basis are fixed and invertible. Define K(X) to be the lexicographic minimum of the full 72-word arrays for X and bar(X). Then K(X)=K(Y) if and only if Y=X or Y=X^(-1). This follows because the two-element sets are either disjoint or identical; selecting one member cannot identify distinct such sets. No hashing, lossy fingerprint or heuristic collision claim is present.

Compute bar(X) by negating exactly the six base-field coordinates in B, preserving A. For one canonical Montgomery residue t, use p-t if t is nonzero and zero otherwise. Six-limb subtraction, zero testing and conditional selection fit below 128 scalar instructions per coefficient, including loads and stores. Six such negations, a full 72-word comparison and full record movement use below 8192 instructions. Charge a full F=2^20 instructions for this local key computation; it needs no dense multiplication or inverse exponentiation. Keep the raw X separately in local scratch so canonicalization does not change the iterated walk value. Copies and all zero cases are paid. Fixed conversion of the public inputs to/from the native tower is setup, not a per-sample operation.

## Centered BSGS algorithm

Validate P,z as in the protected profile, including nonidentity P and subgroup checks. The identity z=1 returns canonical G2 infinity after verification. Otherwise compute h=pairing(Q0,P), using the pinned G2-first library call. Nondegeneracy makes h a generator of GT for every allowed supplied P.

Let m=2^126, S=2m-1=2^127-1, H=ceil(r/S). The exact protected r satisfies r<2^128*S, so H<=2^128. Generate raw h^j for j=0,...,m-1 by repeated dense multiplication. Store 80-word (640-byte) initialized records: 72 words K(h^j), two words j, six zero padding words. Do not replace raw h^j by its key before multiplying by h. Sort all records by their complete keys with bottom-up merge sort, two full arrays and 126 passes. No two baby keys coincide: if h^j=h^(-k), then j+k=0 modulo r; with 0<=j,k<m<r/2 this requires j=k=0.

Let centers c_i=m-1+iS. Compute gamma=z*h^(-(m-1)) and D=h^(-S), using binary exponentiation modulo r in setup. For i=0,...,H-1, compute K(gamma) locally and binary-search the sorted array. On a match j, try x=(c_i+j) mod r and x=(c_i-j) mod r, verifying h^x=z. At least one is correct. For a verified scalar construct Q=[x]Q0, validate G2 membership and verify the full pairing(Q,P)=z before serializing the exact G2 output. The identity gamma matches j=0. If no key matches, update gamma=gamma*D and continue. All at most two terminal scalar/pairing attempts are charged in setup/finish.

## Deterministic success on the complete input distribution

For every integer i, the exponents c_i-j and c_i+j for 0<=j<m fill every integer in the interval [iS,(i+1)S-1]. The intervals are adjacent with no holes: the low endpoint is c_i-(m-1)=iS; the high endpoint is c_i+(m-1)=iS+2m-2=(i+1)S-1. Since H*S>=r, these intervals cover every a in [0,r). For z=h^a, at its interval gamma=h^(a-c_i)=h^(j) or h^(-j) for some baby j, so the exact inverse-orbit key is found. A match always implies one of the two tested exponents is a modulo r, even at a wraparound boundary. Thus h^x=z and bilinearity gives pairing([x]Q0,P)=z for the actual supplied P. Success is 1, not an estimate derived from toy experiments or a random walk model. No random coins are needed.

## Finite storage and transfer convention

The storage requirement is enormous and cannot be allocated on this VM. The challenge permits analytic arguments for exponentially large attacks. This is a finite abstract classical computation, not a claim that the current physical machine has the attack's storage. In particular, a 64-bit register is never assumed to hold a 140-bit address.

Each logical address has three 64-bit limbs. Realize storage as a binary tree of routing controllers and leaf banks. A leaf bank holds at most 2^30 bytes with ordinary 64-bit offsets. Both arrays together contain 1280*2^127 <2^138 bytes; fewer than 2^108 leaf banks suffice after rounding. Choose a padded complete tree with 2^108 leaves and depth 108. Reserve code, controller state and messages separately below. Controllers route one bounded record or comparison request down and back along its three-limb address. Only a constant number of requests need be in flight, as the whole algorithm can be serialized. Work is summed over all controllers; distributed storage does not give free computation or a parallel speedup.

For a record transfer or comparison, charge 2^20 scalar instructions. An elementary hop carries at most a record of 80 words plus a header of 16 words. Per word allow 16 scalar instructions for loading, storing, addressing, acknowledgments and link/control work, i.e. 1536 per hop. A comparison may fetch two records through the root and return a result: at most four depth-108 paths, hence 432 hops. Even charging 1600 instructions per hop gives 691200 instructions. Add 10000 for endpoint comparison/copy, three-limb arithmetic, loop control and message initialization; the result is below 2^20=1048576. Every transferred word, including padding, is charged. Ordinary resident arithmetic operates in a small local bank, never across unpriced global addresses. A control message to create or initialize a bank is charged separately as preprocessing.

This is a conventional finite word/message-storage model, not a measured device implementation. Its material limitation is physical realizability at this scale. The argument assumes links transmit a bounded word with bounded scalar work as above, with no additional distance-dependent energy or latency objective. Latency is not the score. If the normative instruction model disallows any such banked storage, this candidate should be rejected rather than silently treating addresses as single words. The claimed bound includes the concrete routing depth instead of asserting free memory switches.

## Complete work and exact cap certificate

Use the same F=2^20 instruction bound for a routed record comparison or transfer, in the explicit storage model above. Generation costs m dense reference operations plus 4F*m instructions. This includes canonicalization, writing the initialized record, preserving/updating the raw element and all local generation bookkeeping. Each sort pass costs at most 4F*m: m comparisons, at most 2m record moves, plus F*m for merge/run/address control. There are 126 passes. Round generation plus sorting up to 512F*m instructions (4+4*126=508<512).

A complete giant iteration has at most 128 comparisons, record/index retrieval, fixed-size address arithmetic, local key construction and loop updates. Charge 134F instructions, including the extra local key and raw gamma preservation. Its raw update costs one dense reference operation. Hence

 Wmain <= m+H + (512*m+134*H)*2^20/20000.

Every multiplication and every record operation, including failed lookups and all losing key orientations, is included. From m=2^126 and H<=2^128,

 Wmain <= 5*2^126 + 1048*2^146/20000 <2^142.

The setup/initialization allowance below is added before testing the actual cap; it is not inferred merely from two separate terms being <2^142. In exact integers the certificate is

 20000*(m+H)+(512*m+134*H)*2^20+2^140+2^129+2^100 < 20000*2^142,

where H=ceil(r/(2^127-1)). This can be checked by integer arithmetic alone with the public protected r. It certifies the total claimed cap, not an experimentally measured score. The submitted bound is 142 unchanged; ordinary review is not asked to automatically optimize it. Sorting and all table construction are paid cold for this one supplied P,z. No future target shares the tables.

## Setup, representation, output, bank construction

Charge 2^100 scalar instructions for all fixed-size field conversions, canonical encoding validation, P and z subgroup tests, the initial pairing h, stride inversion, final scalar multiplication in G2, output subgroup test, final pairing and serialization, as well as code and constant construction. This extremely conservative bound is justified by elementary bit arithmetic: at most 381-bit base residues, at most 4572-bit field exponents, 12-coefficient polynomial arithmetic and standard projective or affine curve formulas. A 384-bit modular multiplication implemented by shift/add with one bit per byte uses fewer than 2^32 scalar instructions; exponentiation uses fewer than 2^14 such multiplications for the exponents here. Polynomial convolution uses 144 base products and fewer than 256 additions per extension multiplication. G2 uses quadratic arithmetic with u^2=-1. An affine inverse uses binary powering by p-2; a scalar multiplication has at most 510 additions/doublings. A Miller pairing has fewer than 128 loop iterations, each expressible with fewer than 2^16 base/extension operations, and final exponentiation fewer than 2^14 extension multiplications. Even multiplying all these loose upper bounds is below 2^100. Exceptional curve additions test infinity, opposite points and zero doubling denominator before division. These are fixed bounded computations, not exponentially repeated setup hidden under an asymptotic symbol.

For memory construction, initialize both complete arrays (at most 2^138 bytes) by word writes with at most 32 instructions per word, costing <2^140 instructions. Each of the fewer than 2^109 routing tree nodes receives at most 2^20 instructions of initialization, costing <2^129. Code replication at tree nodes is included in the storage/initialization allowance, not uncharged advice. Add all three budgets: (2^140+2^129+2^100)/20000 <2^127 reference equivalents. Together with Wmain<2^142 this is strictly below 2^142, as certified by the exact integer inequality below. The submitted score is the upper cap 142 unchanged; it is not a decimal estimate of an exact optimum. Preprocessing including table generation and sorting is a subset of the fully charged total, so preprocessing_log2=142 is valid. It is never amortized across later targets.

Peak storage is bounded by the two arrays <2^138 bytes, routing nodes fewer than 2^109 each with <=2^20 bytes of code, message and scratch capacity (<2^129 bytes), and <2^30 bytes of coordinator code, arithmetic scratch and public constants. Their sum is below 2^139 and therefore below the declared 2^145 bytes. There are no searched embedded constants, precomputed answers or nonuniform advice. Randomness use is zero.

## Limits and validation status

No production-size attack or storage allocation has run. No candidate file is executed on the host or with credentials. No experiments are declared: the correctness and finite cost certificate are analytic. The ordinary account cannot access the Docker daemon socket, so bounded local experiments were dropped. This is not evidence that proof-only review requires production attack hardware. The finite bank/message computation remains a material modeling assumption; reject it if incompatible with the normative word-instruction model. It is not a free-address or free-transfer oracle. Local intake validates structure only; hosted llm-verifier review and manual promotion are still pending. Acceptance is not formal certification, novelty or a practical cryptographic break.
