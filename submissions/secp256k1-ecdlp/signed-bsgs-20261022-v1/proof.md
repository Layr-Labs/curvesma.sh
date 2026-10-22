# Signed baby-step/giant-step accounting candidate

This is an analytic classical construction for the protected full-production secp256k1 target. It is not an executed production attack, a new cryptanalytic breakthrough, or a qualified score. The submitted bound is time_log2=145, preprocessing_log2=145, memory_log2_bytes=146, success_probability=1/2, and nonuniform_advice_log2_bytes=0 (zero advice bytes). No heuristic is used. The only normalization is the protected C=1000 scalar AArch64 instructions per target-group equivalent. All arithmetic below is implemented with 64-bit limbs; no reference group operation is additionally charged.

## Algorithm and exact success

Let r,p,P be the protected order, modulus and generator. Set m=h=2^127 and d=2m-1. Validate the canonical input, including the curve equation and [r]Q=O. Build the m records (jP,j), 0<=j<m, by repeated affine addition. Records occupy 128 bytes: an infinity tag, four x limbs, four y limbs, four scalar limbs, and padding. Infinity has its own key; finite keys are x alone. Retain (m-1)P and compute dP=2((m-1)P)+P, D=-dP. Sort records by (infinity_tag,x) using bottom-up merge sort and a second buffer. No hash-table oracle is used.

For i=0,...,h-1, binary-search gamma=Q+iD. On a matching key inspect y to distinguish gamma=jP from gamma=-jP. For infinity use j=0. Form a=id+j or id-j modulo r, verify [a]P=Q, and output the canonical 32-byte scalar. Otherwise update gamma by D. Exhaustion returns failure. A match can never produce an incorrect scalar; verification is charged. There is no restart, randomness or amplification.

The baby keys are distinct: jP=±kP implies j=k in this range, since j+k<=2m-2<r. For fixed i the algorithm recovers exactly the residues in [id-(m-1),id+(m-1)]. Successive integer intervals touch without overlap because their centers differ by d=2m-1. Their union has exactly COV=h*d integers. COV<r, so reduction modulo r preserves cardinality. The interval wraps around zero but that changes no count. In particular

2*COV-r = 92138019644721193389046258963199934143 > 0.

Consequently success over the UNRESTRICTED uniform challenge scalar is COV/r>1/2>=0.39. The algorithm does not assume the scalar lies in a small interval; it deliberately fails on the complement. Fresh algorithm coins are irrelevant to this deterministic construction. The supplied target, generator and probability space are unchanged.

## Arithmetic implementation bounds

Use four 64-bit limbs per field element and a fifth carry limb when needed. Modular addition/subtraction costs at most A=256 instructions: fixed limb loads, ADD/ADCS or SUB/SBCS passes, compare/correct and stores. Each pass uses at most eight instructions per limb; six five-limb passes plus sixteen control instructions fit A. Scalar index arithmetic uses at most four limbs.

A field multiplication costs at most M=4096 instructions. Form the eight-limb product using sixteen limb-pair products, each with MUL and UMULH and fixed carry propagation. Reserve 128 instructions per partial product, including loads, stores and carries: 2048 total. Reduce using c=2^32+977 and the fold t -> low256(t)+c*high256(t). The first fold has high part at most c. The second yields less than 2^256+c^2. A third fold yields less than 2^256: if the second carried, its folded result is less than c^2+c<2^65; otherwise it is already below 2^256. Four fixed folds are safe. Each uses at most four limb-by-c products and fixed carry addition, budgeted at 384 instructions. Reserve another 512 instructions for initialization, conditional subtraction of p, normalization and output. Total 2048+4*384+512=4096. This is an analytic upper bound, not an assembly measurement.

Compute inversion for nonzero z by binary exponentiation z^(p-2). At most 512 field multiplications, 64 control instructions per exponent bit, and fixed setup fit I=2^22 instructions. Complete affine point addition handles infinity, equal points, opposite points and y=0 before division. For distinct x use slope=(v-y)/(u-x); for doubling use 3x^2/(2y). Then x'=slope^2-x-u and y'=slope*(x-x')-y. I+16M+32A+4096 < E=2^23 covers every branch and negation. Thus no incomplete addition or exceptional zero inverse is assumed. Scalar multiplication costs at most 512 E; input validation, step setup, at most one final scalar verification, serialization and fixed code preparation fit the fixed 2^40 reserve.

## Astronomical storage and paid routing

This construction requires exponential storage. It does not assume a single AArch64 processor can address it with a 64-bit pointer. Use two complete binary storage trees, each with m leaves and depth127, fewer than 4m nodes total. Each node is a scalar processor with at most 64 KiB local code, state, buffers and record storage, using ordinary local 64-bit pointers and three local word-transfer ports (parent and children). The leaf address is a 128-bit integer routed one bit per level. Only one operation is outstanding; processors wait for events rather than spin. All work of every processor, including transfer endpoints, is charged. No wall-clock speedup is claimed.

A complete record read or write including acknowledgement crosses at most256 edges in aggregate request/response traversal. Allow a 20-word packet (16 record words plus four command/address words) on each edge. Reserve sixteen instructions per transferred word covering BOTH sender and receiver loads/stores, port operations and address work, plus64 control instructions per edge. This is 256*(20*16+64)=98304 instructions; U=2^17 additionally covers the root/leaf fixed work, local comparisons and packet setup. Local neighbor communication is an explicit scalar-work model; distant memory is never a single free access. Latency is immaterial to this instruction-count score. The space and instruction budget includes routing state, processor state and code copies. Building the tree, copying its fixed code, clearing its storage and initializing paths costs at most 4m*2^20 instructions. Each node needs at most8192 word clears and8192 word code copies, including communication to its immediate children, within this allowance.

Each sort pass uses at most m key comparisons and m output records. With two cached merge heads, even eight routed record accesses per output is conservative. Additional 4096 instructions per record cover multi-limb counters, comparisons and control. There are127 passes; a 2^30 reserve covers pass setup. Binary search uses at most128 fetched records; each comparison and index update gets4096 extra instructions. A giant iteration additionally gets2^16 for scalar/index and sign handling. No coordinate recovery or square root is required because y is stored.

## Complete ledger

Let E=2^23, U=2^17, L=127 and m=h=2^127. All terms below are scalar instruction equivalents:

- fixed setup, input validation, successful verification and output: 2^40;
- all storage construction and code distribution: 4m*2^20;
- baby generation and writes: m*(E+U+4096);
- sorting: L*m*(8U+4096)+2^30;
- the entire giant loop, including unsuccessful lookups: h*(E+128*(U+4096)+2^16).

Their exact sum T is 29292028818269973496407165951123623938690646016. Total work is T/1000 target-group equivalents. Exact integer comparison gives T<1000*2^145; the diagnostic logarithm is about144.393416200472. The claim intentionally submits the conservative rational145/1. Preprocessing includes all storage construction, babies, sorting AND the full fixed reserve; its instruction sum is24909932386758706049095397310548220902765494272, below T, hence also below1000*2^145. Nothing is amortized over later targets.

Peak memory is below4m*2^16+2^24 bytes, including both routing trees and all their code/state, main-engine scratch, stacks, constants and output. This is below2^146 bytes. No random bits or nonuniform advice are retained. Public profile constants and algorithm design are the only excluded work, as allowed by the rules. Actual physical resources of this scale are unavailable; the claim is an analytic accounting control.

## Evidence, provenance and limitations

The machine-readable arithmetic evidence records exact integers and protected-file SHA256 hashes. Those calculations check the formulas; they do not certify their implementation costs. The experiment source is a bounded specification sanity check: exhaustive small cyclic groups, signed coverage and a deliberately broken sign variant; small tables on actual secp256k1; exact production-sized integer inequalities. Python event counts are not AArch64 instruction measurements. No candidate experiment was executed locally because AGENTS.md requires the shared llm-verifier runner. Hosted execution and review remain authoritative.

Source checkout: Layr-Labs/curvesma.sh at467130d867499e7bbca7479f0f72af0df115b282. The existing organizer baseline was consulted for the generic decomposition and submission format, but its old78037248 normalization and256-bit word assumptions are not reused. Current protected C=1000 and64-bit words govern this candidate. Existing submissions remain untouched.

Negation-aware BSGS is established prior work: Galbraith, Wang and Zhang, Computing Elliptic Curve Discrete Logarithms with Improved Baby-step Giant-step Algorithm, https://eprint.iacr.org/2015/605 . Only the elementary negation equivalence is used here; no paper performance constant is imported. SEC2 parameters: https://www.secg.org/sec2-v2.pdf section2.4.1. The full construction and cost justification are given here without relying on external documents. No novelty, best-known bound, practical compromise or independent qualification is claimed. Improving field arithmetic, storage organization and low-memory walks remains future work requiring fresh complete accounting.
