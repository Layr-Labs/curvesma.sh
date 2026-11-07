# Chained global dictionary with a hard total traversal cap

The dictionary stores every distinct visited orbit. It changes storage only; the first repeated orbit, transition labels, nonbacktracking law and relation extraction remain those in proof.md. It requires no sorted tree, random key-rank assumption or expected-only work accounting.

## Exact hash and layout

DICT is ASCII FAPI-NB-DICT-v1 followed by one zero byte (16 bytes total). Compute B3(DICT||seed32||zero16||U40) with standard full seven rounds: cache prefix64 len64 flags1 cold, final40 len40 flags10 counter0. The first two little-endian64-bit digest words give an unsigned128-bit bucket index. This additional domain is explicitly included in H. Actual BLAKE3 is computed and charged. It is independent of JUMP/PHAS/retry domains inside that same disclosed bridge, not a second heuristic.

Allocate m=2^128 heads,16 bytes each, and initialize EVERY head to zero during this target's charged cold phase. A head is a128-bit logical payload index,0=nil. Payload indices begin at1. A512-byte payload contains canonical b at0..287,constantzero288..295,zero alignment296..303,A at304..335,B at336..367,k at368..375,zero376..383,U at384..423,zero424..431,next logical128-bit index at432..447,andzero448..511. All paired pointers are16-byte aligned. Every allocated payload is fully initialized before its head is published.

For a query compute its bucket,load head to shared scratch,and traverse the linked payload indices. At each visit FIRST charge the traversal counter,then select the payload bank and compare its five U words to shared current U. If unequal,read next and continue. If equal,compare all36 canonical b words: unequal means fingerprint failure;equal means immediate exact orbit match and terminal recovery. This exact comparison and all finishing work occur at most once and are charged cold. Never accept a truncated key. On a miss allocate a payload,copy b,A,B,k,U and old head as next,then select the head bank and store its index. All buffers needed across bank changes are shared. No memory search oracle or implicit sort is used.

Initialize a128-bit remaining traversal counter floor(K/2)+2. Each visited old payload consumes one before its loads;on zero fail. The initial record is inserted cold. At most K insertions follow. Counter checks and failed dispatch are charged. Thus at most floor(K/2)+2 visits occur on EVERY path,including outside H. At128 scalar instructions per visit this is at most64K+256;the fixed256 is cold. This is a deterministic hard work bound,not an assumed mean chain length.

## Finite failure proof

Condition on the complete ideal group trajectory through its FIRST orbit repeat,and on the tabulation tables. The trajectory uses JUMP/PHAS/retry,not DICT. Except for the separately subtracted fingerprint failure,distinct orbit keys have distinct U. Consequently the DICT values at their distinct inputs are independent uniform buckets inside H. Extend the list to n=K+1 with fresh independent uniform dummy buckets;this can only enlarge the pair count. The extension is a probability argument,not an executed hash or stored record.

Let C=sum_(i<j)1[H_i=H_j] for n independent uniform m-bucket values. All pair indicators are pairwise independent:disjoint pairs are independent,and overlapping pairs have joint probability1/m². Therefore
 E C=mu=n(n-1)/(2m), Var C=mu(1-1/m)<=mu.
Set a=floor(K/4). Exact protected arithmetic gives mu/K<0.101 and a/K-mu/K>0.149. Chebyshev yields
 Pr[C>a]<=mu/(a-mu)²<2^-100.
No independence among the whole collection of pair events is asserted.

Before the final query,each distinct insertion traverses precisely the earlier members of its bucket. Their total is at most C. If the first repeated-key query occurs,it traverses at most all distinct prior members of that bucket. Their count excluding the identical key is at most C (the pairs with that key),plus one for its exact matching record. Thus total visits are at most2C+1. On C<=a,this is at most floor(K/2)+1,below the hard cap. If an earlier scientific failure or work cap stops the algorithm,there are fewer visits. This also covers no repeat before K. The extra failure probability is therefore below2^-100;it is subtracted separately. We do not incorrectly treat a repeated query's dictionary hash as fresh.

## Complete scalar interfaces

Per ordinary transition,head lookup/publication is bounded by192 instructions: bucket extraction/shared copies<=16;two head bank selections<=64;head load/store<=4;128-bit address/allocation/control/context<=64;remaining44 slack. The extra hash costs797,including complete full output,stack preservation and return. Its prefix-copy/arguments/caller work and traversal-cap setup/dispatch fit the enlarged controller848 (prior784 plus64). Hash block padding stays unchanged zero.

Each visited payload is bounded by128:selector32;current/saved U loads<=6;five XOR/four ORR/conditional branch10;next-pointer load1;two-limb counter load/subtract/store/guard<=10;logical cursor/shared copies/control/interface<=32;remaining37 slack. No arbitrary-width primitive is used. All128-bit values use two ordinary64-bit limbs. Exact b confirmation is a once-only terminal comparison,not omitted per candidate;equal U with unequal b terminates immediately as a failure.

Payload initialization/write is224 instructions,including full512-byte zeroing33,bank selector32,18 b load/store pairs36,coefficient/phase/U/next transfers<=32,logical allocation/context/interfaces<=64,and27 slack. Field transfers are scalar instructions. No per-step large-record copy or uncharged allocator exists. Fixed arrays and monotone allocation need no free list. Counter/allocation/bucket-head publication use ordinary paired loads/stores with aligned local offsets.

For head-array initialization,use the explicit zero8192 helper:512 STP xzr,xzr plus RET. Each8192-byte block includes <=24 caller argument,address/loop/context instructions. There are2^119 such blocks,so their work is537*2^119<2^129. Bank selection occurs once per2^40-byte window,at most32*2^92=2^97 instructions. Index arithmetic and bank-loop setup remain below2^100. Add the complete previous cold work<2^120;the sum is below2^130. This counts every zero store and control instruction without a pre-zeroed RAM assumption. The initial payload,DICT prefix,cap and first hash are included. All dictionary work is CPU work;there is no DMA initialization or free memory oracle.
