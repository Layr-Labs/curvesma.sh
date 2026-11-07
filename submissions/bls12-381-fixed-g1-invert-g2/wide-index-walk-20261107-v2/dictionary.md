# Explicit chained dictionary and exact index history

The dictionary stores every visited full Frobenius12 orbit,with exact192-byte canonical S key and NH U320 fingerprint. It terminates at the first U match;terminal cold code compares every exact S word and aborts on a false fingerprint. No unchecked U match becomes an answer. All source is static evidence,never host-executed. dictionary-index.S includes the complete dictionary body and both native selectors,not an unpriced allocator/search routine.

## Independent buckets and finite traversal cap

The standard104-byte JUMP digest has independent disjoint fields:j in low94 bits,bucket in next128 bits. Group transitions inspect only j. Conditional on all ideal raw group values and indices through the first orbit repeat,distinct keys have independent uniform buckets under H,apart from the separately charged NH collision loss. Extend with fresh labels beyond first repeat only for analysis. The final repeated query reuses its saved bucket;it is not treated as a fresh draw.

Let m=2^128,n=K+1 and C count co-bucket pairs in this auxiliary distinct-key prefix. Its mean mu=n(n-1)/(2m)<0.101K. Disjoint pair indicators are independent;overlapping pairs have joint probability1/m^2 and also zero covariance. Thus Var(C)<=mu<K. Retain the conservative old bound<K. At threshold floor(K/8),the gap from mu exceeds0.024K. Chebyshev gives Pr[C>floor(K/8)]<1/(0.024^2*K)<2^-100 since K>2^124.

Distinct insertions traverse at most C old records total. The terminal repeated query traverses at most C additional nonmatching records plus its match. Hence all visits<=2C+1,including the final query. On C<=floor(K/8),this is<=floor(K/4)+1. A global two-limb counter initialized floor(K/4)+2 decrements before every old-node access;fail before underflow. Every path has at most floor(K/4)+2 visits,even outside H. Its extra cap failure is<2^-100. Links point only to earlier allocated indices,so chains are acyclic.

## Initialized layout and A recovered from the index

Initialize all2^128 sixteen-byte heads to zero in charged cold work. Logical payload index0 is nil. A separately initialized32-byte dictionary control at context pointer24 holds next_index128 initially1,and remaining_visits128 as above. Context24..31 is disjoint from label0..15,reservedzero16..23,K32..47,seed pointer48,key64..255,Ublock256..319,prefix320..351,mutable digest352..383,bucket496..511. The fixed window base x27 is initialized by the cold parent.

Each512-byte payload is entirely zeroed by32 explicit STP before publication. Store S0..191,B336..367,k368..375,U384..423,next432..447. Bytes192..335,376..383,424..431,448..511 remain initialized zero padding. In particular no A is stored at the previous304..335 slot. Every byte is initialized;no uninitialized advice,free zeroing or bulk transfer is used. The complete finite record capacity2^126 fits all K+1 states.

At initial time0 allocate index1. Until the first U match,each completed transition visits a new orbit and allocates exactly once;identity,fingerprint/cap failure aborts rather than skipping and continuing. Therefore every allocated payload p has exact original time s=p-1 and exact first history coefficient A_s=s+1=p. This is an integer identity of the serial allocator,not a statistical assumption or a narrow cryptographic input. p<2^126<r,so convert it to a full32-byte Fr integer with upper128 bits zero for terminal h^p. The current coefficient is A_t=t+1 from the preserved two-limb accepted counter. No raw-record A or payload A is read.

The raw state keeps full576-byte GT,B32 at608,k at944 in initialized960-byte buffers;former A576..607 is unused initialized padding. scalar-and-copy.S updates B modulo full r in26 instructions and no longer writes A. All protected group/field/order sizes remain unchanged. Initial B=beta0 and jumps have fixed u1 plus full v_j.

## Exact ABI and native selectors

Call fapi_trace_dictionary(main_context,current_rawrecord),with cold-initialized x27=window base. It preserves x19..x28 and original LR via its96-byte aligned frame. Return0 means inserted,1 means U match,2 means cap failure. On1,x1/x2 are the matched logical p low/high and x3 is its native payload pointer. The selector remains on that payload. The cold terminal handler immediately copies p,B,k and all exact-key confirmation data to initialized shared scratch before any bank switch/call;it uses A_s=p,A_t=t+1. The true three-twist reconstruction and exact h^a=z/G2/pairing verification remain cold paid. It must not read the unused A fields.

The entry frame saves x19..x26 and LR;its slot80 stores the original head128. x19=context,x20=current raw,x21=control,x22=U40,x23=S192,x24/x25=cursor index,x26=node pointer. x27/window and x28 are unchanged. Helpers touch only caller scratch,so every live value survives their BL. Return copies matched index/pointer before restoring callee state. x18 is untouched;no hidden context spill or allocator exists.

The head selector has13 instructions including RET. It forms local36-bit offset<<4,bank low/high of i>>36,sets highword bit31 for region2^95,loads its initialized MMIO pointer,executes both DSB and all three low/high/commit stores,and returns native window pointer. Payload selector has12 instructions including RET:local31-bit offset<<9,bank p>>31 in region0,the same pointer/latches/barriers and native result. All paired offsets are aligned and<=496. A97-bit selector is latched by two ordinary64-bit stores plus commit,not an atomic128-bit store. This explicit finite storage-only model remains a material admissibility condition,not installed/practical VM hardware or a compute/search oracle.

The first DSB of each switch completes prior noncacheable writes. Payload is fully initialized and its next=oldhead set before switching back to publish its head. After return,the next jump selector's DSB completes head publication before selecting a jump. Saved old head lives in shared stack,not a stale window pointer. Every logical index is128 bits;all native addresses stay<2^46 and pointers64 bits.

## Complete scalar count,including every branch/interface

Direct static basic-block counts from dictionary-index.S:

|block|instructions excluding called leaf bodies|
|---|---:|
|entry before chain test|15|
|one full nonmatching chain visit|31|
|match exit before shared return|5|
|cap exit before shared return|2|
|insertion body including32 zero stores and24 key-copy instructions|86|
|shared restore/RET|6|
|head selector complete leaf|13|
|payload selector complete leaf|12|

The fixed empty-chain miss path costs entry15+head13+empty test2+insert86+payload12+head13+return6=147. This includes both selectors,allocation counter increment,full512 zeroing,all S/B/k/U/next loads/stores,head publication,frame,call arguments and RET. One complete nonmatching visit costs31+payload12=43. For a match,the last visit skips its next-load/loop-back2,then adds match5+return6;total is<=147+43*visits. Cap exits also fit that bound because their partial guard and return are cheaper than a complete visit. No path-dependent cost is omitted.

Retain fixed160 per dictionary call (13 spare) and48 per visit (5 spare). There are K ordinary calls plus the initial cold call;hard visits<=floor(K/4)+2. Therefore all dictionary work<=160K+12K+256=172K+256. The initial160 and extra96 traversal work are cold;complete coldD=540*2^119 already includes them and full head initialization. Once-only exact S confirmation/raw reconstruction/output is separately cold paid. No orientation average,free interface or probabilistic work bound is used.

## Whole computation and status

Ordinary scalar subtotal=1264canonical+460NH+776fullBLAKE3+26Bupdate+72controller+160dictionaryfixed+12cappedvisits=2770. Retain2800,with30 additional spare instructions. One protected reference C20000 per accepted transition and complete per-target coldD give work<=(22800K+540*2^119)/20000. Exact certificate(22800K+540*2^119)^2500<20000^2500*2^314543 gives proposed314543/2500=125.8172,above diagnostic125.81719586241536. Success/caps/kernel/interval proof are unchanged;only redundant storage and implementation allowances change. All six production tracks,one computation metric,full p/r,threshold0.39,reference prices,zero advice and memory<2^136 remain preserved. H and finite banked-memory admissibility remain material;this is conditional analytic evidence,not official scoring,qualification,promotion,practical attack execution or a fixed-hash theorem.
