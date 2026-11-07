# Repeated-key dictionary with a weighted variance certificate

Use an explicitly initialized global dictionary of orbit keys,one node per distinct key. The serial chained implementation and all transfers are priced below. On a repeat,confirm all24 exactS words,then compare all four canonical beta words. Equal beta implies equal raw X=h*z^beta and returns CONTINUE without allocating or exponentiating. Different beta gives a useful nontrivial twist and terminates. Copy saved beta/k/key before any terminal bank switch. Payload index is only an address;A=1 for every sample and unused A padding is never read.

## Repeated-key lookup bound

Extend the independent table-index sequence through L for analysis. Let D count all time pairs with equal canonical orbit keys,including the identity's special orbit if it occurs. Equal parameter types occur with probability1/W^2 and always give the same point. For distinct types,the two raw values are independent uniform Fr,so orbit equality probability=(12(r-1)+1)/r^2<=12/r. Thus
E D<=binom(L,2)*(1/W^2+12/r).
Set Dcap=floor(L/2^20). Exact Markov arithmetic gives Pr[D>Dcap]<=E D/(Dcap+1)<2^-38. This is an analysis event,not an unpriced test or extra stopping loop.

Fix any extended key sequence with D<=Dcap,and fix NH seeds with no fingerprint collision. Different key labels have independent uniform128-bit buckets under the separately tagged BLAKE3 bucket function. Let n_u be each distinct key's multiplicity. Then sum n_u=L and sum n_u^2=L+2D. Let Z count co-bucket time pairs with DIFFERENT keys. It is the weighted sum over unordered distinct labels u,v with weight n_u*n_v. Distinct bucket-pair indicators have zero covariance,even when they share one key label. Therefore
E Z<=mu=L(L-1)/(2*2^128)<0.105L,
Var(Z)<=sum_(u<v)n_u^2*n_v^2/2^128<=(L+2Dcap)^2/(2*2^128).
The gap floor(L/8)-mu>0.020L. Exact Chebyshev arithmetic bounds Pr[Z>floor(L/8)|Dgood]<2^-100. This handles repeated keys explicitly;it does NOT treat repeated bucket queries as fresh.

Every nonmatching chain visit corresponds injectively to a different-key co-bucket pair of its query time and the node's original insertion time,hence contributes<=Z. Matching visits are at most the number of repeat-key queries,which is<=D (one earlier equal-key pair witnesses each). Thus total visits<=Z+D+1<=floor(L/4)+2 on both good events. Store that hard global128-bit visit budget;decrement/check BEFORE every visit and abort on exhaustion. It deterministically caps work on every execution,including failures. Dictionary cap failure is bounded by2^-38+2^-100,plus the separately charged NH collision/identity events. No independence between failures is assumed.

## Explicit finite memory and initialization

The storage-only model is material,not installed VM hardware. A noncacheable2^40-byte native window starts at2^45;shared code/data is below2^45,all native addresses<2^46. A97-bit bank selector uses two ordinary64-bit latches plus commit,with DSB/low/high/commit/DSB barriers. All address/index/control/data transfers are counted. There is no arithmetic/search/DMA oracle or flat136-bit native pointer.

Payload p>=1 maps q=p-1 to low31 local bits,offset local*512,bank q>>31<2^95;0 is nil. Head index i128 maps low36 local bits,offset local*16,bank2^95+(i>>36)<2^96. Table t0/1 with92-bit index j maps low30 local bits,offset local*640,bank2^96+t*2^62+(j>>30)<2^97. The table-bank low latch is<2^63,high latch region bit32. Two64-bit registers implement all128/92-bit arithmetic. Each selector including barriers/interface/RET fits32 instructions. Bank ranges are disjoint and no stale window pointer is used.

Zero ALL2^128 sixteen-byte heads by charged stores. Entire512-byte payload is initialized,then stores S192@0,beta32@336,k8@368,U40@384,next128@432;remaining bytes are zero padding. Entire640-byte table record is initialized,then rawGT576 and canonical v/w32@608 are written;unused A padding is never read. Raw960 shared records hold GT576,beta32@608,k8@944. All sampler/PAIR counters,Fr/field/Montgomery/trace scratch,NH1040-byte seed and hash buffers/constants are initialized cold. Complete code is in direct branch/ADRP range and paid initialized data.

On insertion,write the complete node with next=oldhead,complete noncacheable writes before switching back and publishing head. Heads are serial two-word operations;no atomic128 store is assumed. Indices increase,chains point to older nodes,so visits cannot loop. Fixed miss work<=512:three selectors96,full payload zeros32STP,S copies24LDP/24STP,beta4/k2/U6/next2 transfers,head/counter/allocation/control and frame/arguments/publication/RET leave over200 spare. Match paths do no insertion and include fullS/beta confirmation within512. Each visit<=128 including selector32,U comparisons,next128,hard budget and branches. Thus total<=512L+128(floor(L/4)+2)<=544L+256. Retain1024L;256 is cold paid. Ignored equal-beta repeats have no hidden exponentiation or reference cost.

Capacity:2^126 payloads*512=2^135;heads2^128*16=2^132;two tables each2^62 banks*2^40=2^103 including padding;shared/code<2^45. Total<2^136 bytes. All read bytes are initialized;unused capacity supplies no data. Memory is supplementary,not another metric. If this finite banked model is inadmissible,the claim must be rejected;no practical installation is asserted.
