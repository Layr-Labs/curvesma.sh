# Complete current recurrence: shared vectors, one modulus load, private map

All residues are canonical full protected Fp381 Montgomery values Rm=2^384. Write c=A+B*t,X=a+b*t,Y=d+e*t in Fp4/Fp2,t²=xi=1+u,u²=-1. The exact recurrence factor is
A(a-d)+Bxi(b+e)+t[A(b-e)+B(a+d)].
Set V0=a-d,V1=xi*(b+e),V2=b-e,V3=a+d. Compute A*V0+B*V1+t*(A*V2+B*V3),plus the old third history. Only A and B coefficient tables are needed;the complete GT recurrence and all full field coordinates are unchanged.

## Actual fused wrapper and the public register contract

recurrence_fp2_step is the ONLY public PCS entry. Its528-byte frame saves19..29/LR before ANY clobber,and restores ALL after the map. SP/x18 are preserved. Four BOTH transforms construct a±d,b±e in initialized384 bytes atSP144..527. The fifth BOTH uses inputSP240/SP288 and writes differenceSP240/sumSP288,exactlyxi*(b+e). Every input is read before any output store;SUM25..30 remains live after MINUS stores. This exact crossed output alias is safe and touches neither raw histories nor other transformed components. ALL wrapper frame accesses useSP because BOTH clobbers29/LR.

Load the full six p limbs ONCE into19..24 before the first BOTH. No BOTH writes19..24 and no intervening setup touches them. Every full-width CSEL/reduction therefore sees the same correct p. This does not require the external caller to preload anything. Each BOTH body is48 instructions after removing three repeated p loads. Five bodies plus50 actual setup/prologue/epilogue instructions give wrapperown290. Hoisting alone removed16 actual instructions from the preceding306-own wrapper. Full fields/output/carry/alias behavior are identical.

## All ordinary coefficient and fixed-s tables are constructed

For C=A,B,each vector entry is(d*Cr*2^start/Rm modp,d*Ci*2^start/Rm modp),96 bytes. Real input digits add(first,second) to(real,imag);imaginary digits SUBTRACTsecond fromreal and ADDfirst toimag. This is exactly C*(xr+xi*u) foru²=-1. Both output row pairs reuse these same A/B tables. No negated duplicate table or oracle is supplied.

ALL381 significant input bits have widths[35,35,35,35,35,35,35,34,34,34,34],starts0,35,70,105,140,175,210,245,279,313,347. One vector column is96*(7*2^35+4*2^34)=27TiB;two are54TiB. EVERY entry is constructed by actual full-p pinned seed products and existing fill_pair75(n-1)+17<=96n canonical additions,including all zero entries and headers. Full64-bit counters handle2^35/34 sizes. This is ordinary constructed memory,not virtual production arithmetic for free.

Vectors occupy8..62TiB;unchanged fixed-s13.5TiB is moved by its existing base argument to64..77.5TiB. Parent4TiB/private<1TiB/childports80..128TiB remain. All ranges are disjoint and pointers strictly below2^47. Other routing maps and DATA identifiers are unchanged. All new table pages/seeds/constants/MMU work is charged;bothworkers' vector/fixed fills remain<2^50 instructions inside completecold124. Bothworkers'tables<2^48 bytes/globalpeak138;no physical VM memory qualification is claimed.

## Signed accumulators stay nonnegative and below45p

Each real accumulator starts OLDthirdreal+22p. The ordinary seven-limb22p is already constructed in the50-slot header at528+64*22=1936. Exactly22 canonical negative table terms are subtracted. After k<=22 such terms,the accumulator is at least seed+(22-k)p>=0;positive terms cannot violate this. At most22 positive terms give real<seed+22p+22p<45p throughout. Imaginary has44 positive terms plusseed,<45p. All seven-limb add/sub chains are exact with no448-bit wrap. Final real is congruent to the desired native residue because22p=0modp.

The last real SBC intentionally preserves flags;the following imaginary ADDS starts a new carry chain. The shared interpreter implements that exact semantic. Range0..45p<49p preserves every premise of quotient_certificate.md:actual D=floor(p/2^325),computed seven-limbq*p slots0..49,estimate and one correction. Every output is canonical0..p-1. Both old third seed coordinates are loaded before pair writes;the other pair's coordinates remain untouched until read. Exact third/output alias is safe;other partial overlaps are forbidden.

## Private map frame and the unchanged public ABI

recurrence_fp2_map is a PRIVATE helper called ONLY by recurrence_fp2_step. It may clobber19..29 because the public wrapper has already saved them and uses only SP-based restoration after the call. It is NOT advertised as a standalone public PCS function. It uses30 for tablebase,so it still saves incomingLR at entry with STRx30,[SP,#-112]! and restores it with LDRx30,[SP],#112. No nestedBL exists. SP/LR/x18 are preserved. Remove the five redundant callee STP and five redundant LDP:exactly10 instructions,not uncharged caller work.

Frame112 is16-byte aligned:LR0..7,padding8..15,args16..47,seven-limb second-row spill48..103,padding104..111. Atomic offsetmap96->16,112->32,120->40,128->48,144->64,160->80,176->96 preserves every initialized read/write. Biasheader pointer isSP40,outputpointerSP16. No padding or stale208-frame slot is read. Both pre/post112 offsets and every pair/single offset are encodable. Wrapper528+map112=640;key nested maximum816 remains the global stack bound.

## Exact source counts and complete resources

Each COLUMN has6 scalar input LDR,11 extracts,5 crossing EXTR/5MOV,two index/base operations,11 MSTEP and10 advances:270 instructions. Width35 stride96*2^35=3TiB is oneMOVZ0x300<<32;34 advances useLSR1 in the SAME ADD. Column27TiB is oneMOVZ0x1b00<<32. Crossed digits load next limb,EXTR/mask into3 and retain nextlimb in19 before table loads clobber20. No digit is lost and no source read goes beyond six limbs.

Each row pair is9 seed +13 bias +4*270 columns +4 spill +2 header +31 QREDUCE +3 stores +4 reload +31 QREDUCE +3 stores =1180. Main entry/args/context/constants/exit/RET total13 after private-map saving. Thus map2*1180+13=2373. Wrapperown290+map2373=2663 throughRET. Caller16+pointerrotation4=2683,retained2688 with5 spare. The private-map revision saves10 actual and6 previouslyunused recurrence allowance instructions;the inner256 and final75 full-ledger reserves are retained.

Generation2688+824key+376hash+22aggregatescalar+704main=4614. Two local sort pairs720+scan128+cheap1+expensive1+spare256=5720;network933 giveshot6653,retained6728 with75spare. SAMEunequalfullL=88722074495999233395607723401466085376,N=6728L+2^124 satisfies EXACT N^1000<20000^1000*2^124489,proposed124489/1000=124.489. Fullp381/r255,C20000/success0.39/sixtracks/onecomputationmetric,cold124/peak138/advice0 and exact0.3900120659047603 success remain unchanged. All failed/trial/preprocessing/random/transport/verification/output paths remain paid.

Shared-only audit asserts ACTUAL2663publicstep/2373privatemap against an independent DENSE full-p381Fp4 oracle. ABI assertions for19..29/LR/SP/x18 are at PUBLICstepreturn,where the wrapper restores them. Private helper clobbers are explicit. Readonly initialized virtual table fixtures have current two27TiB/11digit geometry and fixed64TiB;they are scoped test oracles,never substitutes for charged production construction. Other full-p/key/369hash/fullY/U/corner/full-rscalar/byte/copy/packet tests remain. These are not productionGTinversion/tablefill/controller timing/probability/hardware qualification. No new pass/score/acceptance/promotion is invented;candidate code was NEVERhostexecuted/imported/compiled. Shared credential-free llm-verifier alone executes/reviews/scores.
