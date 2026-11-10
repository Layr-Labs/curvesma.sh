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

recurrence_fp2_map is a PRIVATE helper called ONLY by recurrence_fp2_step. It may clobber19..29 because the public wrapper has already saved them and uses only SP-based restoration after the call. It is NOT advertised as a standalone public PCS function. It uses30 for tablebase,so it still saves incomingLR at entry with STRx30,[SP,#-48]! and restores it with LDRx30,[SP],#48. No nestedBL exists. SP/LR/x18 are preserved. Remove the five redundant callee STP and five redundant LDP:exactly10 instructions,not uncharged caller work.

Frame48 is16-bytealigned:LR0..7,padding8..15,args16..47. OutputpointerSP16/headerSP40 are unchanged. Both accumulators staylive with scratch19..25;no spill slot or stale112/208-frame read remains. Pre/post48 and all other native offsets are encodable. Wrapper528+map48=576;keyglobal816 isunchanged. QREDUCE loadordering/canonicality/live-register proof below is authoritative.

## Exact source counts and complete resources

Each COLUMN has6 scalar input LDR,11 extracts,5 crossing EXTR/5MOV,two index/base operations,11 MSTEP and10 advances:270 instructions. Width35 stride96*2^35=3TiB is oneMOVZ0x300<<32;34 advances useLSR1 in the SAME ADD. Column27TiB is oneMOVZ0x1b00<<32. Crossed digits load next limb,EXTR/mask into3 and retain nextlimb in19 before table loads clobber20. No digit is lost and no source read goes beyond six limbs.

Each rowpair is9 seed+13bias+1080columns+2header+31realQ+3realstores+31imagQ+3imagstores=1172. Main entry/context/constants/exit/RET13 givesmap2357. Wrapperown290+2357=2647;caller16+rotation4=2667retained2672with5spare. Both accumulatorgroupsremainlive in disjoint registers;no stack spill/reload or copybridge is hidden.

Generation2672+824key+376hash+22aggregatescalar+704main=4598. Two local sort pairs720+scan128+cheap1+expensive1+spare256=5704;network933 giveshot6637,retained6712 with75spare. SAMEunequalfullL=88722074495999233395607723401466085376,N=6712L+2^124 satisfies EXACT N^500<20000^500*2^62243,proposed62243/500=124.486. Fullp381/r255,C20000/success0.39/sixtracks/onecomputationmetric,cold124/peak138/advice0 and exact0.3900120659047603 success remain unchanged. All failed/trial/preprocessing/random/transport/verification/output paths remain paid.

Shared-only audit asserts ACTUAL2647publicstep/2357privatemap against an independent DENSE full-p381Fp4 oracle. ABI assertions for19..29/LR/SP/x18 are at PUBLICstepreturn,where the wrapper restores them. Private helper clobbers are explicit. Readonly initialized virtual table fixtures have current two27TiB/11digit geometry and fixed64TiB;they are scoped test oracles,never substitutes for charged production construction. Other full-p/key/369hash/fullY/U/corner/full-rscalar/byte/copy/packet tests remain. These are not productionGTinversion/tablefill/controller timing/probability/hardware qualification. No new pass/score/acceptance/promotion is invented;candidate code was NEVERhostexecuted/imported/compiled. Shared credential-free llm-verifier alone executes/reviews/scores.

# Keep both seven-limb accumulators live through exact quotient reduction

Only scratch allocation/load ordering in QREDUCE changes. The arithmetic function is identical. For the real seven-limb accumulator4..10 and imaginary11..17,define QREDUCE(a0..a6) parameterized by its accumulator registers. Quotient estimate/slot addressing use1/2/26/28 asbefore. Load seven ordinary q*p limbs into19..25,subtract from a0..a6 with SUBS/SBCS. ONLYTHEN reload six fullp limbs into19..24 and CSEL them to p/0 using ORIGINAL finalborrow. LDP does NOT alter flags,so those three delayed loads preserve the condition exactly. ADDS/ADCS/ADC restorespifborrowed,identical to the prior proven correction. Allq*p registers are dead before p replaces them.

The macro does NOT write11..17 when reducing4..10. Its scratch set1/19..26 is disjoint from both accumulator groups;2/28 are readonly context and0/output29/input30/tablebase are preserved. Imaginary accumulator can therefore stay live while real is reduced/stored. Invoke the SAME31-instruction macro directly on11..17 next and store11..16. No spill,no reload,no seven MOV bridge,not an uncharged copy. Canonical high limb/estimated-q/one-correction proof remains exact because0<=Y<45p<49p and all50ordinaryq*p slots are unchanged.

Delete4STP/STRspills and4LDP/LDRreloads PERROWPAIR,eightinstructions each,sixteentotal. Newrowpair9seed+13bias+1080columns+2header+31realQ+3realstores+31imagQ+3imagstores=1172. Entry/context/constants/exit/RET13,so map2*1172+13=2357. Fullwrapperown290+map2357=2647throughRET. Caller16+pointerrotation4=2667retained2672with5spare. Thisrevision saves16actualand16retainedinstructions;allpriorABI/source/table/key/cap/parametercontracts remain.

Private map'sframe48 containsLR0..7,padding8..15,args16..47. No spillregion exists. Pre STRLR[SP,#-48]! andpostLDRLR[SP],#48 are encodable and preserve16-bytealignment. Outputpointer remainsSP16,headerSP40. Paddingisneverread. RemovealloldSP48/64/80/96spillreferences. Mapstillprivateclobbers19..29;thepublicwrapper's528-byteframefullyrestoresALL19..29/LR/SP/x18. Nestedrecurrencestack528+48=576;globalkey816unchanged. No external standalonePCSmapcall is advertised.

Generation2672recurrence+824key+376NH+22aggregatescalar+704main=4598. Two localsortpairs720+scan128+cheap1+expensive1+spare256=5704;network933giveshot6637,retained6712with75spare. CurrentL=88722074495999233395607723401466085376,N=6712L+2^124 satisfiesEXACT N^500<20000^500*2^62243,proposed62243/500=124.486. Fullp381/r255,C20000/success0.39/sixtracks/onecomputationmetric,cold124/peak138/advice0andexactsuccess0.3900120659047603remainunchanged. Everypreprocessing/trial/failedprefix/transfer/random/input/output/verification operation remains paid.

Shared-onlyauditassertsACTUAL2647publicstep/2357privatemapandfullpubliccallee/LR/SP/x18restorationagainstindependentdensefullp381Fp4oracle,currentvirtualreadonlytwo27TiBvectorgeometryandfixed64TiBbase. It checks initializedreads,writes,carry/borrow,aliases andencodableimmediates. Thisisnotproductiontablefill/GTinversion/probability/timing/hardwarequalification. Static31-opmacrocount/sourceconstructionwasusedonhost;NOcandidateinstructionwasexecuted/imported/compiled/emulatedthere. Allactualcandidateexecution/review/scoringbelongtosharedcredentialfreellm-verifier. No new pass,score,acceptance or promotion is invented.
