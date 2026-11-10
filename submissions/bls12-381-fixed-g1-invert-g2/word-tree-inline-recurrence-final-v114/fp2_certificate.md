# Complete current recurrence: shared vectors, one modulus load, one public leaf

All residues are canonical full protected Fp381 Montgomery values Rm=2^384. Write c=A+B*t,X=a+b*t,Y=d+e*t in Fp4/Fp2,t²=xi=1+u,u²=-1. The exact recurrence factor is
A(a-d)+Bxi(b+e)+t[A(b-e)+B(a+d)].
Set V0=a-d,V1=xi*(b+e),V2=b-e,V3=a+d. Compute A*V0+B*V1+t*(A*V2+B*V3),plus the old third history. Only A and B coefficient tables are needed;the complete GT recurrence and all full field coordinates are unchanged.

## Actual fused wrapper and the public register contract

recurrence_fp2_step is the ONLY public PCS entry. Its528-byte frame saves19..29/LR before ANY clobber,and restores ALL after the map. SP/x18 are preserved. Four BOTH transforms construct a±d,b±e in initialized384 bytes atSP144..527. The fifth BOTH uses inputSP240/SP288 and writes differenceSP240/sumSP288,exactlyxi*(b+e). Every input is read before any output store;SUM25..30 remains live after MINUS stores. This exact crossed output alias is safe and touches neither raw histories nor other transformed components. ALL wrapper frame accesses useSP because BOTH clobbers29/LR.

Load the full six p limbs ONCE into19..24 before the first BOTH. No BOTH writes19..24 and no intervening setup touches them. Every full-width CSEL/reduction therefore sees the same correct p. This does not require the external caller to preload anything. Each BOTH body is48 instructions after removing three repeated p loads. Fivebodiesremain48instructionseach;theinlinedpublicleafstructure/countareprovedbelow. Thepriorseparatewrapperown290isnotacurrentstandalonecost. Hoisting alone removed16 actual instructions from the preceding306-own wrapper. Full fields/output/carry/alias behavior are identical.

## All ordinary coefficient and fixed-s tables are constructed

For C=A,B,each vector entry is(d*Cr*2^start/Rm modp,d*Ci*2^start/Rm modp),96 bytes. Real input digits add(first,second) to(real,imag);imaginary digits SUBTRACTsecond fromreal and ADDfirst toimag. This is exactly C*(xr+xi*u) foru²=-1. Both output row pairs reuse these same A/B tables. No negated duplicate table or oracle is supplied.

ALL381 significant input bits have widths[35,35,35,35,35,35,35,34,34,34,34],starts0,35,70,105,140,175,210,245,279,313,347. One vector column is96*(7*2^35+4*2^34)=27TiB;two are54TiB. EVERY entry is constructed by actual full-p pinned seed products and existing fill_pair75(n-1)+17<=96n canonical additions,including all zero entries and headers. Full64-bit counters handle2^35/34 sizes. This is ordinary constructed memory,not virtual production arithmetic for free.

Vectors occupy8..62TiB;unchanged fixed-s13.5TiB is moved by its existing base argument to64..77.5TiB. Parent4TiB/private<1TiB/childports80..128TiB remain. All ranges are disjoint and pointers strictly below2^47. Other routing maps and DATA identifiers are unchanged. All new table pages/seeds/constants/MMU work is charged;bothworkers' vector/fixed fills remain<2^50 instructions inside completecold124. Bothworkers'tables<2^48 bytes/globalpeak138;no physical VM memory qualification is claimed.

## Signed accumulators stay nonnegative and below45p

Each real accumulator starts OLDthirdreal+22p. The ordinary seven-limb22p is already constructed in the50-slot header at528+64*22=1936. Exactly22 canonical negative table terms are subtracted. After k<=22 such terms,the accumulator is at least seed+(22-k)p>=0;positive terms cannot violate this. At most22 positive terms give real<seed+22p+22p<45p throughout. Imaginary has44 positive terms plusseed,<45p. All seven-limb add/sub chains are exact with no448-bit wrap. Final real is congruent to the desired native residue because22p=0modp.

The last real SBC intentionally preserves flags;the following imaginary ADDS starts a new carry chain. The shared interpreter implements that exact semantic. Range0..45p<49p preserves every premise of quotient_certificate.md:actual D=floor(p/2^325),computed seven-limbq*p slots0..49,estimate and one correction. Every output is canonical0..p-1. Both old third seed coordinates are loaded before pair writes;the other pair's coordinates remain untouched until read. Exact third/output alias is safe;other partial overlaps are forbidden.

## One full public frame and two live accumulators

ThereisNOstandalonemaphelper.Themacrofileisexpandedintorecurrence_fp2_stepusingitsoriginal528-byteframe/outputSP96/headerSP128/tempSP144. Originalcaller19..29/LRvaluesaresavedonceandfullyrestored. Both7-limbaccumulators4..10/11..17staylive;QREDUCE31usesqp19..25,subtracts,THENreloadsp19..24withflag-preservingLDP. Theotheraccumulatorisuntouched;reducetheseconddirectlyin11..17. No spill/reload/copybridge or privateLRframeisneeded. Wrapperstack528;globalkey816unchanged. Fullregistered/canonicality/sourceproofbelow.

## Exact source counts and complete resources

Each COLUMN has6 scalar input LDR,11 extracts,5 crossing EXTR/5MOV,two index/base operations,11 MSTEP and10 advances:270 instructions. Width35 stride96*2^35=3TiB is oneMOVZ0x300<<32;34 advances useLSR1 in the SAME ADD. Column27TiB is oneMOVZ0x1b00<<32. Crossed digits load next limb,EXTR/mask into3 and retain nextlimb in19 before table loads clobber20. No digit is lost and no source read goes beyond six limbs.

Eachrowpairis9seed+13bias+1080columns+2header+31realQ+3realstores+31imagQ+3imagstores=1172. Inlinecontext/constants9givescore2353;remainingwrapper286givesonepublicleaf2639. Caller20=2659retained2664with5spare. Everyinput/table/fieldcorrectionandoutputread/writeisretained.

Generation2664+824key+376hash+22aggregatescalar+704main=4590. Two local sort pairs720+scan128+cheap1+expensive1+spare256=5696;network933 giveshot6629,retained6704 with75spare. SAMEunequalfullL=88722074495999233395607723401466085376,N=6704L+2^124 satisfies EXACT N^250<20000^250*2^31121,proposed31121/250=124.484. Fullp381/r255,C20000/success0.39/sixtracks/onecomputationmetric,cold124/peak138/advice0 and exact0.3900120659047603 success remain unchanged. All failed/trial/preprocessing/random/transport/verification/output paths remain paid.

Shared-only audit asserts ACTUAL2639publicleaf/nohelpercall against an independent DENSE full-p381Fp4 oracle. ABI assertions for19..29/LR/SP/x18 are at PUBLICstepreturn,where the wrapper restores them. The macrocore is within the public frame;no standalone helper ABI is asserted. Readonly initialized virtual table fixtures have current two27TiB/11digit geometry and fixed64TiB;they are scoped test oracles,never substitutes for charged production construction. Other full-p/key/369hash/fullY/U/corner/full-rscalar/byte/copy/packet tests remain. These are not productionGTinversion/tablefill/controller timing/probability/hardware qualification. No new pass/score/acceptance/promotion is invented;candidate code was NEVERhostexecuted/imported/compiled. Shared credential-free llm-verifier alone executes/reviews/scores.


# Inline the recurrence core in its only public entry

recurrence_fp2_step is now one LEAF function with noBL. Its existing528-byteframe saves/restoresALL19..29/LR/SP/x18 asbefore. ThefifthBOTH mayclobber29/LR,buttheoriginalcaller'svaluesarealreadyatSP0/8andrestoredatend. No extra private helper frame or return address exists.

recurrence_fp2_map.S is now MACRO DEFINITIONS ONLY,with no global function entry or standalone PCS promise. recurrence_fp2_step.S explicitly includes that file. The shared-only parser prepends exactlythisknownadjacentinclude and checks its literal name;it does not invent arithmeticor omit runtimeinstructions. NativeassemblyandtheinterpreterexpandtheSAME MSTEP/COLUMN/QREDUCE/ROWPAIR bodies. Codeimage construction/copyingremainsinthechargedcoldbootstrap,notaper-recordfreecompiler.

Thewrapper'soriginalargsalreadyholdoutputpointerSP96,tablebaseSP120,constant3728pointerSP128. ROWPAIRthereforeloadstheexistingSP96/SP128slotsinstead ofprivateSP16/SP40. After thefiveBOTH transforms,LDRx0,[SP96],ADDx29,SP144,LDRx30,[SP120] establishoutput/temp/tablecontextDIRECTLY. Twooriginalthree-instructionstride/entrysize/columnsizeinitializationsandtwoROWPAIRsruninplace. EverytemporaryisstillSP144..527fullyinitialized;noarg,rawhistoryoroutputaliaschanges. ThelatepublicepiloguerestoresoriginalLRfromSP8despiteSUM/tableusingregister30.

Previouslywrapper4arg/BLinstructions+privatecore13context/entry/exitinstructions=17outsidebothROWPAIRs. Nowdirectcontext3+two3constantinitializations=9. Eightactualinstructionsaresaved:BL/privateLRentry/privateLRexit/privateRET/duplicatedargstores/contextcopies. No arithmetic,fieldload/tableentry/quotientcorrection/flagbarrierisremoved. IncomingflagsareirrelevantbecausefirstbiasADDSandallfullfieldchainsstartwiththeir ownflagsettinginstruction.

EachROWPAIRremains1172instructions. Inlinecore9+2*1172=2353. Thepublicwrapperportionafterremovingitsoldfourmap-callsetupinstructionsis290-4=286. Completeleaf286+2353=2639throughRET. Caller16+pointerrotation4=2659,retained2664with5spare. Thisrevisionreducesactualandretainedcostby8exactly;theinner256/final75reservesremain. Stackonly528,globalkeymaximum816unchanged.

Generation2664recurrence+824key+376hash+22aggregatescalar+704main=4590. Sort720+scan128+cheap1+expensive1+spare256=5696. Unchangedexplicitlyjustifiednotification/streamnetwork933giveshot6629,retained6704with75spare. SAMEfullL88722074495999233395607723401466085376,N=6704L+2^124 satisfiesEXACTN^250<20000^250*2^31121,proposed31121/250=124.484. Fullp381/r255,success0.39/C20000/sixtracks/onecomputationmetric,cold124/peak138/advice0andexact0.3900120659047603successremainunchanged.

Shared-onlyauditassertsACTUALpublic2639instructions,NONESTEDmapcall,allpubliccallee/LR/SP/x18preservation,andindependentDENSEfullp381Fp4oracle/canonicality/aliases/initializedreadwrite/immediates. Itreadsactualmacroincludeassetswithcurrenttwo27TiB/11digitgeometry andfixeds64TiB. Other369hash/fullY/U/scalar/key/copy/packetandexplicitnotificationfixturesretainactualscope. Noactualcandidatecodewasexecuted/imported/compiled/emulatedonhost;agent-ownedSTATICmacrocount2639isnotCPUexecution. No newpass/score/qualification/acceptance/promotionorproductionattackisinventedbeforegenuinesharedreview.
