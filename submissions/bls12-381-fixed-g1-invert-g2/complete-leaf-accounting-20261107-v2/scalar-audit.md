# Complete scalar call-graph bound

All counts are scalar AArch64 instructions,not cycles. No SIMD,crypto extension,free oracle,arbitrary-width arithmetic primitive or uncharged transfer is used. All code is analytic source evidence and has not executed on the host or with credentials. One blst dense Fp12 multiplication per transition is charged as the protected reference operation; its internals are excluded from R to avoid double counting. Everything else ordinary is in the five submitted assembly files.

## Complete path certificate

hot-loop.S ordinary main path:62 instructions outside callees,including all BL/arguments,pointer swap,index/counters,block clearing,CV copies,label extraction/store,prefix cache,lookup flags and final branch. A last top termination test may run once more; it is fixed exit work inside2^90. Checkpoint/found/exhausted handlers are cold bounded events,not omitted ordinary work.

fapi_scalar_update:46 instructions (5 pointer/constant setup,20 per residue,RET). Current/jump operands are canonical four-limb residues. r<2^255 prevents a fifth sum limb. SUBS/SBCS compute a separate candidate reduction;CSEL uses final carry to choose candidate/original,then stores to alternate record. Both residues are canonical. All source/destination bases are distinct,offsets576..632 aligned; only caller-saved x0..x16 are used. Its interface is already in the62 main instructions.

fapi_copy_key:73 instructions (36 LDR/STR pairs+RET). Source B_i offset48*floor(i/6)+8*(5-(i mod6)),output8*i. It preserves all callee state,uses x0,x1,x3 only,and copies36 numerical b words. Sign is separately stored at928,prefix at936;copy never overwrites them. All field coordinates are canonical native Montgomery residues from the reference arithmetic,with cold input/table conversion paid.

fapi_tabulate:892 instructions exactly (5 accumulator MOV,73*12 body,5 REV/STR pairs,RET). Each body has digit load,table-base load,two ADD for x40 address,two LDP plus one LDR for the five-word entry,and five EOR. Sign digit is0/1;other72 digits are native32-bit halves at offsets0..284. Tables and40-byte records are8-byte aligned,and all addresses fit64 bits below2^50. Inputs x0/key,x1/table-base-array,x2/output;accumulators x3..x7,temporary x8..x12. No internal call,stack or callee-saved register. Every input read,output store,address and serialization is counted. The serialized40 bytes go directly into the cleared final block.

fapi_sign:89 instructions bound by the entire forward-only body,including both mutually exclusive terminal paths. Six coefficient tests use6 LDR,5 ORR andCBNZ each,72 total. On all zero,MOV/STR/RET adds3. Selected coefficient usesADRP/ADD2,three LDP of H=(p+1)/2,six SUBS/SBCS,CSET,STR,RET:14. The full89-instruction listing bounds every acyclic path (actual worst selected path86). No backward edge,stack,call or callee-saved register exists. A valid sign depends on the first nonzero coefficient:canonical x compared with p-x gives sign1 iff x>=H. Therefore the six-limb comparison computes the specified lexicographic sign exactly,including all-zero case.

fapi_blake3_final:797 instructions exactly,including its96-byte frame andall preloads/returns. State W0..W15;ten message words W16,W17,W19..W26;X27 output/CV,X28 block,X29 IV. x18 untouched,x19..x30 saved/restored. Mixing56*12-2+70=740;fixed57 consists6 saves,4 pointer/IV addressing,10 message loads,14 state initialization,16 output and7 epilogue. All7 standard rounds/permutation/rotations/root flags and32 output bytes remain.42 zero message injections are identity omissions. First d12,d13 counter zeros feed directly into ROR,so their first EOR and initial zero MOV can be folded exactly. No message word is modified after its preload;G destinations are only state registers. No runtime schedule,hidden spill,internal call or additional input normalization occurs. The absolute input/output pointers remain valid through the leaf's frame adjustment.

fapi_prefix_lookup ordinary part:<=128 per invocation,plus globally bounded full-confirmation work below. Entry4;at most10 lower_bound iterations11 each;last loop test2;terminal bounds/prefix miss6;return3=125. Extra3 for final zero-budget check fits128 (also conservatively covered by confirmation slack). Equal-prefix scans confirm every37-word group key,never accept a prefix alone. The complete caller-saved leaf uses x0..x16 only. Prefix order and insertion preserve contiguous ties;no true match can be skipped.

For each full confirmation,scan6,budget6,cursors2,exact equality112 andadvance2/return3 cost<=129. Budget checks LDP/ORR/CBZ then128-bit SUBS/SBC/STP are explicit;cap initialization B_cap=floor(K/256) is cold. All attempted confirmations decrement once;before exceeding cap,the leaf returns flag2 and main aborts. Including a final exhausted check,confirmation work<=129B_cap+3<=256B_cap<=K (B_cap>=1 at protected parameters). Thus charge at most one extra scalar instruction per transition across the entire attack,not a false claim of worst-case no prefix collisions. Probability.md proves the added cap-failure bound2^-43 separately.

Complete nonreference bound per K steps:

 62+46+73+892+89+797+128 =2087 ordinary instructions per step;
 all capped full confirmations <=K additional instructions;
 total <=2088K <=2112K.

The declared2112 leaves24 instructions/step beyond these conservative maxima. Every ordinary call interface is in the62 main instructions,so none is assumed free or duplicated in a missing wrapper. All helpers preserve persistent x19..x28,with BLAKE3 saving/restoring the ones it reuses. No ordinary allocation or uncontrolled spill occurs. Reference BL/interface is in main62 while its multiplication internals are counted once as H=1. Function/control placement can use direct calls;code/constants fit a small linked image with no required PLT/veneer per call. Cold linking/initialization is fully in2^90.

## Exact key and probability contracts

For X=a+bw in GT,p^6=-1 mod r gives a²-vb²=1. Given b there are at most two a values,a and-a;the lexicographic sign distinguishes them. Thus36 b words plus sign are injective without square roots. Physical record has576-byte X,64-byte A/B,296-byte exact key and8-byte cached prefix=944. Prefix is the first8 U bytes loaded natively;it is only an ordering aid. Key equality excludes prefix. All bytes are initialized by direct arithmetic/key/sign/prefix writes or cold setup;there is no missing per-transition copy.

Simple tabulation on72 32-bit halves and one sign digit has pairwise full collision2^-320 andprefix collision2^-64 under H. A cached-group-label surrogate independent of tables permits the full-K union/expectation bounds despite adaptive actual states. At most1023K unequal query/checkpoint pairs give E confirmations<=1023K/2^64+1. With B_cap>=K/512,Markov bounds exhaustion<2^-43. The hash bridge H remains the only cryptographic idealization;pairwise universality does not assert full prefix independence.

## Cold paths,storage and ABI

Cold setup validates/converts inputs,computes h,constructs every jump/table entry andinitial key/prefix/label,sets t=0,next marker=s,n=1,K andB_cap. Full construction,failures,randomness,arithmetic andmemory transfers are paid. The tables occupy716056947589200 bytes together,leaving over400 trillion bytes below2^50 for code,1023 records,stack andscratch. No production arrays are allocated on this VM.

At most1022 marker hits plus initial record copy need1023*118 copy instructions. All sorted-pointer insertions,marker/count updates andhandlers stay below2^30 instructions,inside2^90. The top t/K bound andnext-marker checks are explicit;1023s>K andK<2^128 at protected constants. Found handler recovers/validates scalar andcanonical G2 output for supplied P;exhaustion aborts. Two pairings,field/group checks,output,all cold arithmetic andconstant construction retain the complete2^90 bound in proof.md.

An aligned288-byte controller frame saves original LR at48 andx19..x28 at192..271,budget272..287. Controller enters loop by branch;terminal handlers restore original LR/frame;checkpoint handlers resume by branch. BLAKE3's separate96-byte frame explicitly restores x19..x30. Other leaves never alter callee state. No return through the last helper's LR is assumed.

This certificate is static analytic instruction accounting,not a measured trace,security theorem orproduction execution. Shared llm-verifier owns scientific review and scoring. The independently checked total inequality in probability.md establishes only the proposed bound conditional on the complete stated algorithm and H.
