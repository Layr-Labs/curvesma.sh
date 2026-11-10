# Centered382-bit difference codes and exact full-field recurrence

For native canonical Montgomery input residues Rm=2^384,the exact field identity is cX-c^p²Y=A(a-d)+Bxi(b+e)+t*(A(b-e)+B(a+d)),where c=A+B*t,t²=xi=1+u,u²=-1. p² fixesFp2/negates t. Actual c=3Phi(g) is constructed cold from the actual fullGT increment g=h orh^-b;norm1 gives third coefficient1. All final raw-history coordinates remain canonical after unchanged quotient reduction and third seed. No structured c,input conditioning or reduced curve is assumed.

## Raw differences and exact signed-digit interpretation

Six-limb SUB/SBC stores the FULL384-bit two's-complement difference E=(a-d) mod2^384. Canonical inputs a,d in[0,p) imply|a-d|<p<2^381. Thus its low382 bits form a UNIQUE signed382-bit code:if bit381 iszero,value iscode;otherwisevalue iscode-2^382. Stored highbits382/383 are exactly redundant sign extension (bothzero orbothone),not independent input bits. No source curve/order parameter or meaningful difference information is dropped. Four sums remain ordinary rawx+rawy<2p<2^382.

Shared A columns alone consume differences. Their final35-bit top digit atstart347 is interpreted as t=d ifd<2^34,else t=d-2^35. Store canonical t*C_native*2^347/Rm. Other digit tables are canonical unsigned linear terms. Summing digits then equals C_native*(a-d)/Rm modulo p for BOTH positive/negative differences. This is an exact bounded signed-integer identity,not a probabilistic sign heuristic. ALL2^35 top entries are constructed;there is no hot branch,biasADD or mask. B/Bxi sum tables retain ordinary full382-bit interpretation.

## Complete paid construction and memory

A columns0/4 share(Ar,Ai),1/5 share(-Ai,Ar). Reuse only these exact identities within this ONE cold attack. Two shared A columns use seven35/three34/final35-bit digits atstarts0,35,70,105,140,175,210,245,279,313,347,total382;each28.5TiB. Four unique B/Bxi columns2/3/6/7 use eleven32/final30,total382;each4.21875TiB. Every96-byte entry contains two canonical fullFp values. All table widths include EVERY entry,not only observed codes.

Construct each real coefficient/seed with pinned_mont_core.S/wrapper and actual full-field arithmetic. All ordinary unsigned tables use fill_pair.S fromzero. For each signed A top table,construct actualincrement seed=C_native*2^347/Rm. Fill its firstFULL2^34 entries fromzero with fill_pair. Construct canonical offset=-2^34*seed by34 full-p doublings andzero-minus negation ofbothcoordinates,all paid. Fill its secondFULL2^34 entries atbase+2^34*96 fromthat offset using fill_pair_offset in the SAME sourcefile. At index2^34+j this is(-2^34+j)*seed=(d-2^35)*seed. Thus the two half-fills construct EVERY top entry without a runtime sign operation or free advice. Both explicit scalar loops use two full-p additions per later entry. Every seed,power,offset,half-prefix,fill,page,write andtransport is charged.256 scalar instructions perentry plusfixedwork fitscold2^92 subset;globalcold<2^124,algorithmdata<2^110,peak<2^138,advice0 remain.

Vectorbase4.03125TiB=0x408*2^32,size73.875TiB ends77.90625. Fixed-s base77.90625TiB=0x4de8*2^32,size2.0859375 ends79.9921875. Parent4 has32GiB gap beforevectors;children[80,128) have8GiB gap afterfixed. Half-open regions aredisjoint,1GiB-page aligned andbelow2^47. Both base constants cost oneMOVZ. Ordinary page/processor/edge initialization remains fully charged;no vector table is hidden in private4MiB bootstrap or a parent mailbox. No physical allocation on this VM is claimed.

## Exact embedded field body2368

Perrowpair thefirst Atransformloads actualtablebaseSP56 andextractsall382 differencebits;itcosts264. The secondA transformadvances its finaltable by3TiB/resetstride andleavesnormalBbase DIRECTLYinx26. Row1 usesphysicalnormalcolumns0/1,so neitherneedsx2 blockconstant:ci1 compile-timebranch omitsitsMOV. Row2 usesphysical2/3,so ci5 keepsit. Thisconstant compile-time distinctionchangesnoproductiondata. A pairsnowcost531/532,total1063.

Firstnormalcolumn0 usesalready-correct x26 andneedsNO baseMOV. Nextphysical1base isprevious finaltablebase+96*2^30 viaADDx26,x26,x27,LSR2. Physical2 startsfromdirectnormalbase plus2*block;physical3usesexactsameadjacent increment. Thusnormalcost281+3*282=1127. No x30 registerisreadorwrittenanywhereinactualleaf;RET's callerlinkremainsuntouched. Fourreductions108,seedloads/initializations16,secondaccumulator/const20,outputstores12,stride/elementsetups4 givefieldbody1063+1127+108+16+20+12+4=2350. Allsignedtopvalues,coefficientmatrices,canonicalterm/49p bounds andoutputbytesareunchanged.

## Complete fused leaf2465 and current ledger

One528-byteframe:sixusedcallees19/21at0,26/27at16,28/29at32;out/tableat48/56,constproductsbase64. Fiveinitialstorespreserveallcallees/arguments;callees20/22..25/LRareuntouched. Workspace80..135,padding136..143,temp384at144..527aredisjoint. ADDx29,SP144keepstheSAMEoriginalSP-384temporarybase. ActualfieldmacrosloadtableSP56/constSP64;outreloadSP48matchesitssave.

MSTEP/WSTEPinterleaveeachactualLDPx3,x21withthetwoADDS/ADCSusingthatpair,thenloadthenextpair. OrdinaryLDPdoesnotmodifyNZCV,soeverycarrypropagatesexactlyasin thepreviousload-all-firstsequence. TherearestillthreeLDP/sevenaddinstructionspervector,withallactualtablewordsloadedbeforeuse. Onlyx3/x21holdtablelimbs;callees22..25areneverused. Inputlimbx19/adjacentx2andaddressesstayliveunchanged. Everyoneofthetwofullaccumulatorsstillreceivesthesamecanonicalterm,withoutdroppedreads/carries,bulkcopy ornewtableoracle.

Exactcount:SUBSP+fivecallee/argumentstores+ADDtempbase7;fourBOTH96;sixinputpointerADD;oneoutLDR;fieldbody2350;threecalleeLDP+ADDSP+RET5. Total7+96+6+1+2350+5=2465. Sixcalleesinstead oftenremovetwosaves/tworestores,4actualinstructions. Thefieldbody/alldigit/quotientarithmeticcountsareunchanged. Frame528aligned16;allSPimmediates/offsetsencode. Standalonemapreferenceisuncalled/notadditionallyhotcharged.


Caller16+fourrotations gives2485,retained2506 with21spare. Nativekernel/key/NH pricesareunchanged forbothfull-r streams. Non-scalar generation4202L=(2506+984+712)L. Actual scalar allowance12m+28n charges atmostm baby andn giant records onEVERY success/failure path;no balanced20L average isused. Local720L+scan128L+confirmation1L+inner spare256L,network920L andouter spare80L give N=6307L+12m+28n+2^124. EXACT N^5000<20000^5000*2^621767 gives124.3534. Bothfull bucket constructions,cold2^124/peak2^138/advice0,referenceC20000/threshold0.39/fullfields/sixtracks/solemetric remainunchanged.

## Scoped independent checks

ALL16 shared cases interpret thecomplete2465-instruction fusedleaf,assert no standalone mapcall,andcompare canonical finaloutputs to independent densefull-p Fp4. They also inspect all eight actual temporary integers. Raw differences equal exact(a-d) mod2^384;low382 signed decode equalsa-d andtwo highbits equaltheir requiredsignextension. Case3 nativep-1 forcesall sums2p-2>=2^381. Case4 rawx=p-1/rawy0 forcesALLpositivedifferencesp-1 in firstsigned top halves. Case5 rawx0/rawyp-1 forcesALLnegative differences-(p-1),stored2^384-p+1,insecondsigned halves. Equal inputs producezero differences directly. Bothsigns,zero andextended sums are genuinely exercised.

For each16 full-p seedcases,the actual nativezero-fill andoffset-fill symbols are interpreted on2..5 entries (32 helperexecutions total),includingzero/maxseeds. Everycanonical entry equalsindependent initialoffset+digit*seed,withallbytes/initializedreads/ABI/coldcount checks. These bounded fixtures areNOT production2^34 half-fills orfree tables;full loops areconstructed/paid and justifiedbycanonical-add induction. Virtual accessedcells are scoped aids. All16 native/hash/key/scalar/protocol cases,80 copies,all-path parent andexactledger checks remain. Finite evidence doesnot provegroup-event theorem,GTsampling,physicalhardware,productioninversion ortotalcomplexity.


Current quotient refinement:actualQREDUCE27 selectsanactualtableaddressone64-byteslotearlierONLYontrialborrow,preservingoriginalaccumulatorandreloadingtheexactordinaryq*p. Fourreducerssave8instructionsrelativetoQ30;fieldbody2350,wholeleaf2465,caller20gives2485,retained2506,non-scalar6307L. Standalonemap31/fixed28 remainunchanged. Allconstruction/layout/tablevaluesareunchanged;theexact28-instructionproofisinquotient_certificate.md.


CURRENT product-base argument contract:afterconstructingtheSAME3728-byteblock,coldcomputesproductsbase=blockbase+528 onceandpassesitin theexistingconstargumenttox4forfused/key,x3forfixed. Parent/oddpassitunchanged;theirseparatep48argisstillsuppliedforFpadd/sub/negation. No additionalhotargument/ADD ishidden. ThisonecoldADDperworker ispaidintheunchangedcoldallowance. Allheader/padding/fiftyslotsremainfullyconstructed/delivered.

Atproductbase-16,anencodableLDPx2,xzr readsD@block512andinitializedzeropad@520;oneordinaryinstruction replacesLDRD. Negativepairdisplacement-16iswithintheordinaryAArch64signedscaledimmediate. Thediscardedxzrwordisactuallyreadandinitialized;noout-of-boundsreadoccurs. BothhelpersnowformqhatproductaddressinONEADDfromproductsbase,eliminatingonebaseADDperreduction. Trialborrowselectstheprecedingslotexactlyasproved. Ifqhat0,theunselectedtemporaryaddressisblock464andneverread. FusedQ27/fixedQ28savefourfusedinstructionsandoneperfixedcall. Nativeleaf2465 retains2506,fixed192,odd504,key655 retains666;non-scalar4202generation/6307total plus12m+28n+cold2^124. Main712 isunchanged. Fulltablegeometry/field/probability/inputdistribution/verificationareunchanged. Standalonemapreferencekeepsitsoldblockbasecontractandisnevercalledbytheactualfusedleaf.
