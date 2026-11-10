# Centered382-bit difference codes and exact full-field recurrence

For native canonical Montgomery input residues Rm=2^384,the exact field identity is cX-c^p²Y=A(a-d)+Bxi(b+e)+t*(A(b-e)+B(a+d)),where c=A+B*t,t²=xi=1+u,u²=-1. p² fixesFp2/negates t. Actual c=3Phi(g) is constructed cold from the actual fullGT increment g=h orh^-b;norm1 gives third coefficient1. All final raw-history coordinates remain canonical after unchanged quotient reduction and third seed. No structured c,input conditioning or reduced curve is assumed.

## Raw differences and exact signed-digit interpretation

Six-limb SUB/SBC stores the FULL384-bit two's-complement difference E=(a-d) mod2^384. Canonical inputs a,d in[0,p) imply|a-d|<p<2^381. Thus its low382 bits form a UNIQUE signed382-bit code:if bit381 iszero,value iscode;otherwisevalue iscode-2^382. Stored highbits382/383 are exactly redundant sign extension (bothzero orbothone),not independent input bits. No source curve/order parameter or meaningful difference information is dropped. Four sums remain ordinary rawx+rawy<2p<2^382.

Shared A columns alone consume differences. Their final35-bit top digit atstart347 is interpreted as t=d ifd<2^34,else t=d-2^35. Store canonical t*C_native*2^347/Rm. Other digit tables are canonical unsigned linear terms. Summing digits then equals C_native*(a-d)/Rm modulo p for BOTH positive/negative differences. This is an exact bounded signed-integer identity,not a probabilistic sign heuristic. ALL2^35 top entries are constructed;there is no hot branch,biasADD or mask. B/Bxi sum tables retain ordinary full382-bit interpretation.

## Complete paid construction and memory

A columns0/4 share(Ar,Ai),1/5 share(-Ai,Ar). Reuse only these exact identities within this ONE cold attack. Two shared A columns use seven35/three34/final35-bit digits atstarts0,35,70,105,140,175,210,245,279,313,347,total382;each28.5TiB. Four unique B/Bxi columns2/3/6/7 use eleven32/final30,total382;each4.21875TiB. Every96-byte entry contains two canonical fullFp values. All table widths include EVERY entry,not only observed codes.

Construct each real coefficient/seed with pinned_mont_core.S/wrapper and actual full-field arithmetic. All ordinary unsigned tables use fill_pair.S fromzero. For each signed A top table,construct actualincrement seed=C_native*2^347/Rm. Fill its firstFULL2^34 entries fromzero with fill_pair. Construct canonical offset=-2^34*seed by34 full-p doublings andzero-minus negation ofbothcoordinates,all paid. Fill its secondFULL2^34 entries atbase+2^34*96 fromthat offset using fill_pair_offset in the SAME sourcefile. At index2^34+j this is(-2^34+j)*seed=(d-2^35)*seed. Thus the two half-fills construct EVERY top entry without a runtime sign operation or free advice. Both explicit scalar loops use two full-p additions per later entry. Every seed,power,offset,half-prefix,fill,page,write andtransport is charged.256 scalar instructions perentry plusfixedwork fitscold2^92 subset;globalcold<2^124,algorithmdata<2^110,peak<2^138,advice0 remain.

Vectorbase4.03125TiB=0x408*2^32,size73.875TiB ends77.90625. Fixed-s base77.90625TiB=0x4de8*2^32,size13.5 ends91.40625. Parent4 has32GiB gap beforevectors;children[92,128) have608GiB gap afterfixed. Half-open regions aredisjoint,1GiB-page aligned andbelow2^47. Both base constants cost oneMOVZ. Ordinary page/processor/edge initialization remains fully charged;no vector table is hidden in private4MiB bootstrap or a parent mailbox. No physical allocation on this VM is claimed.

## Complete current fused leaf and exact all-input count

One528-byte aligned frame:19/21 saved0,26/27saved16,28/29saved32,out/table48/56,workspace80..135,padding136..143,temp384at144..527. ADD29,SP144 givesentrySP-384tempbase. SP64stores/restorescallee22,notaproduct-pointerstore. Products24/D25/entry96in28arePAIDcold/read-only;worker27=48,23LUTrootunchanged. Fieldtemporarilyuses26fortableaddresses/29tempbase andACTUALLYrestoresboth(andallothersavedcallees);noBL/no18/LRwrite. Rawhistory/transformedinputs/fullcanonicaloutputs unchanged.

MSTEP/WSTEPinterleaveLDP3/21withtheactualtwoADDS/ADCSusingthatpair. LDPdoesnotmodifyNZCV,soeachcarryreachesnextlimbunchanged;ALLsixactualtablelimbsareloaded. Entry28=96 drives ordinaryMADD/shiftedADD. Advancesareseven35/three34/finalnonadvance35,normal32/final30. Columnbase1/3=28LSL30,base2actual2LSL1;wideboundaries28LSL35. SAMEcomplete96-bytevectorgeometry/nohotfreeconstants. FirstLDP0,26loadsout/table;staticfirstcolumnskipsredundanttableload,latercolumnretainsitafterreducerclobber. BOTHsumhighsaved21avoidsclobberingcoldproducts24. Allreadsinitialized/aligned/encodable.

Current exact bound N=6005L+12m+28n+2^124;N^2500<20000^2500*2^310649 gives124.2596. Fieldnormal2345/MAX2421/ret2386,keynormal806/MAX842/ret827;allmargins/globalrare1Lunchanged. Detailedledgerproof.md.

## Scoped independent full-p checks

For each16 full-p seedcases,the actual nativezero-fill andoffset-fill symbols are interpreted on2..5 entries (32 helperexecutions total),includingzero/maxseeds. Everycanonical entry equalsindependent initialoffset+digit*seed,withallbytes/initializedreads/ABI/coldcount checks. These bounded fixtures areNOT production2^34 half-fills orfree tables;full loops areconstructed/paid and justifiedbycanonical-add induction. Virtual accessedcells are scoped aids. All16 native/hash/key/scalar/protocol cases,80 copies,all-path parent andexactledger checks remain. Finite evidence doesnot provegroup-event theorem,GTsampling,physicalhardware,productioninversion ortotalcomplexity.

## Retained carry and quotient proofs

Current carry liveness proof:9p<2^384<10p. Eachseed/EVERYtabletermcanonical[0,p),includingsignedtopoffsetfills. Seed+firstk<=8terms<(k+1)p<=9p<2^384,soeachsix-limbcarryOUTzero. Physical10/17DEADuntilNINTHfirst-columntermADC10/17,xzr,xzr initializescurrentcarry0/1. ALLlaterhighADCsreadinitializedhighnormally;secondwide/normalcolumnsretainallupdates. Staticinit=1-secondONLYninthWSTEP,nohotbranch. First8highADCsperaccumulatoromitted,twoaccumulators/tworows=32;firstMOV10/bothMOV17unneeded. Fixedscalarfirsttermcanonical;terms2..9sum<9p,first8highADCs omitted/11DEAD untilterm10ADC11,xzr,xzr,last2readinitialized11. InitialMOV11also removed. Everylowlimb/digit/tableloadretained,nextlowchainstartsADDS/resetsflags. Native48omitted-siteassertionscheckACTUALcarry0,notirrelevantdeadregister0. Completefullporacles/ABI/initializedreadsretained;deterministicWORSTCASEproof,notexpectedwork.

Current quotient-reduction dead high:FULLseven-limbTRIAL productload/subtraction andborrowselectqhat/qhat-1MUSTremain. Approxquotientbounds/50ordinaryq*p slots/D=floorp/2^325/actualsmallerrorproofinquotient_certificate.md;allfiftyslots/D/paddingconstructedPAID. SelectedproductFINALcorrectionneedsSIXstoredlimbs,sinceexactcanonicalremainder<p<2^381. RemoveONLYsecondproducthighreload/finalhighsubtract(4field+2keyreducers):doesn'tchangeLOW384result. FlagsareDEADuntilnextADDS/SUBS/CMP. Fieldfirst/thirdresultthenreload10fromnextaccumulator;secondresultthen10DEADuntilninthADC-from-XZR;final10caller-temp,NOTdefinedstatus. Keyfirstbody11DEADuntilnextterm10ADC-from-XZR,secondbody11overwrittenmodulusLDP10/11orNHprefixbeforeuseonallzeropaths. Actualsix-limboutput<P/densefullporacle/STPbytes/ABIchecksremain. Obsoletephysical10zeroassertretiredbecausevalueunused,NOTbecausecanonicaloutputrequirementchanged. No freecorrection/nooriginalhighborrowcomparisonomitted.

Paid public-divisor ABI25=ACTUALDfrominitializerLDP25,xzr,[productsarg9,-16],bothD/discardedwordINITIALIZED;EIGHTcoldinstructions/ALLcontextcode/delivery2^16ADDEDinside2^124outsideR. FieldandfixedUDIV25removeheaderreloads. Products24/D25readonly,LUTroot23/size48in27/96in28/actualp29preservedABI;keypaired25/26restoreD/contextBEFOREfixedbodies. NHactualkeypointer624readctx26+112,SAMEonehotload;complete128limit96..111neveraliasesit. Fullarrays/construction/sampler/new125projectionproofinnh_certificate.md/randomness.md. No freeadvice/entropy/powers/future-targetamortization;freshscientificreviewrequired. FinitefixturesnotproductionGT/theorem/hardwarequalification.

## Uniform nonzero scaling and bounded six-reducer corrections

Production usesONEsharednonzeroFpLambda/UNSCALEDc/homogeneousrecurrence;normativeexactscaling/sampling/constructionproof in proof.md and randomness.md. No tablebias/advice/extraGTpower.

The complete shared-scaling/sign-folded8epsilon/128cap/Markov argument is normative in proof.md,section Nonzero shared scaling and complete rare-reduction bound. It applies unchanged to these actual reducers;all381bits/entrysizes retained.

The complete14-optrial/branch/19field-or18keycorrection/counter/partialabort/ABI proof is in proof.md and actual sources. Bothkernels sharefull128counter/context;fieldTBNZ/keyCBZabortENTIREattack. Actualnewfield22slot proof follows in fp2_certificate.md.


Nativeevidence:ALL16originalfullpboundary/rawtransform/48carryflags tests retainLambda=Rm,soactualunscaledoutputsunchanged. NewfieldfixtureforcesallFOURremaindersp-1/2421ops,allZEROremainders/2345ops,partialabortsafterfirstandthirdstore,counter64carry/128edges/repeatedcallswithoutreset. KeyfixturesforceEACHoddBcoordinatep-1throughM^-1=-M andboundedpublicother-coordinate search preservingfirstsign,thencheckcanonicalhashandactualzeroCBZabortwithNOmetadatawrites. Everywordagainstindependentfullporacle;scopedfinitefixturesnotproductionGT/success/physicalqualification. Oldoffsetfillprefix/signedtop/tabletests retained;102-specificbiasedCOL6prefixretiredbecauseNO biasedtableis used.


## Paired wide-column input loads with a paid saved register

EachofFOURCOLUMNWIDEusesLDP19/2offset0,LDP19/22offset16,LDP19/2offset32,insteadofoneLDPplusfourLDR. Atmiddlepair,old2remainsliveforEXTR41;new22suppliesEXTR47/UBFX18. Atlastpair,old22remainsliveforEXTR53;new2suppliesEXTR57/UBFX27. Input19/2/22areNOTclobberedbyWSTEP(only3/21loads/digittemp;accumulators4..17). Everycross-limbbit/35or34width/mask/columnstride/fullvalueunchanged. TwoREALinputinstructionsremovedperCOLUMNWIDE,fourcolumns=>EIGHT. Callee22waspreviouslypreservedimplicitly;NOWactualSTR22SP64atentry/LDR22SP64incommonEPILOGUEcostTWO,includingALLabortpaths. NetSIXactualfieldinstructionsremoved. Frame528unchanged,SP64initialized8-byteslot/disjointout48/table56/workspace80/temporary144. LR/SP/19..29exactABIchecked. Sourceworkerordinal22thereforepreservedbeforekeycall. Products24/D25/size27/28andLUT23remainreadonly. Fieldnormal2345/rawmax2421/rare19each;retain2386=2345+caller20+spare21. NewN6005Lwithallpreviousmargins/globalrare1Lunchanged;seeproof.md forexactformula/inequality. Fullrawtransform/48carryflags/forcedfourcorrections/partialcap/ABI/stack/readinitialization goldensretain.



## Ordinary column pairs reuse the already-preserved register

EachofFOURordinaryCOLUMNreadersnowusesTHREELDP19/22atinputoffset0/16/32 insteadofSIXLDR19. MSTEPsourceisanexplicitcompile-timeregister19or22,twodigits0/32perword;finalwordhighdigitwidth30stillcoversALL382bits. MSTEPwritesonly3/21tableloads/entry1/accumulators4..17/address26;input19/22andstrideconstant2remainlive. ThusALLtwelvefulltableterms/loads/lowcarry/highcarry/stridecostsareUNCHANGED;onlyTHREErealinputloadsremovedperCOLUMN,fourcolumns=>TWELVE. Extra22alreadyhasPAIDSTR/LDRSP64incommonentry/exit,NOadditionalABIwork/frame/stack/abortchange. Fieldnormal2345/rawmax2421/ret2386=2345+caller20+spare21,key806/ret827/Main720=>3933generation. Other2071=>6004L,PLUSsameglobalrare1L=>N6005L+12m+28n+2^124. ALLpriorrawtransform/fullp/48carryflags/forced4correction/partialabort/initializedreads/ordinal22ABItestsretain. Samefullvector73.875TiB/fixed13.5TiB/workerports9KiB/2^32children/control/47bitgeometry. No field/order/tablewidth/threshold/referenceprice/marginshrinks.
