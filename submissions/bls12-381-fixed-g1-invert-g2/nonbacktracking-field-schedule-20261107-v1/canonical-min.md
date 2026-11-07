# Pointer-based scalar minimum and tighter hard canonical cap

canonical-min.S specifies the complete leaf selecting the minimum of12 precomputed canonical Fp2 candidates. It is static evidence,neverhost-executed. Candidate k contains sigma^k of the selected nonzero b0/b2 pair,real component first48bytes,imaginary next48;each component uses six little-endian64-bit native canonical Montgomery limbs. Output96-bytepair and64-bitk destinations aredisjoint from the candidate arena. The shared arena1152bytes andcode fit the alreadycharged scratch/code allowance.

The total order is numeric real first,then numeric imaginary. View a pair as integer real*2^384+imaginary. Each comparison subtracts oldminimum-candidate in radix2^64:imaginary limbs first,then real limbs,withoneSUBS andelevenSBCS. Six pairs per operand require12 LDP loads. Loads preservecarry. The final C bit means oldminimum>=candidate;two CSEL updateonly minimum pointer andk. They preserveflags,andthe second uses no altered comparison operand. Never use final Z as multiword equality. No copying occurs on each new minimum.

One comparison costs28:ADDcandidatepointer,MOVcandidatek,twelveLDP,twelveSUBS/SBCS,twoCSEL. Eleven comparisons cost308. EntryfourMOV,finalsixLDP/sixSTP,STRk,MOVreturnk,RET add19. Internal cost327. Caller moves,BL andresultuse fit16,giving343;charge384. Allpaired offsets0,16,32,48,64,80 arealignedandinrange. Candidatebase offsets96k<=1056 fitADD. No x18,callee-savedregister,frame,internalcall,SIMD orwideunitprimitiveisused. Long-lived context ispreserved bythe existingcanonicalwrapper;its512budget isretainedseparately.

The minimum/copy allowance is384. Inadditiontothe327internalinstructionsand16interfaceinstructions,itpaystheinitialxcopyintoslot0:12LDP/STPdatainstructionsandtwoaddressinstructions,total357<384. Othercandidatefieldsarewritten directlybythefield-schedule.md product/rotation operations. Sixteenrotationfieldoperationsandninestructured-productfieldoperationsareexplicitlycountedthere;allremainingkey/control768 andwrapper512remainpaid. The minimumorderand k(sigma X)=k(X)-1 remainexact.

Thus worstcanonical7024 replaces7696,andtheorientation-dependent helper is
 f(k)=5360+832c(k), c=0fork0/6,1forothereven,2forodd.
The explicitnonmultiplicationconstantis2864,withnewbase5360;field-schedule.md givesallgeneration/copycosts. Initialize the remaining classbudget B=floor(1027K/768),ratherthan theoldfloor75K/52. B/K is4/3+1/256 up tofloorrounding,andfits128bits. The same16instructions/check arepaid. Existingstationarity/mixing stillgivevariance260K+4*2^-127K² andmean4/3+2/(3r). ExactChebyshev evaluationgives capfailure<2^-100 for thissmallerbudget. Itisnota newheuristic. Theexplicitcapbounds totalcanonicalwork
 5360K+832B<=6474K,
plusonefixedpartiallyfailedhelpercold. The counterischeckedbeforeoptionalproducts. Group-onlybudgetdecisions donotinspectcoefficienthistories.

Also replace theloosecold2^130 allowance byexplicitD_cold=540*2^119. Completeheadinitialization costs537*2^119,selectors<=32*2^92,index/setup<2^100,andallpreviouscompleteconstruction/finish<2^120. Their sum is strictly below540*2^119. This counts everyheadzero andincludesall failures,initialkey/hash/payload,code/scratch,validation,pairings,outputandverification. No pre-zeroedRAM oramortizationisassumed.

The completeordinaryledger is6474canonical+807tabulation+797jointlabelhash+46scalar+848controller+192head+224payload+32cappedtraversals+1retry=9421. Retain9440,19slack. Onereferencemultiplicationperstep atprotectedC20000 gives
 W<=(29440K+540*2^119)/20000.
Exactintegerinequality(29440K+540*2^119)^1000<20000^1000*2^126239 supportsproposed126.239,abovediagnostic126.23814349539336. Thejointphasecoupling/dictionaryvarianceandallpreviousscientificfailurelossesareunchanged,includingcanonicalcap<2^-100. Successstill>0.39 underH andtheexplicitfinitebankedstoragemodel. Noofficialscore,qualification orproductionexecutionisclaimed.
