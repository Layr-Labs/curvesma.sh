# Explicit field-call dispatch and256 canonical key/control instructions

This is a concrete implementation of the existing formulas with fixed ordinary scalar descriptors. There is no computation oracle or executable input from the candidate author:the public formulas specify every descriptor,andthe algorithm constructs them cold from its initialized shared addresses andpublic constants. Candidate evidence itself isneverhost-executed. Thegroup/canonical/hashvalues andall failureprobabilities areunchanged.

## Fixed initialized shared buffers and programs

Use fixedshared RAWB288 bytes,CANONB288 bytes,candidatearena1152 bytes,andboundedfieldtemporaries. Before eachcanonicalizationcopycurrentrawb intoRAWB. Every descriptoroperand points into these fixed shared buffers orcoldconstants;no bankedwindow address isretained. Candidate outputs arewritten directly beforeuse;unallocated/unwrittenbytes areneverread. The existingphysicalshared2^45 allowance easily includes these buffers,programs,framesandconstants. Everydescriptor/code/constant isconstructedandchargedcoldwithin theunchanged2^120 completepriorcoldterm.

There aretwo seedprograms,selected bywhetherb0 or b2 is nonzero. Each has24 fieldcalls:four forT(x) (twoadd/sub,twoFpM),four forbeta*x/beta*y (fourFpM),and16rotationfieldoperations. The initialxcopytocandidate0 remains separatelypricedinside384minimum/copy. y,beta*x,beta*y andeightderivedpairswriteallothercandidates. The authoritative lists arefield-schedule.md andfolded-conjugation.md.

There are24 remainingprograms,selected byappliedk0..11 andselectedbpair0/2. Theywrite the other two b pairs into CANONB. Oddk hasatmost8fieldcalls (fourFpM,fouradd/sub);cheapk0/6 hasfourcopy/negationcalls;otherevenk hasatmostfourcalls. The minimumpairwasalreadywrittenintoits CANONB slotbycanonical-min.S. Thus atmost32 fieldcalls occuracrossthetwoprograms. Programlengthsarepositive fixedsmallintegers,notruntimeunknown loops.

Each32-byte descriptorstoresfourordinary64-bitwords:functionaddress,arg0,arg1,arg2. FunctionaddressesarethepinnedFpM orprovidedcanonicaladd/sub/neg/copy leaves. Encodearguments intheactualhelperABI:blstFpM uses(output,a,b),field-ops.S add/sub use(a,b,output),copy/neg use(input,output,unused). NoadapterfunctionorassumedcommonABIisneeded. Unusedarg2iszeroandneverreadbythecopy/neg leaf. Thefixeddescriptorarrayisordinaryinitialized data,notasearchedtableorfreeadvice. Its readsanddispatch arepaid.

field-program.S gives thecompleteinterpreter. TwoLDP loadthefourdescriptorwords andadvancecursor;BLRcallsthefixedfieldfunction;SUBS/B.NE advancethepositivecount. Five scalarinstructions percall,plussevenentry/exitinstructions. It preservesx19,x20,originalLR and16-byte stackalignment;allcalledleavespreservex19/x20 undertheirordinaryABI. x18untouched. NoouterLR/calleecontextislost. Atmost32callsandtwo interpreterentriescost5*32+2*7=174;twocallerinterfacescost<=6,total180. Parentcallee-preservation/frameandanyfieldwrapperworkfittheretained256 wrapperallowance. EachnativeFp416 andfieldoperation48 allowanceisstillchargedseparatelywithitsowninterface,so some overheadisconservativelycountedtwice. Nothingisfree.

## Explicit remaining key/control budget

Replace theold768key/controlallowanceby256. Thefollowingcompletebounddoesnotdiscountthefieldprogramwrapperor384minimum/copy:

- Copyrawb toRAWB:18LDP/18STP=36,plusatmost3addressinstructions,39.
- Testb0 thenb2 forzero:eachsixLDP,elevenORR,oneCBNZ=18;both36,plus4choice/identitydispatch,40.
- Selectseedprogram,pointers/count andchosenpairmode:atmost10ordinaryword/addressinstructions.
- Selectremainingprogramfrom(k,pairmode):atmost16,includingtwo-limb-free smallindex arithmetic,bases/lengthloads andbranches.
- SerializeCANONB tocurrentexactkey:fapi_copy_key37 plusfourcallerinstructions,41. Nativekeyword36remainscoldzero.
- Save/readparentcurrentrecord,kpointer,mode/scratchcontext outsidefield-programwrappers:atmost16.
- Finalkey/phase-pointerresultselection,return/failuredispatch outsideotherchargedinterfaces:atmost12.

Total39+40+10+16+41+16+12=174<256,82slack. Theminimum/copy384separatelypays327internal+16interface+14initialcopy=357. Theclasscounter16remainsinsidecontroller848. Helperprograminterpreter180andparentframepreservationareinside256wrappers. Nofieldconstruction,multiplication,rotation,keycopy,argumentmapping,pointerload,ABIworkorpartialfailedpathisunpriced. Allhelpercallsareordinaryboundedscalarcalls.

## Complete revised bound

Nonmultiplication1952=6*48structuredfieldops+16*48rotationfieldops+384minimum/copy+256key/control+256wrappers. Worsthelper10*416+1952=6112;orientationboundf(k)=4448+832c(k). TheunchangedclasscounterB=floor1027K/768 bounds
 4448K+832B<=5562K,
withthepreviouscanonical-capfailure<2^-100. Thefullordinaryledgeris
 5562+807tabulation+797jointlabelhash+46scalar+848controller
 +192head+224payload+32cappedtraversals+1retry=8509.
Retain8528K,19slack. ColdD=540*2^119,peakmemory<2^136bytes,zeroadvice andallscientificlossesareunchanged. OnereferenceperstepatC20000 gives
 W<=(28528K+540*2^119)/20000.
Theexactpublicintegerinequalitysupportingtheproposedscoreisrecordedinproof.md. ThisisanalyticconditionalevidenceunderH andtheexplicitfinitebankedRAMmodel,notproductionexecution,officialscore orqualification.

## Explicit parent wrapper bound256

Thecanonicalparentneeds onlyx19=currentrecord andx20=selectedpairmode acrosscalls. Save themandoriginalLR inanaligned32-byteframe:STPx19/x20,STRLR;restoreLDRLR,LDPx19/x20,RET. MOVx19,inputaddsone. Parententry/exitthereforeusesatmostsixinstructions;reserve16includinganyframe/addressdetails. Modeiscomputedbythealreadychargedzerotest/choiceandkisstoredbytheminimumleafatcurrentrecord+944;reload/selectuses256key/control. No thirdpersistentregisterorhiddenordinaryspill isneeded. field-program.S preservesitsownx19/x20andLR,so nestedcallsdonotdestroyparentcontext.

Twofieldprogramentriesand32fieldcallswithcallerinterfacescost180;parentcontext<=16;top-levelarguments/result/interface<=16. Sum212<256,44slack. Minimum327andits16interface+14inputcopyareseparatewithin384. Keycopy37anditsfourarguments/callareseparatewithin256key/control. Weightedclasscounter16remainsincontroller848. NativeFp416/fieldoperation48owninterfacesarealsoretained. Thus256wrapperbudgetreplaces512withallactualABI/dispatchworkexplicit. Nonmultiplication1952,base4448,worst6112,hardcanonical5562K andordinary8509retained8528follow. Exactscorecertificate(28528K+540*2^119)^1000<20000^1000*2^126193 supports126.193. No group/hashvalueorprobabilityboundchanges.
