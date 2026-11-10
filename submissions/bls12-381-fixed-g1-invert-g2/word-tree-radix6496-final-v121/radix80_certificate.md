# Stable low64 then high96 radix passes with concrete node reuse

The filename radix80_certificate.md is historical. CURRENT algorithm sorts bits0..63 in a direct one-edge pass,then bits64..159 in a two-edge96-bit pass. It retains the exact160-bit fingerprint,ALL128ordinalbits,the unchanged40-byte rows/full24tag88-byte payload and actual production L/nb/ng. The exact group/hash/cap/success proof is unchanged because the final array is the SAME stable permutation by U. No parameter or hash truncation is introduced. All counter/pending arrays,including empty buckets,are constructed and scanned,not replaced with sparse production maps.

## Tree roles, complete transitions and native addresses

Keep the SAMEdepth3/fanout2^32 processor tree and all2^96 data leaves. Root children0/1 produce baby/giant streams ONLYduring GEN. GEN ends after BOTH fullproduction streams finish,ALL payload/meta frame ACKs return,pending flags clear,and final buffer0 is fully initialized. OnlyTHEN Root's children H,ALL0<=H<2^32 INCLUDINGworkers0/1,become direct low32 shards for pass0. There is no offset2+H;H maps directly tochildH and never exceeds the existing fanout. Full source tables are no longer used after GEN:all failed candidate GTpowers are dense ordinary full-field arithmetic atRoot,not recurrence LUTlookups atworkers.

Root scans buffer0 sequentially;pass0 writes buffer1. ALLpass0count/prefix/scatter/finish/queue/tail ACKs return before pass1 starts. All output positions0..L-1 are written exactlyonce,so buffer1 is complete. OnlyTHEN the same ALLRootchildren become middle relays for pass1,and ALLtheirdepth2children become low32 shards. Count/next/pending/cache arrays are reinitialized before that phase's requests. No record data bank is reclaimed:payload/metadata live atdepth3 leaves in their original banks;only now-dead cryptoLUT/code scratch RAM is reused. Pass1 readsbuffer1/writesbuffer0;final verification readsbuffer0. There is no simultaneous crypto/sort role or hidden free migration.

Root arrays:128-bit Htotal[2^32] at1TiB(64GiB),cursor64[2^32] at2TiB(32GiB),4096-byte cache[2^32] at[8,24)TiB. Cache is reused between count/scatter and both passes. Code/source/reply/stack/residentbuffers are below1TiB;childportmapping[64,112)TiB. All fit47-bit pointers. RootHcache construction/zeroing/fullprefixscans are charged cold,not free extra RAM.

For NORMALnonworkerRootchildren:parent32TiB,childports[64,112)TiB;pass0count64GiB at1TiB,pending2520*2^32 at[2,11.84375)TiB. Pass1 replaces pending with16TiB middlecache[8,24),middle128totals64GiB at1,cursors32GiB at2,completedqueue32MiB at24. Each phase's arrays andports are disjoint. Temporary aliasing BETWEENPHASES is intentional onlyafter all ACKs and before reinitialization.

WorkerRootchildren0/1 retain their GENportmappings:parent4TiB,children[80,128)TiB,cryptoA/Btables[8,62),fixed-s[64,77.5). In pass0 these two directshards put count64GiB at1TiB and pending2520*2^32 at[8,17.84375)TiB,overwriting DEADA/Btabledata. That avoids the parent4 port;their childports and fixed-sregion are disjoint. In pass1 they use SAMErelaymap as otherRootchildren:middlecache[8,24),totals1,cursors2,queue24;parent4 and children80..128remainvalid. No MMUport remapping/TLBshootdown is required. The arrays can be mapped from bootstrap,then filled when their phase begins. Initialization,including every newstore,remains charged even when afterGEN. Updating phase/base constants at one serialized boundary is paid cold.

ALLdepth2 nodes are bootstrapped as future pass1shards:count1TiB,pending[2,11.84375),parent32TiB,children[64,112). Their generic DATAforwarding always uses those mappings. A parent4 mapping would overlap pending and is explicitly NOTused. All otherordinaryinternalnodes retain parent4/children8..56. Leafbanks remain their original complete1TiB mappings. Every address is below2^47;code/MMUdescriptors/stack/dispatch/phase contexts fit4MiB perPE,large algorithmarrays separately priced. Ordinary1GiBblock descriptors suffice;no unbounded nativepointer or free globalphysicalRAM is assumed.

## Exact digit extraction and full128-bit prefixes

Let initialized metadata words0/1/2 be w0,w1,w2. Stage0(low64):H=w0>>32,B=w0&0xffffffff;the directshard handles B,there is NOmiddle relay. RootCOUNTsends B as one64-bit word with high32ZERO. Stage1(high96):H=w2&0xffffffff,M=w1>>32,B=w1&0xffffffff;RootCOUNTsends w1 torelayH,which sends B toshardM. Color begins w2bit32 and is NEVERincluded. For eachstage the concatenation H||M||B is exactly the selected full digit;stage0M is absent,not extra zero buckets. Extraction costs atmost2Root/2relay/1shardinstructions as before,stage0no relay.

Root maintains full128totals[H]and4096byteCOUNTcache[H]. Eachfull512wordframe goes directly toshardH atstage0;atstage1 relayH maintains full128totals[M],4096bytecache[M],and forwards512Bwords toshardM. Shards update128count[B]viaordinaryADD/ADC andmemoryoperations. Counters/prefixes stay belowL<2^127 but BOTH64-bit words andALL128ordinalbits remain present. Root scans ALL2^32 H in ascendingorder;stage1eachrelay scansALL2^32 M;eachshard scansALL2^32 B,reads oldcount,overwrites the SAME16bytecell withnext[B],adds oldcount tofullprefix. Stage0initializes/scans ALL2^64 counts/pendingbuckets;stage1 ALL2^96. Emptybucketsareincluded. Adjacent globalintervals have exactcountlengths and cover0..L-1.

## Stable scatter, bounded messages and saved contexts

Root caches64 complete RAW40-byte records perH in sourcearrivalorder. Stage0sends320words directly toshardH. Stage1relayH caches64records perM and forwards320words toshardM inarrivalorder. Allcachetailsareexplicitlyflushed. Shard's SAME63/64pendingalgorithm:ifemptypending andnext not64aligned,emit singletonuntilalignment(atmost63);otherwiseappendfirst63records into2520bytes,byte21storesoccupancy;64threstoresbyte21zero,assembles2560bytes,emitsdestnext-63,clearsoccupancy. Incrementnext EXACTLYonce perrecord. Finalq<=63records usefull128destnext-q+i withbyte21restored. Color20/all128ordinal24..39areuntouched.

A320wordinputbatch has64records and can produceatmost64completedmessages at a directshard. Their queue64*2576fitsprivatebufferreserve. A stage1Rootbatch mayfill64middlecaches of64records,so relayqueue<=4096messages,4096*2576<32MiB at24TiB. A completedfullreply is16bytedest+2560bytes=322words;32byteheader+reply=2608<6144half. Allheaders/data useflag0/header16/data2048physicaloffsets specified below. Fullreplydest is128bits,not64. Full64recordblocks startat64alignment andcannotcross2^34recordbanks. Singletonpayload16+40=56bytes/7words;oddwordcopyroundup is chargedcold.

Root drains each queue before its nextbatch. A shard/relay COMPLETESitsparentACKbeforewaitingforqueueGET. AfterACK it resumes genericDATAdispatch,retaining phase,cursors,caches,pendingandqueueddata in explicitstate. Root may perform genericwrites throughTHATsame nodebetweenqueueGETs;ordinaryDATAopcode is checkedfirst. OneactiveRootRPCisforwardedsequentiallydepth<=3;Root neverinterleaves unrelatedRPCs inside a stable F-frame stream. No nodeholds unfinishedACKwhilewaitingforRoot,so noqueue/recordwrite deadlock. This same explicit continuation discipline applies to worker-role nodes0/1;they do not waitforGEN afterits completed phase.

Within eachfull digit,Root/relaycache order andnext[B] preserve inputorder. Stage0therefore stablysorts exactlow64bits;stage1stablysorts exacthigh96bits while preserving low64order insideequalhigh96groups. The result is exactfull160sorting,identical to sorting the full20keybytes numerically. EqualU groups andallunorderedpairenumeration are unchanged. Noextracolor/ordinalsortkey is introduced.

## Finish, all-bucket initialization and enlarged cold/memory bounds

Finish usespersistentcursors andboundedsteps:eachshardstep scansatmost256emptybuckets oremitsatmost64tailrecords;it maypausewithin abucketq<=63. Stage1relayforwardsoneboundedstep,returnsboundedqueueandACK;Rootdrainsbeforecontinuing. Stage0usesoneedgeequivalent. Fullprefix/zeroing/finish visits EVERYallocatedcounter/pendingelement;the fixture'ssparseinitializedzeromapsare ONLYscopedtestrepresentations.

Exceptions<=126*(2^64+2^96)records forinitialmisalignment/finaltails acrossbothpasses. Even4096instructions perexception is<2^116. Include fullpendingzeroing(2520*(2^64+2^96)<2^108bytes),allcounters,fullprefixloops,Root/relaycachetails,queues,pageconstruction,boundedemptyfinish/phasecallbacks/roleinitializations;retainSORTCOLD2^120,inplaceoftheformer2^104. This deliberateincrease is required for2^96buckets and is NOTsilentlykeptattheoldallowance. For scale,atmost2^108initialbytes withwordstores/ordinaryaddressloops<2^111;all2^96prefix/finishbinsatboundedconstant<2^110;the2^116exceptionbound leavesfourbitsforremainingrouting/callbackconstants. Fullprocessor/code/MMU/bootstrap<2^121,edge/flag/controller<2^118,otherfixedfield/tables/coins/input/output<2^92;sumSTILL<2^124. Everycoldworkphase,evenafterGEN,iscountedonce.

Allalgorithm-specificpending<2^108bytes,counters<2^101,relayscaches2^76,relayqueues2^57,Root<2^45;combined<2^109. Genericdata-leafcapacity2^136,privatecode/state2^119,edges2^111stilldominate,peak<2^138. No data-leafcapacity shrinkage orphysicalVMallocationclaim. Advice0;tablesandallarraysareactuallyconstructedduringthissingleattack.

## Ordinary scalar instructions and complete traffic

Prefixcontrol160perendpoint remainsconservative:save/restore16,header8,validation12,address/context12,notification32,continuation16,pointeradvance16,counter/finish16,return8=>136<=160. OrdinaryDATAcontrol120<=128. The concrete selective-monitor construction below accountsALLwaiting,waking,dispatch,earlyheaderstores,staleevents,coldpark/STOPandstableheaders. Noexternalfree scheduler or depth-only wake argument is used.

Oneprefixedgecost N+(Nmod2)+328;twoedges2*(N+(Nmod2))+656. Actualnativecopy320/322/512costN+RETthroughRET;callerRETcoveredcontroller. Allhotframesareeven;oddexceptionrepliesroundupandarecoldcharged. StagesuseSAMEordinarygeneric20framemetadatasourcecost17.16796875/record andWRITE320cost27.375/record. Stage0COUNT840/512+17.16796875=18.80859375<=20;SCATTER648/64+650/64+17.16796875+27.375+32pending=96.82421875<=104. Stage1COUNT1680/512+17.16796875=20.44921875<=23;SCATTER1296/64+1300/64+17.16796875+27.375+32pending=117.10546875<=128. ALLdestinationwordcopies/headers/relayhops/queuedreplytrafficareincluded.

Stage0localCOUNT96 coversRootformation64+Rootdigit/buffer/totals18+shardcounter12=94;SCATTER192coversRootformation96+partition24+shardindex/counter32+pendingappend/fullformation32=184. ThereisNOrelaypartition/copy atstage0. Stage1retainsCOUNT120(Root64+18+relay18+shard12=112);SCATTER240(Root96+24+relay24+queuecopy9+shard32+pending32=217). Fivewordrecordcopyneedssixscalarinstructions,coveredformationbudgets;queuesarenotfree. Pending32isALSOretainedintrafficasconservativeoverpayment. Sumlocalsort96+192+120+240=648;sortnetwork20+104+23+128=275.

Generation4560unchanged;sort648+scan128+cheap1+expensive1+spare256=>base5594. Network595generation+275sort+36scan=906;hot6500retained6578with78spare. SAMEfullL88722074495999233395607723401466085376,N=6578L+2^124satisfiesEXACTN^1000<20000^1000*2^124457,proposed124457/1000=124.457. Fullp381/r255/uniformGT/suppliedP/success0.39/C20000/sixtracks/onemetric/exact0.3900120659047603success/advice0areunchanged.

Shared-only4096recordfixture uses direct64shards atstage0,two-level96relay/shardsatstage1,full160keys/ALL128ordinals/40byterows,512wordcount/320wordscatter/322wordreply,11+fiveframes/20bulkframes,initialmisalignment/tails,queue/padding/carry/bankboundaries. Independentfull160sort comparesEVERYrawrecord andstability. SparsemapsdoNOTclaimproductionarrayconstruction/probability/controllertiming/hardwareavailability/qualification. Existingactual2635recurrence/211fixed/570odd/fused1162/full24wordtrace/369hash/fullY/scalars/copy/180quotientboundary/notificationfixturesremain. CandidatecodeisNEVERhostexecuted/imported/compiled/emulated;sharedcredentialfreellm-verifieraloneexecutes/reviews/scores. Noactualnewpass/score/acceptance/promotionisassertedbeforegenuinereview.

# Explicit bounded selective notification and serialized dispatch

This section replaces the formerly unspecified 'finite notification premise'. It responds directly to the rejection of public submission7f7bacd8-cbb9-4cb3-b13d-607d91ef7406,commit1b1bc8ff565fc043dcabd5914fdaf7c5a1111d00. That review supported the mathematical/local fixture components but rejected the COMPLETE cost because an idle-processor wake/dispatch construction was missing. Acceptance of earlier sources is not used to fill this gap.

## Execution environment and ordinary A64 instructions

Use a bare-metal scalar A64 implementation with real WFE clock gating and ordinary shareable-memory exclusive monitors. A blocked WFE executes no further instructions until an event;it is not implemented as a polling loop or trapped hypervisor call. This is a stated choice of implementation,not a claim that every permitted implementation of the WFE hint blocks. Timer event streams,interrupt-driven scheduling,periodic timers,debug event streams and WFE traps are disabled during the attack. There is no operating-system thread scheduler. Reliable memory/no asynchronous hardware faults is the same ordinary computational reliability convention used by the field arithmetic.

The required mechanism is architectural:clearing a PE's global exclusive monitor generates an event for THAT PE,equivalent to a local event;software need not issue a global SEV. The monitor tracks one shareable reservation at a time. A peer's store to its marked reservation granule clears that reservation. We choose the ordinary address-tagged per-PE monitor implementation retaining its one reservation until a matching store/explicit clear;there is no capacity eviction table or background monitor-clear task. WFE is implemented to gate instruction issue while its event bit is zero. The one-bit event latches a matching store even when it arrives before WFE. No SEV,global event broadcast,wake-all scan,interrupt-routing oracle or atomic-retry loop is used.

These requirements are supported by Arm's documented global-monitor event behavior and a real clock-gating WFE implementation;we do NOT claim that an arbitrary emulator treating WFE asNOP meets them. References:
- Arm A-profile architecture,global-monitor clear event is local to the affected PE:https://support.arm.com/documentation/ddi0487/md/-Part-E-The-AArch32-Application-Level-Architecture/-E2-10-Synchronization-and-semaphores/-E2-10-6-Use-of-WFE-and-SEV-instructions-by-spin-locks
- Arm Cortex-X4 TRM,WFE gating and latched event behavior:https://support.arm.com/documentation/102484/0001/Power-management/Architectural-clock-gating-modes/Wait-for-Interrupt-and-Wait-for-Event?lang=en
- Arm A-profile architecture,the reservation granule is at most2048bytes:https://support.arm.com/documentation/ddi0487/mb/-Part-B-The-AArch64-Application-Level-Architecture/-Chapter-B2-The-AArch64-Application-Level-Memory-Model/-B2-12-Synchronization-and-semaphores/-B2-12-5-Load-Exclusive-and-Store-Exclusive-instruction-usage-restrictions

All LDAXR/STLR/WFE/CLREX/DMB and control instructions count as ordinary scalar A64 instructions at the fixed protected price1. Wake latency is not a wall-time speedup;the metric counts SUMMED executed instructions. Hardware monitor/event state is included in the existing per-processor private state. Mailboxes are the existing neighbor-shared buffers under the endpoints'local native mappings;this is NOT a claim of one globally addressed2^138-byte physical AArch64 RAM. The analytical distributed memory/ordinary load-store model is explicit. No extra unlisted wake primitive is priced or assumed.

## Exact12KiB edge layout,including early header wakes

Each existing tree edge has two6144-byte half-mailboxes,request and response. Their bases are multiples2048 because all12KiB edge bases/strides are multiples2048. In EACH half:
- flag64 at0..7,initial0;
- initialized padding8..15;
- header32 at16..47,two STP stores;
- unused48..2047,NEVERread/written while parked;
- data4096 at2048..6143,exactly one full frame.

The flag and header share ONE possible maximum2048-byte reservation granule. Header stores may therefore cause EARLY wakes;they are NOT ignored. There are at most TWO such header stores per published frame. Data stores touch only the other granules and cause NO flag wake. Different half-mailbox flag granules and different edges are disjoint,so unrelated idle processors' reservations are not invalidated. The layout works even at the maximum architectural granule;no small-cacheline assumption is needed. Full512-word frames and320/322-word prefix data fit. Headers and data are copied separately at these offsets;the intervening padding is not transmitted or charged as copied data. Address formation changes constants,not a record-sized or huge-width instruction.

## Actual waiter and lost-wake proof

metadata_io.S contains the actual notify_wait_flag leaf:

    wait:
        ldaxr x1,[x0]
        cbnz x1,ready
        wfe
        b wait
    ready:
        clrex
        ret

x0 identifies ONE expected flag,not an array. x1 is the only scratch register;SP/x18/callees/LR are unchanged. If a producer publishes before LDAXR,the acquire load sees1. If publication is after LDAXR and before WFE,the reserved store sets the latched local event,so WFE does not lose it. If publication is after WFE blocks,the same store releases that wait. An early header event only causes another acquire check;it cannot falsely accept a zero flag.

At most two early header stores,one real flag publication and ONE stale local event can cause empty iterations of a wait. Events coalesce in the single bit;coalescing cannot add iterations. There is no other event source in the specified execution environment. A stale event from the preceding CLREX/flag phase is at most one bit,not a new unbounded process. Thus at most FOUR empty trips,each4instructions,followed by LDAXR/CBNZ/CLREX/RET4 gives at most20 instructions per waiter call,regardless of time spent parked. A ready-at-entry call costs4. No count depends on the number2^96 of idle nodes. After a zero check the processor can remain asleep arbitrarily long without executing more instructions.

## Request/response ownership and one outstanding transaction

Flags are level bits,not64-bit sequence numbers that could overflow during this attack. Producer writes the complete header/data,executes DMB ISH then STLR1 to publish. Consumer waits with LDAXR/acquire,reads/processes data,STLR0 clears its consumed incoming flag,and publishes the response flag only AFTER clearing request. The requester consumes/clears response before sending the next request. Release/acquire ordering ensures the preceding frame is consumed and its flag cleared before reuse. No producer overwrites an outstanding frame. Clear/publication helpers in metadata_io.S are literal STLR/RET and DMB/STLR/RET.

One Root transaction is active. A forwarding node can wait only for the ONE child ACK/response belonging to its saved continuation. No unrelated parent request can arrive while that transaction is outstanding. On return it replies to its parent and resumes waiting on its parent request flag. Prefix/queue work also follows this discipline:completed blocks enter the explicit bounded queue;Root drains them with subsequent requests. There is NO unsolicited fan-in notification from all children. A node never scans2^32 child flags. Batch bursts/tails use saved state and bounded continuation requests exactly as in the routing certificate. Source GEN returns its finite batch through repeated frames to the same requester,with the pending flag preventing double generation.

This gives a well-founded depth-at-most3 RPC stack,not a cyclic wait. Each handler's arithmetic/array chunk/queue operation is finite and charged. Root may make a generic DATA write between queue GETs only after the preceding ACK,when all affected nodes have resumed ordinary parent dispatch. Neither an unfinished parent ACK nor an arbitrary asynchronous interruption is hidden. On failure/return,optional STOP traversal uses the same protocol and costs less than2^107 total,cold-paid;no free shutdown broadcast is needed.

## All notifications and dispatch fit the existing numerical allowances

Use the supplied leaf helpers literally. Per endpoint/frame:
- wait argument MOV+BL+max20 waiter =22;
- publish two argument MOV+BL+DMB/STLR/RET =6;
- clear argument MOV+BL+STLR/RET =4.
Total notification/flag/barrier/caller work is at most32 instructions. This deliberately overpays endpoints that need fewer phases. The former16-unit category is REPLACED by32,not silently retained. Existing larger controller budgets already have room:

Ordinary DATA controller:save/restore16,header8,validation12,address/context12,notification32,continuation16,opcode/routing16,return8 =120<=128. Copy data words are separate. Prefix controller:save/restore16,header8,validation12,address12,notification32,continuation16,pointer advances16,counter/finish16,return8 =136<=160. Their remaining8/24 instruction reserves include the copy-kernel call/RET and bounded opcode choices. Parent/child IDs are existing fixed three32-bit digits;address/routing operations use ordinary64/128-bit words,not arbitrary integers. No array-wide dispatch search appears in either path.

Consequently ordinary edge cost remains N+(Nmod2)+264,three-edge cost3*(N+(Nmod2))+792. Prefix edge remains N+(Nmod2)+328. Full source/bulk frame streams still have separately charged64/192 per continuation frame. All published network906 and complete cheap8192 bounds therefore remain conservative. This is an explicit justification of their notification/dispatch component,not inheritance from a source's acceptance. Data header bytes still cost their separately retained8 copy instructions;including them also in controller categories is conservative overpayment.

Each PE's one monitor/event bit and state fit4MiB private code/state. Cold initialization masks/configures the described event sources,initializes each flag/header/page and starts each PE at its parent waiter. Fewer than2^97 nodes park/configure at most64 instructions each:less than2^103. Edge initialization,all pages/code/contexts/counters and optionalSTOP remain inside the previously charged2^121 processor/2^118 edge/2^124 total cold bounds. Setup occurs BEFORE dependent frame publication;publication before a PE's first wait is safe by the acquire check. Unrequested PEs execute only this cold park and no hot polling. Memory peak2^138/advice0 and native47-bit layouts remain unchanged.

## Scoped executable fixture and exact limits

The shared-only audit interprets the ACTUAL waiter opcodes under explicit per-PE monitor/event-bit rules,enumerates producer races before load/between load and WFE/while blocked,includes both early header stores/stale events,and checks <=20 receiver instructions. It models several OTHER parked PEs on disjoint flag granules and verifies that they execute zero hot instructions. It also verifies physical12KiB half layout/header/data/flag ordering and helper opcode bodies. Its sparse small event map is a SCOPED protocol fixture;it does not allocate2^96 hardware processors or prove physical implementation availability. The mathematical arbitrary-size bound follows from unique granules,one outstanding waiter per PE,and at most three producer stores per wait,not test extrapolation.

The attack explicitly chooses a real-blocking WFE implementation with the specified ordinary monitor semantics and disabled extraneous event sources. The claim is not a portable bound for WFE-asNOP emulators or periodically interrupted systems. This limitation is visible rather than folded into 'finite notification'. No production execution,physical VM qualification or new review outcome is claimed. Fullp381/r255,success0.39/C20000/sixtracks/onecomputationmetric and complete existing bound N=6578L+2^124 are preserved.

## Stable-header continuations preserve the64-instruction streaming bound

The first frame of an F-frame stream writes its header with two STP operations and pays the full128/160 endpoint controller allowances above. Its32-byte header includes opcode,96-bitbankID,initial local offset and frame countF,all in ordinary bounded words. Header construction/validation/path choice and saved context are charged once. Every actual bulk stream stays within ONE bank:512 payload records and2048 metadata records align to their2^33/2^34 record banks. Source GEN payload11 and metadata5 are separate streams with stable source state. Root does not start an unrelated transaction until the current F-frame stream finishes.

For the remainingF-1 frames,the header is NOT REWRITTEN. Keep the current source/destination pointers,Fremaining and the known peer flag addresses in saved/live controller context. Read/write the next data frame at offset2048,clear consumed flag,and publish the next expected flag. Each endpoint already knows its pending-stream opcode/path and cursor;it does not revalidate or scan the header on continuation. The header bytes and flag0 remain in their isolated reservation granule,data alone is modified outside it. A sender issues the next level-bit request only after consuming/clearing the prior response. A forwarding node handles parent request,child request/response,copy and parent response sequentially;it never waits on all children or accepts an unrelated parent RPC halfway through this stream.

Consequently a continuation has ZERO early header stores. A waiter has at most one stale event and one true publication:at most2empty4-instruction trips plus4ready instructions=12 INCLUDING RET. With literal helper calls,the per-endpoint notification budget is waitMOV/BL+12=14,publish6,clear4,total24. Remaining streaming control is at most8:two MOV copy arguments+BL/RET4,two pointerADDs+SUBSremaining+conditionalbackbranch4. Source/destination/flag pointers live in registers19..29 across helpers that clobber only0..3;their initial context save/restore is already paid in the first controller allowance. Four distinct parent/child flag addresses can be precomputed there;no hidden6144-byteoffset reconstruction occurs each frame. Where an endpoint does not itself copy data,the same8 is conservative overpayment.

Thus each continuation endpoint is at most32,and each edge's two endpoints at most64 instructions. For a three-edge generic path the allowance is192 per continuation. One data-word copy per routing hop is counted separately asN+(Nmod2);Root resident formation/copies stay in the existing local allowances. Header-copy8 is paid ONLY on the first frame because all continuation headers remain identical. This explicitly justifies
source F frames:512F+264+64(F-1),
generic three-edge F frames:1536F+792+192(F-1),
without adding early-header wake work to an already full64 allowance. No timer,global wake or free dispatcher is concealed. Prefix320/322-word transactions are single-frame and retain their full160endpoint controller budget.

The scoped shared event fixture separately enumerates first-frame3producer-store races with max20 wait and continuation1publication races with max12 wait. A byte fixture sends full512-word streams of F=1,5,11,20 while writing the header ONCE,checks every data frame/cursor/flag clear and verifies that no continuation header byte changes. These are protocol-semantics fixtures,not physical WFE/production2^96PE certification. Arm documents the acquire/release ordering used here:https://documentation-service.arm.com/static/62d6777531ea211bb6627683 . No ordering property is inferred from the Python fixture alone.

### Exact shared mapping and32-byte DATA header

Map each shared mailbox as Normal,Shareable,WriteBack-cacheable memory,with coherent aliases at the two endpoints and caches/coherence enabled before parking. All translation descriptors are installed and synchronized BEFORE hot waits;there is no hot TLB/cache maintenance or page-fault scheduler that could clear reservations unexpectedly. The monitored wait loop fits resident code and executes no unrelated data loads while parked. These are explicit cold boot/configuration requirements of the chosen ordinary A64 implementation,not an inference about this VM or a grant of a free software scheduler.

The generic DATA header is exactly four64-bit words:IDlo64,IDhi32 paddedwithzero32,localinitialoffset64(<2^40),descriptor64. Descriptor bits0..7 are opcode,8..23 are positiveframecountF,24..33 are datawordcountN<=512,34..63 zero. ActualbulkF<=20 andN=512. Prefix-specific data remain in their payload and use their own boundedopcode;no128-bitordinal is truncated into thisheader. BankID has ALL96 bits;offset+F*4096<=2^40 for these alignedbulkstreams. Both path and bank stayconstant forFframes. Private full-word source/destination pointers and remaining-frame counter advance instead of rewriting H2/H3. Source GEN payload/meta streams use distinct initialopcodes with F11/F5 and persistent512-record batchstate;their counters do not advance a second time during frame replies.

Header formation uses two LDP/STP pairs and is included in the first-frame header/control cost. Continuations have NO new header stores/copies. This is why all early header wake counts and the64/192 continuation inequalities apply to the actual wire format rather than a conceptual oversized message. Everyunusedpadding byte is neither read nor copied.

Currentproductionkey+hashisONEretained1174call,not824+376. FullS192islogicalA96||B96readentirelywhilekey384frameislive,no192-bytecopy/materializedscratch. Publicsixargsincludehashkeys;identityreturnsbeforeNH;payloadscalar32andrawhistoriesareuntouched. Exactfunction/byteorder/lifetime/allpath/ALL24actualreadwordproofsareinkey_path_certificate.md. Generation2660+1174+22+704=4560.
