> Current unequal nb/ng counts,ordinalmap,two verification caps and complete probability/work are in proof.md and scalar_stream_certificate.md. Native/table/metadata/sort constructions below retain their scoped contracts. FullS192 is logicalA96||B96,hashedwhilekeyframeislive;stored Ctag24 requires verified continuation.

# Two full80-bit stable radix passes with explicit two-level routing

This changes only the metadata sort implementation. The supplied-P reduction, complete production nb/ng/L, actual trace recurrence, fused1164 key/hash retained1176,verified-tag continuation, finite coins, all group events and success lower bound remain unchanged. Every metadata byte is retained. Two LSD stable passes process fingerprint bits0..79 and80..159, never color byte20. This is a finite sum of ordinary64-bit scalar work, not an elapsed-time or parallel speed claim.

## Existing tree, exact addresses and complete arrays

Keep the existing depth3/fanout2^32 ordinary processor tree and all2^96 data leaves. Root children0/1 remain the baby/giant workers. Root children2+H,0<=H<2^16, serve as relays. Each relay's child M,0<=M<2^32, serves as a counter/pending shard. These are existing depth2 processors, not extra unpriced processors. They continue ordinary storage routing to their depth3 leaves. The remaining tree retains its existing routing roles.

Root algorithm arrays are a256MiB cache (2^16*4096 bytes) at1TiB,1MiB of128-bit H totals at2TiB and512KiB of64-bit cursors immediately above those totals. The same cache is reused across count/scatter phases. Source/reply buffers and code/stack/dispatch are below1TiB; Root child ports remain[64,112)TiB. All ranges are disjoint.

Each selected relay has16TiB of middle caches (2^32*4096) at[8,24)TiB,64GiB of128-bit middle totals at1TiB,32GiB of cursors at2TiB and a32MiB completed-message queue at24TiB. Its parent edge is at32TiB and child ports at[64,112)TiB. Code/stack/dispatch is below1TiB. Each selected depth2 shard has a64GiB count/next array (2^32*16) at1TiB and2520-byte pending slots for ALL2^32 low digits at[2,11.84375)TiB; parent32TiB and child ports[64,112)TiB avoid every algorithm array. A parent at4TiB would overlap pending RAM and is explicitly NOT used for these nodes.

Other internal nodes retain parent4TiB/children[8,56)TiB. Workers retain parent4TiB, vector tables[8,62), fixed-s tables[64,77.5) and child ports[80,128)TiB exclusive. Every native pointer is strictly below2^47. Changing a node's port-base context changes a constant value, not an arithmetic instruction width. Ordinary DATA dispatch remains the first opcode branch, so its existing hot routing price is unchanged; new prefix opcodes use the larger controller allowance below.

Code, MMU descriptors, stack and dispatch state remain within the retained4MiB per-processor bootstrap budget. The LARGE ALGORITHM ARRAYS and queues listed above are separately constructed data, as were the original Root count/pending arrays; they are not falsely asserted to fit4MiB. Maps can use the same ordinary1GiB block scheme, with all descriptors, code and source/controller setup charged. No large address, DATA identifier or array element is a unit-cost arbitrary-width pointer.

## Full-word digit extraction and128-bit prefixes

Represent an80-bit digit as H16,M32,B32. If U words are w0,w1,w2, pass0 has low64=w0,H=w1&65535. Pass1 has low64=EXTR(w2,w1,16),H=UBFX(w2,16,16). In either pass M=high32(low64),B=low32(low64). Thus Root needs at most two extraction instructions, relay at most two and low shard at most one. Pass1 excludes word2 bit32, where color byte20 begins. No160-bit native arithmetic or color-dependent comparison is used.

Root caches512 packed low64 words per H during COUNT. Each full4096-byte frame is sent to relay2+H. A relay increments its128-bit total[M], caches512 B words per M, then forwards a full4096-byte frame to shard M. Each shard increments count[B] as ordinary low/high64-bit words. Root likewise maintains128-bit total[H]. All counters are below L<2^127; ordinary ADD/ADC pairs and all memory accesses are priced.

ALL2^80 bucket counters and pending slots are initialized. At prefix time Root scans all2^16 H totals in ascending order. Each relay scans ALL2^32 M totals in ascending order and supplies each shard its global128-bit starting prefix. Each shard scans ALL2^32 B counts, reads the old count, overwrites that same16-byte cell with next[B], then adds the old count to its running128-bit prefix. A single count/next array suffices. Empty buckets are not skipped in the production construction. The sparse maps used by the bounded test are only initialized-zero experiment fixtures, not the production algorithm.

The resulting global intervals are ordered by H,M,B, exactly numeric80-bit digit order. Each interval length is its exact count, and their sum is L. All frame lengths, indexes, DATA masks and wide counters are ordinary64-bit or paired64-bit data. Metadata IDs remain buffer*2^93+(i>>34), local40*low34(i); payload IDs remain2^94+(i>>32), local256*low32(i). All IDs fit96 data bits and all actual pointers fit47 bits.

## Stable scatter, full messages and bounded queues

Root caches64 complete RAW40-byte metadata records per H in arrival order. A full320-word frame is sent to that relay. The relay caches64 records per M, preserving their order, and forwards a full320-word frame to the corresponding shard. Root and relay final cache tails are explicitly flushed; none are discarded. Input order within every full digit is therefore unchanged.

Each shard uses the original63/64 pending algorithm with exact128-bit destinations. When pending is empty and next[B] is not64-record aligned, emit a singleton and advance next until aligned, at most63 such records. Thereafter store first63 records in the2520-byte slot; the first record's padding byte21 stores occupancy. On the64th arrival restore that byte to zero, assemble exactly2560 metadata bytes, emit destination next-63 and advance next. Clear occupancy. All128 ordinal bits at24..39 and color byte20 are untouched. Final pending q<=63 records use next-q+i and restore padding before emission.

Every full block starts at a multiple64, so cannot cross a2^34-record metadata bank because2^34 is divisible64. The final full sort is a stable permutation of every original RAW40-byte record. Two passes preserve increasing lower80 bits inside the final upper80-bit intervals, giving exact full160-bit sorting. Generation starts buffer0, pass0 writes1, pass1 writes0; the final group scan explicitly reads buffer0.

A shard receives at most64 incoming records per batch, so emits at most64 completed messages. One relay Root batch can fill at most64 middle caches, each forwarding64 records. Therefore its maximum queue is4096 completed messages, each at most2576 bytes: under32MiB. Each shard queue is at most64*2576 bytes, below its private buffer reserve. A full reply is322 words:128-bit destination followed by2560 metadata bytes. Adding the32-byte header gives2608 bytes, below the existing6144-byte half-mailbox. No oversized implicit frame is used.

Root drains a relay's queue before sending its next batch. A queued full block is copied to the relay's explicit queue; this additional copy is priced below. After each ACK both relay and shard resume ordinary DATA dispatch with their queue/cursor state saved. Root may perform generic destination writes through any routing node between queue reads, even when that node has saved algorithm data. No node waits for Root to drain its queue while holding an unfinished parent ACK, and no live arithmetic is arbitrarily interrupted.

There is one active Root transaction, with sequential forwarded subrequests along its path as in the existing multi-edge RPC. Root does not start an unrelated RPC before that transaction returns. This is not simultaneous free coordination or a parallel scheduler. Every forward, return, flag, saved continuation and copied word is included. A batch can have a burst of completions, but the TOTAL number of full blocks is at most L/64; the traffic bound is global amortization, not a false per-batch constant-time bound.

## Bounded finish and complete cold work

Root flushes H caches, relays flush M caches, and shards flush pending records. Finish uses persistent cursors and bounded steps: a shard step scans at most256 empty buckets or emits at most64 tail records, and may pause within one bucket's q<=63 tail. Relays forward one bounded step and return a bounded queue; Root drains it before continuing. Thus finish does not materialize all2^32 tails into a supposedly tiny queue. Empty scans, callbacks, Root control and all initial/tail exceptions are charged.

Across each pass, initial/tail exceptions are at most126*2^80 records. Two passes, even at4096 instructions per exceptional operation plus routing/finish overhead, remain below the retained2^104 cold subset. This also pays ALL count/pending zero construction, Root/relay cache tails, full prefixes over2^80 counters, queue resets, all page stores and bounded empty-finish calls. For scale, pending data is below2^92 bytes, prefix scans below a constant times2^80, and exceptional-record traffic below2^101 scalar instructions. The2^104 allowance is deliberately conservative and is included, not added outside the complete bound.

Processor code/state/bootstrap remains below2^121 total; edge/flag/controller initialization below2^118; source R and all other fixed field/key/phase work below2^92 (with table fills within their earlier allowance). Adding the new2^104 subset still gives complete cold<2^124. No work is amortized over another target or supplied as advice.

All selected shards' pending data totals<2^92 bytes, their counters2^84, relay caches2^60, relay queues2^41 and Root data<2^47. Complete algorithm-specific extra data is<2^93, negligible beside generic data-leaf capacity2^136, private code/state2^119 and edge pages2^111. Peak remains<2^138 bytes. These are analytical bounds, not an assertion that the review VM physically allocates these arrays.

## Scalar copy and controller accounting

New prefix control is allowed160 instructions at each endpoint, rather than the original128. A conservative fixed-body ledger includes save/restore16, header8, validation12, address/context12, flag/wake/barrier32, continuation16, bounded pointer advances16, counter/finish16 and return8:136<=160. The spare24 covers new opcode/queue selection and both caller/RET instructions of the copy kernels. These controls are separate from per-record processing. The explicit selective-monitor waiter,early-header-wake bound32/endpoint and ordinary120<=128/prefix136<=160 dispatch construction appended below REPLACE the missing finite-notification premise. No globalSEV or unboundedpolling occurs.

One prefix edge costs N+(N mod2)+328 (two160 controls+eight header-copy instructions), and two edges cost2*(N+(N mod2))+656. For hot full frames N=320,322,512 are even,so the exact formulas reduce to N+328 and2N+656. Odd singleton/tail frames include the extra scalar final-word copy and are separately cold charged. Actual native metadata_io.S (prefix_copy320)/metadata_io.S (prefix_copy322) supply data copying with160/161 unrolled LDP/STP postincrement16 pairs, costing N+RET1 throughRET. Post16 is encodable; no LDP offset exceeds504. Caller/RET is already within the control allowance, not double counted. The routines alter only caller-saved x0..3, preserve SP/x18/callees/LR, and support disjoint or exact-alias buffers.

COUNT full-frame prefix traffic stays3.28125 per record. With new ordinary20-frame source read17.16796875 it is20.44921875,retained23. SCATTER forward320-word traffic is20.25 per record and322-word reply traffic20.3125. With source17.16796875,WRITE320 output27.375 and retained pending32,total117.10546875<128. All reply destinations,headers,relay hops and resident copies are charged;partial tails stay cold. Native copy odd5-word semantics and exact source/bank sizes are in metadata40_certificate.md.

Count local120 includes Root formation64, Root digit/buffer/totals<=18, relay digit/buffer/totals<=18 and shard index/counter<=12:112<=120. Scatter local240 includes Root formation96, Root partition<=24, relay partition<=24, relay queue copy<=9 amortized, shard indexing/counter<=32 and pending append/full formation<=32 amortized:217<=240. The pending32 is also retained in the network allowance as conservative overpayment. Every intermediate resident copy is included; no free Root repack, relay queue or wide arithmetic is assumed.

These are scalar instruction-equivalent upper bounds and sum all processors' work. They are not CPU cycle or wall-clock timings. The exact native copy counts and bounded packet-state checks do not replace this analytical controller/total-work argument.

## Complete original claim and scoped evidence

Generation4566 is unchanged. Two local count120/scatter240 passes720,scan128,cheap_confirmation1,expensive_verification1 and spare256 give base5672. Generation network595+two count23/scatter128 passes302+scan36 gives933. Hot6605 is retained6680 (75 spare). For full L,N=6680L+2^124 and EXACT N^1000<20000^1000*2^124479,proposed124479/1000=124.479. Full parameters,success0.39,C20000,six tracks/one metric are unchanged.

The stable permutation preserves all exact keys, fingerprint groups and ordinals. The verified continuation,phase recovery,two caps and finite event proofs follow proof.md at its current nb/ng counts. Hence the same exact rational success bound0.3900120659047603 applies; there is no new hash, Poisson, random-walk or record-independence premise.

The shared-only audit retains all accepted-source arithmetic tests and adds actual7/15/321/323/513 copy and byte-access interpretation plus a bounded4096-record packet/state model using full160-bit fingerprints,128-bit ordinals, actual40/88-byte records,11 payload/five metadata frames per source batch, queue limits, stable sorting, padding restoration and counter/bank boundaries. The model's sparse initialized-zero dictionaries are explicit test fixtures; it does not execute production2^80 array construction or certify dispatcher timing, success probability, physical hardware or qualification. No new pass, score, acceptance or promotion is claimed before actual shared review. Candidate code was never host executed or compiled.

Metadata representation is now40 bytes U24+full128-bit ordinal16,with temporary occupancy at padding byte21. metadata40_certificate.md gives every exact byte size,source five-frame assembly,20-frame bulk,bank/confirmation bound and native byte/copy opcode. Its layout definitions apply throughout this topology.
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

Consequently ordinary edge cost remains N+(Nmod2)+264,three-edge cost3*(N+(Nmod2))+792. Prefix edge remains N+(Nmod2)+328. Full source/bulk frame streams still have separately charged64/192 per continuation frame. All published network933 and complete cheap8192 bounds therefore remain conservative. This is an explicit justification of their notification/dispatch component,not inheritance from a source's acceptance. Data header bytes still cost their separately retained8 copy instructions;including them also in controller categories is conservative overpayment.

Each PE's one monitor/event bit and state fit4MiB private code/state. Cold initialization masks/configures the described event sources,initializes each flag/header/page and starts each PE at its parent waiter. Fewer than2^97 nodes park/configure at most64 instructions each:less than2^103. Edge initialization,all pages/code/contexts/counters and optionalSTOP remain inside the previously charged2^121 processor/2^118 edge/2^124 total cold bounds. Setup occurs BEFORE dependent frame publication;publication before a PE's first wait is safe by the acquire check. Unrequested PEs execute only this cold park and no hot polling. Memory peak2^138/advice0 and native47-bit layouts remain unchanged.

## Scoped executable fixture and exact limits

The shared-only audit interprets the ACTUAL waiter opcodes under explicit per-PE monitor/event-bit rules,enumerates producer races before load/between load and WFE/while blocked,includes both early header stores/stale events,and checks <=20 receiver instructions. It models several OTHER parked PEs on disjoint flag granules and verifies that they execute zero hot instructions. It also verifies physical12KiB half layout/header/data/flag ordering and helper opcode bodies. Its sparse small event map is a SCOPED protocol fixture;it does not allocate2^96 hardware processors or prove physical implementation availability. The mathematical arbitrary-size bound follows from unique granules,one outstanding waiter per PE,and at most three producer stores per wait,not test extrapolation.

The attack explicitly chooses a real-blocking WFE implementation with the specified ordinary monitor semantics and disabled extraneous event sources. The claim is not a portable bound for WFE-asNOP emulators or periodically interrupted systems. This limitation is visible rather than folded into 'finite notification'. No production execution,physical VM qualification or new review outcome is claimed. Fullp381/r255,success0.39/C20000/sixtracks/onecomputationmetric and complete existing bound N=6680L+2^124 are preserved.

## Stable-header continuations preserve the64-instruction streaming bound

The first frame of an F-frame stream writes its header with two STP operations and pays the full128/160 endpoint controller allowances above. Its32-byte header includes opcode,96-bitbankID,initial local offset and frame countF,all in ordinary bounded words. Header construction/validation/path choice and saved context are charged once. Every actual bulk stream stays within ONE bank:512 payload records and2048 metadata records align to their2^33/2^34 record banks. Source GEN payload11 and metadata5 are separate streams with stable source state. Root does not start an unrelated transaction until the current F-frame stream finishes.

For the remainingF-1 frames,the header is NOT REWRITTEN. Keep the current source/destination pointers,Fremaining and the known peer flag addresses in saved/live controller context. Read/write the next data frame at offset2048,clear consumed flag,and publish the next expected flag. Each endpoint already knows its pending-stream opcode/path and cursor;it does not revalidate or scan the header on continuation. The header bytes and flag0 remain in their isolated reservation granule,data alone is modified outside it. A sender issues the next level-bit request only after consuming/clearing the prior response. A forwarding node handles parent request,child request/response,copy and parent response sequentially;it never waits on all children or accepts an unrelated parent RPC halfway through this stream.

Consequently a continuation has ZERO early header stores. A waiter has at most one stale event and one true publication:at most2empty4-instruction trips plus4ready instructions=12 INCLUDING RET. With literal helper calls,the per-endpoint notification budget is waitMOV/BL+12=14,publish6,clear4,total24. Remaining streaming control is at most8:two MOV copy arguments+BL/RET4,two pointerADDs+SUBSremaining+conditionalbackbranch4. Source/destination/flag pointers live in registers19..29 across helpers that clobber only0..3;their initial context save/restore is already paid in the first controller allowance. Four distinct parent/child flag addresses can be precomputed there;no hidden6144-byteoffset reconstruction occurs each frame. Where an endpoint does not itself copy data,the same8 is conservative overpayment.

Thus each continuation endpoint is at most32,and each edge's two endpoints at most64 instructions. For a three-edge generic path the allowance is192 per continuation. One data-word copy per routing hop is counted separately asN+(Nmod2);Root resident formation/copies stay in the existing local allowances. Header-copy8 is paid ONLY on the first frame because all continuation headers remain identical. This explicitly justifies
source F frames:512F+264+64(F-1),
generic three-edge F frames:1536F+792+192(F-1),
without adding early-header wake work to an already full64 allowance. No timer,global wake or free dispatcher is concealed. Prefix320/322-word transactions are single-frame and retain their full160endpoint controller budget.

The scoped shared event fixture separately enumerates first-frame3producer-store races with max20 wait and continuation1publication races with max12 wait. A byte fixture sends full512-word streams of F=1,5,11,20 while writing the header ONCE,checks every data frame/cursor/flag clear and verifies that no continuation header byte changes. These are protocol-semantics fixtures,not physical WFE/production2^96PE certification. Arm documents the acquire/release ordering used here:https://documentation-service.arm.com/static/62d6777531ea212bb6627683 . No ordering property is inferred from the Python fixture alone.

### Exact shared mapping and32-byte DATA header

Map each shared mailbox as Normal,Shareable,WriteBack-cacheable memory,with coherent aliases at the two endpoints and caches/coherence enabled before parking. All translation descriptors are installed and synchronized BEFORE hot waits;there is no hot TLB/cache maintenance or page-fault scheduler that could clear reservations unexpectedly. The monitored wait loop fits resident code and executes no unrelated data loads while parked. These are explicit cold boot/configuration requirements of the chosen ordinary A64 implementation,not an inference about this VM or a grant of a free software scheduler.

The generic DATA header is exactly four64-bit words:IDlo64,IDhi32 paddedwithzero32,localinitialoffset64(<2^40),descriptor64. Descriptor bits0..7 are opcode,8..23 are positiveframecountF,24..33 are datawordcountN<=512,34..63 zero. ActualbulkF<=20 andN=512. Prefix-specific data remain in their payload and use their own boundedopcode;no128-bitordinal is truncated into thisheader. BankID has ALL96 bits;offset+F*4096<=2^40 for these alignedbulkstreams. Both path and bank stayconstant forFframes. Private full-word source/destination pointers and remaining-frame counter advance instead of rewriting H2/H3. Source GEN payload/meta streams use distinct initialopcodes with F11/F5 and persistent512-record batchstate;their counters do not advance a second time during frame replies.

Header formation uses two LDP/STP pairs and is included in the first-frame header/control cost. Continuations have NO new header stores/copies. This is why all early header wake counts and the64/192 continuation inequalities apply to the actual wire format rather than a conceptual oversized message. Everyunusedpadding byte is neither read nor copied.

Currentproductionkey+hashisONEretained1176call,not824+376. FullS192islogicalA96||B96readentirelywhilekey384frameislive,no192-bytecopy/materializedscratch. Publicsixargsincludehashkeys;identityreturnsbeforeNH;payloadscalar32andrawhistoriesareuntouched. Exactfunction/byteorder/lifetime/allpath/ALL24actualreadwordproofsareinkey_path_certificate.md. Generation2664+1176+22+704=4566.
