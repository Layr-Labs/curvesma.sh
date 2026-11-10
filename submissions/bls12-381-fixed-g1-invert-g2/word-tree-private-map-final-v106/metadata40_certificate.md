> Current unequal nb/ng counts,ordinalmap,two verification caps and complete probability/work are in proof.md and scalar_stream_certificate.md. Native/table/metadata/sort constructions below retain their scoped contracts. FullS192 is temporary;stored Ctag24 requires verified continuation.

# Full40-byte metadata, byte occupancy and all exact transfer sizes

Metadata40 is retained exactly. Payload changes to88 bytes(U32,fullscalar32,Ctag24),with complete new two-cap probability/verification in proof.md and nh_certificate.md. Full p381/r255,GT distribution,success0.39,C20000,sixtracks/one metric are unchanged. Raw S192 is fully constructed before hashing but not stored in persistent payload. No fingerprint or ordinal bit is truncated.

## Complete layout and reversible occupancy

Metadata is exactly40 bytes:
- bytes0..19: all160 fingerprint bits,unchanged;
- byte20: baby/giant color,unchanged;
- byte21: initially zero,temporarily pending occupancy0..63;
- bytes22..23: zero;
- bytes24..39: complete little-endian128-bit ordinal.

This is U24 followed by ordinal16. U24 is the first24 bytes of the original32-byte hash/color buffer. Only eight known-zero U bytes and all sixteen trailing zero bytes of the former64-byte metadata are omitted; the pending flag no longer needs a separate word. The160-bit fingerprint,independent color and ALL128 ordinal bits remain present. Actual ordinals still satisfy0<=t<L<2^127,but no top ordinal bit is stolen for color or occupancy.

The hash initializes its upper96 bits to zero before byte20 is tagged,so byte21..23 start zero. A pending bucket's first record uses LDRB/STRB at offset21 for occupancy. Since at most63 records are pending,the byte is sufficient without overflow. Incrementing it never touches fingerprint0..19,color20 or ordinal24..39. On full completion and every final tail restore byte21 to zero BEFORE emitting the record. Every other pending record remains zero at21. Thus the final RAW40-byte row equals the original byte-for-byte.

Metadata record starts are multiples40,hence8-byte aligned. Ordinal at24 is also8-byte aligned; ordinary non-atomic AArch64 LDP24 is encodable and needs no16-byte atomic alignment assumption. All paired row copies and odd scalar final words are charged. There is no unaligned4-byte64-bit-field shortcut.

The actual metadata_pending_get/set/clear.S wrappers contain one unsigned-offset LDRB/STRB byte opcode plusRET. Production inlines that single opcode,as previously it inlined a word flag access. The test wrappers expose exact byte semantics and preservation; theirRET/caller costs are not hidden in the production one-opcode claim. Offset21 is within unsigned-byte addressing range. Set uses q<=63 and clear uses WZR. No SP/x18/callee/LR state changes.

## Source replies: five complete frames,with byte assembly before parsing

A source batch has512 full88-byte payloads=45056 bytes,EXACTLY11 full4096-byte frames. Neither88 nor40 divides4096. Root concatenates ALL11 payload frames before parsing88-byte records,as it concatenates all five metadata frames before parsing40-byte records. It now has512*40=20480 metadata bytes,EXACTLY five4096-byte frames. A frame need not contain an integral number of40-byte records. Root concatenates all FIVE frames into the complete20KiB metadata buffer before record parsing. It does not parse each frame as a separate record array,drop fragments or invent another half-batch.

The same pending flag prevents regeneration/advance between payload and metadata replies,and clears only after all five metadata frames and ACKs. All512 complete ordinals and rows are delivered in their original order. Root aggregates four512-record metadata buffers into2048 records=81920 bytes=20 frames for the ordinary bulk output. Actual production m is still a multiple2048,so all record/color boundaries align. This is a new representation,not the old erroneous four-frame reply for64-byte rows. That old defect remains acknowledged;64-byte rows needed eight frames,whereas these40-byte rows require exactly five.

Source-edgeF-framecostis512F+264+64(F-1). Payload11costs6536;metadata5costs3080,total9616/512=18.78125perrecord. Three-edgegenericbulkcost1536F+792+192(F-1). Payload11costs19608/512=38.296875;metadata20remains35160/2048=17.16796875,total55.46484375. Endpoint200/resident320remain. Generationnetwork18.78125+55.46484375+200+320=594.24609375,rounded595. Allbyteassembly/copiesremainchargedbyoldlargerreserves. Metadatasortingequationsbelowareunchanged.

All native frame routines operate on byte/word counts,not a false fixed number of records per4096-byte frame. Full source/bulk buffers are assembled before the respective record loops. Those buffers are smaller than before and fit their existing explicitly charged RAM reserves. Header offsets advance in bytes by4096,with final metadata-byte lengths20480/81920;no implicit record-sized frame step is used.

## Generic banks and exact confirmation reads

Payload banks now use ID2^94+(i>>33),local88*low33(i). All2^33 records occupy704GiB inside a1TiB bank;the trailing320GiB is NEVER read. Metadata IDs remain below2^94,payload IDs below2^95. All512/64 aligned groups fit record banks because2^33 is divisible512/64. READ11 at the final record reads exactly88 initialized bytes. Full128-bit ordinal arithmetic and independent scalar32 are retained.

Keep the same metadata DATA-ID formula buffer*2^93+(i>>34),but local byte offset is40*low34(i),computed by ordinary scalar MADD. A bank contains2^34 records occupying640GiB,strictly below its existing1TiB capacity. Unused trailing384GiB is never read as initialized data. Full64-record blocks and2048-record bulk groups cannot cross a record bank because2^34 is divisible64 and2048. At the last bank record,READ5 covers exactly its40 initialized bytes.

Confirmation fetches two READ5 metadata rows and two READ11 payloads. Five words are odd: a copy uses two LDP/STP pairs plus a final LDR/STR,so data copying costs N+1 instead of N. The actual metadata_io.S (prefix_copy5) wrapper costs7 throughRET (data6 plusRET),reads no sixth word and preserves SP/x18/callees/LR. Across three generic edges READ5 is3*(5+1)+792=810;READ11 is3*(11+1)+792=828. Total confirmation transport3276 is below the retained3408,so the original8192 complete-confirmation budget and floor(L/8192) cap remain conservative. No overread into the unused bank tail is assumed.

The payload lookup still decodes ALL128 ordinal bits at24..39. Group equality uses only first20 fingerprint bytes,so color20 and temporary flag21 are excluded. Ctag24,fullscalar32 and color are fetched from88-byte payloads. Ctag equality is probabilistic,so actual GT phase verification/continuation and a SECOND expensive-work cap are REQUIRED;see proof.md. Group events remain unchanged,but a separate rigorously bounded tagcap loss is subtracted.

## Two80-bit passes and all new message sizes

Retain the complete two-level80-bit sort construction and all initialized2^80 counters/pending buckets from radix80_certificate.md. Root and relay caches still reserve4096 bytes each,large enough for512 count words or64 scatter rows. COUNT frames remain512 words. A full SCATTER frame is now64*40=2560 bytes=320 words. A full completed reply is128-bit destination plus2560 metadata bytes=2576 bytes=322 words. Adding32-byte header fits the6144-byte half-mailbox. Initial/final singleton replies are16+40=56 bytes and are separately cold charged.

Each pending slot is63*40=2520 bytes;at2TiB all2^32 slots end11.84375TiB,disjoint from parent32TiB and child ports64..112TiB. Relay completed queues remain conservatively32MiB;4096*2576 bytes is smaller than that. Shard64*2576 remains below its reserve. Root/relay queue and finish invariants,all full128-bit destinations,bounded tail steps and no-discard stable permutation remain identical.

Prefix_copy320.S and metadata_io.S (prefix_copy322) have160/161 unrolled postincrement16 LDP/STP pairs plusRET,cost321/323. Prefix_copy512.S remains513 for count frames. All immediates are encodable. No copying relies on a large LDP immediate or odd-word overread. The actual odd5-word routine is also supplied. Caller/RET are within the existing controller budgets;data-word counts are not double priced.

Two prefix edges retain control2N+656. COUNT traffic remains1680/512=3.28125 per item;with new bulk source17.16796875 this is20.44921875,retained23. SCATTER forward is(2*320+656)/64=20.25;reply(2*322+656)/64=20.3125. With source17.16796875,ordinary WRITE320 output(3*320+792)/64=27.375 and deliberately retained pending32,total117.10546875<128. This includes all relay hops,full reply destinations,headers and the conservative pending allowance. Relay queue resident copies and every Root/relay/local record copy remain separately within the unchanged scatter-local240. Copying five record words needs six scalar instructions,not an uncharged half-pair.

Controller160 per endpoint,count-local120 and scatter-local240 are retained rather than lowered. The byte flag replaces one former word access without adding a call. All two-level cache/prefix/queue/finish operations and Root formations stay charged. Generic scan retains local128/network36 even though its actual metadata input is smaller.

## Complete original claim and evidence boundaries

Native arithmetic is unchanged:recurrence2663 retained2688,key813 retained824,NH369 retained376,scalar aggregate22 per L/main704 give generation4614. Two local count120/scatter240 passes720,scan128,cheap confirmation1,expensive verification1 and spare256 give base5720. New complete network is generation595+2*(count23+scatter128)+scan36=933. Complete hot6653 is retained6728 (75 spare).

For the SAME full L,N=6728L+2^124 satisfies EXACT N^1000<20000^1000*2^124489,proposed124489/1000=124.489. The fully constructed tables,counter/pending construction,all bounded finish/exception work remain within the earlier2^104 sort subset and total cold2^124;shrinking metadata adds no larger construction. Peak2^138/advice0 remain conservative. FixedC20000,success0.39 and the exact0.3900120659047603 success lower bound are unchanged.

The shared-only audit retains actual accepted-source full-p/key/NH checks and extends the protocol fixture to full40-byte rows,byte21 occupancy,all128 ordinal bits,source11/five-frame assembly,20-frame metadata bulk sizes,320/322 prefix messages,bank/counter boundaries and complete stable160-bit sorting. It interprets actual7/13/321/323/513 copy wrappers and2-instruction byte access wrappers against independent byte/word oracles. Sparse initialized-zero maps are scoped test fixtures;production ALL2^80 arrays and their full construction/scans are explicitly charged,not replaced by lazy dictionaries.

These checks do not execute production table/bucket construction,certify controller timing,prove GT success experimentally or establish physical hardware qualification. No new pass,acceptance,score or promotion is claimed before genuine shared review. Candidate source remained inert on the host;all execution/review/scoring belong only in shared credential-free llm-verifier.

Actual metadata_io.S (metadata_color) loads only byte20;metadata_io.S (metadata_ordinal) uses LDP24 to return ALL128 ordinal bits without writeback or overread. Their test wrappers cost2 includingRET;production inlines their single loads within existing Root formation reserves. Shared tests check both while occupancy21 is nonzero,including arbitrary full128-bit ordinal values.

For arbitrary prefix frame sizes,the complete data/controller price is N+(N mod2)+328 per edge. A singleton reply contains7 words,so two-edge traffic is2*8+656=672 rather than silently treating its odd word as a pair. Its generic exact WRITE5 costs810. Even with local/header/flag work the stated4096-per-exception allowance and total2^104 cold subset remain conservative. Hot N320/322/512 are even and every displayed hot price is unchanged.

The scalar preparation refinement in scalar_stream_certificate.md changes no metadata format or transfer price. Full32-byte payload scalars are preserved;new aggregate22L replaces the prior uniform32L generation allowance.