# Finite native arithmetic, complete two-pass metadata sort and all transfers

## Accepted-source arithmetic retained

The full trace recurrence, actual constructed vector/fixed-s tables,full-width quotient argument and complete native key are retained from the public review-valid124.739 source d6ca984e-c5cf-4b6e-bfa0-ae630cfa4d0a. The divisor shift is325 bits as specified in quotient_certificate.md; no arithmetic parameter changes here. The map is2461, fused wrapper250 own gives2711 retained2752; complete key813 retained824; NH375 retained384; scalar aggregate20 per L and main704 give generation4684. All complete field sizes, original success claims and finite coins remain unchanged. No physical production execution or promotion is inherited.

The vector tables occupy[8,42)TiB on both workers and every used entry is filled. Fixed-s tables occupy[56,69.5), with all eleven35/34-bit digit tables filled. The common3728-byte constant block contains ordinary p/D and all fifty actual ordinary products0p..49p. Recurrence reads every needed0..49 slot; fixed-s only0..11. The proofs, alias/ABI and native counts remain in fp2_certificate.md,quotient_certificate.md,constant_lookup.md,key_path_certificate.md and key_contract.md.

## Complete source streams and generic storage

metadata40_certificate.md defines the full40-byte metadata:U24 (fingerprint0..19,color20,temporary occupancy21 initially zero,padding22..23),then ALL128 ordinal bits24..39. Payload remains256 bytes. No fingerprint or ordinal bit is removed. Source512 payload records128KiB use32 frames;source512 metadata records20KiB useFIVE4096-byte frames. Root concatenates all five before parsing records,so frame-boundary fragments are preserved. The pending flag clears only after all metadata ACKs. Four512-record buffers give2048 metadata records81920 bytes=20 frames. All full-size m/record/color boundaries remain aligned.

Source payload18632+metadata3080 gives21712/512=42.40625. Generic metadata bulk20 costs35160/2048=17.16796875;payload bulk32 remains109.171875. Generic output126.33984375 plus source42.40625,endpoint200 and resident320 gives688.74609375,rounded689. All byte assembly and resident copies remain charged.

Keep the depth3/fanout2^32 tree,2^96 data leaves and<2^97 processors. Metadata ID stays buffer*2^93+(i>>34),but local byte offset is40*low34(i). Each2^34-record bank occupies640GiB inside its1TiB allocation;unused tail is never read. Payload ID/stride256 are unchanged. All identifiers are ordinary96-bit DATA,not pointers,and all real addresses fit47bits. Full64-record blocks and2048-record bulk groups cannot cross a record bank.

Ordinary data price is3N+792 for even-word copies. For oddN,charge one extra scalar copy instruction per edge. Thus metadata READ5 is810,payload READ32 is888;confirmation four RPCs total3396 is below retained3408. The complete8192/cap budget stays conservative. No sixth metadata word or unused-bank byte is read.

## Two80-bit passes and exact ordinary addresses

radix80_certificate.md supplies the complete construction and controller/local/transport ledger. In brief, Root groups each80-bit digit by its high16 bits into65536 existing relays (Root children2+H). Relays group by middle32 into existing depth2 count/pending shards; those group by low32. No processor or edge count exceeds the original complete tree. All new algorithm arrays are actually initialized and scanned; no sparse/lazy production shortcut is used.

Root caches256MiB at1TiB and128-bit totals/cursors near2TiB. Relays cache16TiB at[8,24),totals64GiB@1TiB,cursors32GiB@2TiB and queue32MiB@24TiB; parent32TiB/children64..112. Count/pending shards have count/next64GiB@1TiB,pending[2,11.84375)TiB,parent32TiB/children64..112. Other internal nodes retain parent4/children8..56. The workers' table/port map is unchanged. Every range is disjoint and below2^47. Code/MMU/stack/dispatch remains4MiB per node; large algorithm arrays are separately listed as they were for the earlier Root arrays.

Only one Root transaction is active; its forwarded subrequests are sequential. Each node saves its bounded queue/cursors and returns its parent ACK before resuming ordinary DATA routing. Root drains completed queues before the next source batch and may route generic writes between queue reads. No unfinished parent ACK waits on Root output, so saved pending queues do not deadlock the generic tree. There is no free coordination or arbitrary interruption of live field arithmetic. The same paid finite notification mechanism remains.

Pass0 sorts fingerprint bits0..79 and pass1 sorts80..159. Root/relay/shard use high16/middle32/low32 extraction. Pass1 lower64=EXTR(w2,w1,16),higher16=w2bits16..31; byte20 color atw2bit32 is ignored. Input bytes0..19 alone define equality; full S192 remains in the payload. Counts and prefixes are128-bit, bounded by L<2^127. Global prefixes scan ALL high16,middle32,low32 arrays, overwriting each16-byte low count with its exact stable destination.

Root count caches512 words per H,relay512 words per M; scatter caches64 full40-byte metadata records at each level. Both levels preserve arrival order. Low shards retain first63 records per bucket and use the64th to emit aligned2560-byte metadata blocks, restoring byte21 before output. Initial unaligned records and final tails are all charged and delivered. A full aligned64-record block cannot cross a2^34-record metadata bank. The two passes are stable permutations of ALL L records, preserving full128 ordinals at24..39 and color. Generationbuffer0->pass0buffer1->pass1buffer0; final grouping scans buffer0.

## Queues, copies and bounded finish

A shard batch yields at most64 completed messages; a relay Root batch can trigger at most64 middle-cache flushes,so at most4096 messages. Maximum relay data4096*2576<32MiB is explicitly allocated; shard64*2576 fits its reserve. A full completion reply is322 words (destination128 plus2560 metadata bytes),with32-byte header total2608<6144. Root drains each relay before its next batch. Globally full blocks<=L/64, so bursts are accounted by valid total amortization rather than a false per-batch constant.

Finish is bounded and persistent: shard scans at most256 empty buckets or emits64 tail records per step,possibly pausing within one bucket; relay forwards one bounded step and Root drains it before continuing. Thus all2^80 buckets are finished without claiming their tails fit a tiny queue. Root/relay cache tails, initial alignment exceptions, all prefix/scans, initialization and empty-finish calls are included cold. radix80_certificate.md gives the explicit finite bounds and no-discard permutation proof.

The supplied prefix_copy320/322/512.S use ordinary postincrement16 pairs and cost321/323/513 throughRET. metadata_io.S (prefix_copy5) uses two pairs and a final scalar word,cost7 includingRET,to avoid metadata-bank overread. metadata_pending_get/set/clear.S expose one unsigned-offset byte opcode at21 plusRET;production inlines that opcode. All instructions are encodable and preserve SP/x18/callees/LR. Caller/RET is in the separate controls,not double priced.

## Complete per-record and cold ledger

Prefix control160 per endpoint gives2*(N+(N mod2))+656 over two edges. All hot full-frame N320/322/512 are even;odd tail frames include the extra scalar word and stay in the cold exception allowance. COUNT512 plus new metadata source is20.44921875<23. SCATTER320/322,source/output/pending is117.10546875<128. Count-local120 and scatter-local240 remain conservative,including all byte flags,odd record copies,wide counters,resident queues and Root formations.

Generation4684+two local passes720+scan128+confirmation1+spare256 gives base5789. New network689+two count23/scatter128 passes302+scan36 gives1027. Complete hot6816 is retained6896 (80 spare). Confirmation transport3396 is below old3408;Rootformation256 plus comparisons/copies768 remain below8192. All false/ignored pairs consume the same cap and at most12 phase attempts remain cold.

All processor/code/page/bootstrap work<2^121,edge/flag/controllers<2^118. New complete bucket/pending construction,prefixes,cache tails,initial/final exceptions and bounded finish are<2^104,including ALL2^80 bins. Other actual field/table/source/validation/phase work remains within its prior fixed allowance. Their sum is<2^124 cold. Algorithm-specific extra data<2^93,private119,edges111,Root47 and data-leaf capacity136 give peak<2^138. Unused generic capacity is never read as initialized data. No physical hardware qualification is asserted.

For actual full L,N=6896L+2^124 satisfies EXACT N^125<20000^125*2^15568,proposed15568/125=124.544. Full fields,target distribution,success0.39,C20000,six tracks/one metric remain unchanged. Stable sorting preserves every event of the same0.3900268131920969 lower-bound proof. Scoped shared tests cover native copies and bounded full-word packet/state models,not production array construction,controller timing,GT probability or qualification. No new pass or score is claimed before actual shared review.

Exact metadata40,source five-frame assembly,20-frame bulk,320/322 prefix packets,byte21 occupancy and exact READ5 are specified in metadata40_certificate.md. All128 ordinal bits remain separate from color/occupancy.

Scalar stream pricing now uses the actual baby7-instruction public128-bit j encoder/counter and full-r giant21-instruction addition. Including callers,retain12 per baby/28 per giant;each count<=m gives20L across all executions,including aborts. Metadata and payload fields remain full sized;see scalar_stream_certificate.md.
