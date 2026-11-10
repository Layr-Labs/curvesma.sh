> Current unequal nb/ng counts,ordinalmap,two verification caps and complete probability/work are in proof.md and scalar_stream_certificate.md. Native/table/metadata/sort constructions below retain their scoped contracts. FullS192 is logicalA96||B96,hashedwhilekeyframeislive;stored Ctag24 requires verified continuation.

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

The same pending flag prevents regeneration/advance between payload and metadata replies,and clears only after all five metadata frames and ACKs. All512 complete ordinals and rows are delivered in their original order. Root aggregates four512-record metadata buffers into2048 records=81920 bytes=20 frames for the ordinary bulk output. Actual production nb and ng are both multiples2048,so all record/color boundaries align. This is a new representation,not the old erroneous four-frame reply for64-byte rows. That old defect remains acknowledged;64-byte rows needed eight frames,whereas these40-byte rows require exactly five.

Source-edgeF-framecostis512F+264+64(F-1). Payload11costs6536;metadata5costs3080,total9616/512=18.78125perrecord. Three-edgegenericbulkcost1536F+792+192(F-1). Payload11costs19608/512=38.296875;metadata20remains35160/2048=17.16796875,total55.46484375. Endpoint200/resident320remain. Generationnetwork18.78125+55.46484375+200+320=594.24609375,rounded595. Allbyteassembly/copiesremainchargedbyoldlargerreserves. Metadatasortingequationsbelowareunchanged.

All native frame routines operate on byte/word counts,not a false fixed number of records per4096-byte frame. Full source/bulk buffers are assembled before the respective record loops. Those buffers are smaller than before and fit their existing explicitly charged RAM reserves. Header offsets advance in bytes by4096,with final metadata-byte lengths20480/81920;no implicit record-sized frame step is used.

## Generic banks and exact confirmation reads

Payload banks now use ID2^94+(i>>33),local88*low33(i). All2^33 records occupy704GiB inside a1TiB bank;the trailing320GiB is NEVER read. Metadata IDs remain below2^94,payload IDs below2^95. All512/64 aligned groups fit record banks because2^33 is divisible512/64. READ11 at the final record reads exactly88 initialized bytes. Full128-bit ordinal arithmetic and independent scalar32 are retained.

Keep the same metadata DATA-ID formula buffer*2^93+(i>>34),but local byte offset is40*low34(i),computed by ordinary scalar MADD. A bank contains2^34 records occupying640GiB,strictly below its existing1TiB capacity. Unused trailing384GiB is never read as initialized data. Full64-record blocks and2048-record bulk groups cannot cross a record bank because2^34 is divisible64 and2048. At the last bank record,READ5 covers exactly its40 initialized bytes.

Confirmation fetches two READ5 metadata rows and two READ11 payloads. Five words are odd: a copy uses two LDP/STP pairs plus a final LDR/STR,so data copying costs N+1 instead of N. The actual metadata_io.S (prefix_copy5) wrapper costs7 throughRET (data6 plusRET),reads no sixth word and preserves SP/x18/callees/LR. Across three generic edges READ5 is3*(5+1)+792=810;READ11 is3*(11+1)+792=828. Total confirmation transport3276 is below the retained3408,so the original8192 complete-confirmation budget and floor(L/8192) cap remain conservative. No overread into the unused bank tail is assumed.

The payload lookup still decodes ALL128 ordinal bits at24..39. Group equality uses only first20 fingerprint bytes,so color20 and temporary flag21 are excluded. Ctag24,fullscalar32 and color are fetched from88-byte payloads. Ctag equality is probabilistic,so actual GT phase verification/continuation and a SECOND expensive-work cap are REQUIRED;see proof.md. Group events remain unchanged,but a separate rigorously bounded tagcap loss is subtracted.

## Direct64 then relayed96-bit stable passes

Thehistoricalfilename radix80_certificate.md nowdescribesCURRENTlow64directshard/high96relay-shardconstruction. ALL2^64and2^96bucketarraysareinitialized/scanned,includingemptybuckets. RootchildrenHdirectly,includingworkers0/1 AFTERGENandallACKs;separateworkerphaseRAMconstructedonlyafterlastgenerationread. Full128ordinals/40byterows/byte21occupancy/pending2520/nativecopy/64alignment/bankboundariesareunchanged. Pass0buffer0->1;pass1buffer1->0;exactfull160sort.

Framesremain512COUNTwords/320SCATTERwords/322replywords plus32byteheader. OneprefixedgeN+(Nmod2)+328;twoedges2*(N+(Nmod2))+656. Stage0count20/scatter104/local96+192;stage1count23/scatter128/local120+240. Generic20framesource17.16796875/WRITE32027.375 included. Localtotal648,network275;explicitcontrollers/queues/finish/roles/addressranges/coldconstructionarefullyitemizedinradixcertificate. Nohiddenrelaycopyatdirectstage0;allrelayedstage1copiesremainpaid.

Generation4266+sort648+scan128+cheap1+expensive1+spare256=5300. Network595+275+36=906;hot6206retained6284with78spare. EXACTN=6284L+2^124satisfiesN^1000<20000^1000*2^124391,proposed124391/1000=124.391. SORTCOLDincreasesfromold2^104to2^120tochargeALL2^96buckets/tails;processor121+edge121+sort120+other92STILLtotalcold<2^124. Peak138/advice0 unchanged;algorithm-specificarrays<2^109,dominatedbyexistinggenericdata2^136. Samefullp381/r255,suppliedP/uniformGT/noauxpowers/success0.39/C20000/exactsuccess0.3900120659047603.

Shared-only4096recordfixtureusesREALdirect64shardsthen96relays/shardswithALL160bits/128ordinals/40rows,11+fiveframes/20bulk,fullcountframes,queue/carry/bank/tails/padding/stablepermutation. Sparsemapsarescopedfixtures,notproductionlazyallocationorproofsofcontrollercost/probability/hardwarequalification. Actual7/13/321/323/513copiesandbytehelpers/fullpfield/hashfixturesexecuteONLYinsharedcredentialfreeDocker. Noactualnewpass/score/acceptance/promotionclaimedbeforegenuinereview.

Actual metadata_io.S (metadata_color) loads only byte20;metadata_io.S (metadata_ordinal) uses LDP24 to return ALL128 ordinal bits without writeback or overread. Their test wrappers cost2 includingRET;production inlines their single loads within existing Root formation reserves. Shared tests check both while occupancy21 is nonzero,including arbitrary full128-bit ordinal values.

For arbitrary prefix frame sizes,the complete data/controller price is N+(N mod2)+328 per edge. A singleton reply contains7 words,so maximum two-edge traffic is2*8+656=672;direct-stage0one-edgecostis8+328=336 rather than silently treating its odd word as a pair. Its generic exact WRITE5 costs810. Even with local/header/flag work the stated4096-per-exception allowance and total2^120 sortcold subset remain conservative. Hot N320/322/512 are even and every displayed hot price is unchanged.

The scalar preparation refinement in scalar_stream_certificate.md changes no metadata format or transfer price. Full32-byte payload scalars are preserved;new aggregate20L replaces the prior uniform32L generation allowance.
Physical6144-bytehalf nowplacesheader32at16/dataframeat2048/flagat0,asexplicitlyspecifiedinradix80_certificate.md. The existingconceptualheader+datacopycountsareunchanged;unusedgapbytesarenotcopiedorread. Earlyheaderwakesareexplicitlycountedwithnotify32/endpoint.

Currentproductionkey+hashisONEretained1125call,not824+376. FullS192islogicalA96||B96readentirelywhilekey384frameislive,no192-bytecopy/materializedscratch. Publicsixargsincludehashkeys;identityreturnsbeforeNH;payloadscalar32andrawhistoriesareuntouched. Exactfunction/byteorder/lifetime/allpath/ALL24actualreadwordproofsareinkey_path_certificate.md. Generation2417+1125+20+704=4266.

CURRENTtenfull39/38digitsretainALL381bits;LUTvectors128..656TiB/fixed672..804,childports1024..1280,64KiBedgebackingwith12KiBactivehalves/NormalShareable/FEAT_LVA+LPA/64KiBgranule/52-bitlocalMMUfullyconstructedinlarge_address_certificate.md. Edgecold121(memory113)/processor121/sort120/other92sumcold124,peak138/advice0. Actualrec2392ret2417/fixed194/odd521/fused1113ret1125/gen4266/base5300/net906/hot6206ret6284spare78. NoVMhardware/productionqualificationclaimed.
