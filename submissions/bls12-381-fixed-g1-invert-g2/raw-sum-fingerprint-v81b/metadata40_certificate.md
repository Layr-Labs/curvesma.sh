# Full-bit packed40 metadata with scalar-only32 payload

The production payload remains the full canonical32-byte Fr scalar. The fused native key/hash leaf still reads every conceptual S192 word and returns all183 hash bits in U32. No stored S or equality oracle is added. Every informative fingerprint pair is checked by actual GT powers and false matches continue.

## All bits and exact offsets

Bits0..182 are the complete three61-bit hash residues. Bit183 is color0/1: byte22 bit7. Byte23 contains the pending occupancy0..63 in its low6 bits; its top2 bits are zero. Bytes24..39 hold ALL128 ordinal bits in two aligned64-bit words. Thus183+1+6+128=318 meaningful bits fit40 bytes with2 reserved zero bits. No hash/scalar/ordinal bit is truncated.

The ACTUAL combined leaf now returns all183 fingerprint bits in registers:word0=x1,word1=x0,word2=x5(low55),statusx7 (zero=identity,nonnull=success). It writes no output memory. Worker ABI-preserved19=output,20=color,21/22=ordinal feed the actual inlined body STP(x1,x0,at0),ORR(x5,colorLSL55),STR(x5,at16),STP(ordinal,at24). These four instructions initialize all40 bytes;ordinalSTP was already charged in main712,so NEW sourcework is3. No argumentMOV or extraBL is hidden. The integrated actual-key fixture checkszero leaf output writes/full registerhash,andthenexecutesallfour sourceinstructionswithcallerregistersuntouched. Identityx7=0 skips packing. No hashbit is truncated.

For sorting and equality,read the first three words and MASK word2 by0x7fffffffffffff (55 low bits). This removes color and occupancy entirely. It is mandatory before upper96 radix extraction and every equality comparison. The logical192-bit key is the masked word triple,with9 high zero padding bits. Neither color nor pending contributes entropy,sorting order or equality. The native high-word helper explicitly LDRs word2 then ANDs the encodable55-one immediate. The full 183-bit collision bound and all exact group-event probabilities remain unchanged.

Every emitted row restores pending23 to zero. Color22bit7,all lower183 hash bits and both ordinal limbs remain untouched by that byte write. Stable sorting preserves all40 raw bytes except temporary occupancy restored before emission. Confirmation READ5 returns color,full ordinal and hash bits;READ4 returns the complete scalar.

## Real addresses and complete source frames

Metadata ID=buffer*2^93+(i>>34),local40*low34(i). A2^34-row bank uses640GiB within its existing1TiB allocation. Addressing is one ordinary64-bit MADD with40. Payload ID=2^94+(i>>32),local32*low32(i),128GiB per bank. IDs remain96-bit software data and pointers below2^47. Full64-row blocks cannot cross either bank. All counters/ordinals retain128 bits.

Source batches contain512 scalar payloads=16384 bytes,four complete4096-byte frames. Their512 metadata rows are20480 bytes,five complete frames. Rows may straddle frames:Root concatenates ALL five before parsing and never drops fragments. Source cannot advance/regenerate until all four/five frames and ACKs complete. Starts are multiples512,so both payload and metadata starts are4096-aligned and within their banks. Public m/n remain multiples2048 and source readiness stays within1..32.

Four metadata batches make2048 rows=81920 bytes=20 full bulk frames. Source/control/resident copies and every tail/partial row remain charged. No incomplete-frame bandwidth shortcut is used.

Source F-frame edge=512F+264+64(F-1):payload4 costs2504,metadata5 costs3080,total5584/512=10.90625 per record.

Generic three-edge bulk=3*512F+792+192(F-1). Payload4 costs7512/512=14.671875;metadata20 costs35160/2048=17.16796875. Together with source10.90625,endpoint200,resident/control320,this is562.74609375,roundedUP563. All base advances,progress,dispatch,notifications,ACKs and copies retain their ordinary scalar charges from handoff.md.

## Native copies and pending blocks

Full64 metadata rows occupy2560 bytes=320 words. A full reply adds128-bit destination:2576 bytes=322 words. Header32 gives2608<6144 half-mailbox. Singleton40-byte row plus destination16 is56 bytes=7 words;odd-word padding/copy is explicitly paid by N+(N mod2). Root/relay4096-byte caches still hold64 rows. Pending63 uses2520 bytes INSIDE the existing3024-byte slot stride;retain the full old slot allocation/zeroing and all old array addresses. No memory saving is needed for the score.

metadata_io.S has actual unrolled copy4/copy5/copy320/copy322/copy512. Even sizes use N data instructions;copy5 uses two LDP/STP pairs plus LDR/STR tail,six data instructions. All addRET for standalone fixture counts. They support exact alias or disjoint buffers,read/write exactlyN words,and preserve ABI/SP. No vector/DMA or uninitialized tail is assumed.

Two READ5 metadata and two READ4 scalars cost2*(3*(5+1)+792)+2*(3*4+792)=3228,including odd-word alignment,identical to the previous generous per-pair allowance. New color extraction/masking is included in R<=2^25. All capped confirmations remain belowL numerator and final fullG2/pairing occurs once,cold.

## Explicit supplementary scalar allowance

Retain old main712 plus EIGHT additional aggregate scalar instructions per record:worker source hash/color pack3,upper-pass COUNT/SCATTER high-word masks2,final-scan two operand masks2,spare1. These are summed across their actual processors;the name main720 is an accounting allowance,not free execution on Root. Existing local COUNT/SCATTER budgets also have slack;no reduction in those local budgets is taken. Metadata color extraction at confirmation is inside the unchanged generous per-pair bound. Ordinal LDP stays one instruction at its new aligned offset.

Native field2465/retained2506 and combined key/hash899/retained920 are unchanged. Aggregate generation/supplementary allowance2506+920+720=4146. Local sort720,scan128,confirmation1,inner spare256,network886,outer spare80 gives full non-scalar6217L. Actual12m+28n and cold2^124 remain. Current exact score inequality is in proof.md.

All full2^96 bucket sets are still constructed in both passes,including unreachable padded values. Old larger pending slots,caches,full counters,all processors/pages/queues/edges remain initialized and charged. Cold<2^124,extra data<2^110,peak<2^138,advice0 remain conservative. The bounded isolated fixture tests every frame byte,full183 fingerprint bits,color,full128 ordinals,pending restoration,masked radix digits,stable ordering and native exact copies. It is scoped evidence,not productionGT success,controller timing,hardware qualification or proof acceptance.
