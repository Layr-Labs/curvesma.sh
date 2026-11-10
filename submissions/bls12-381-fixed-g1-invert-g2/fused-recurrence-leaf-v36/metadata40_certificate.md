# Current scalar payload32 and metadata48 contract

The production payload contains ONLY the full canonical32-byte Fr scalar j or bi. Complete S192 is computed transiently in worker scratch and hashed before reuse;it is NOT stored and is NOT fetched at confirmation. Each informative fingerprint pair is checked by actualGT exponentiation for up to12 phases,continuing after false matches. All checks and cap failures are charged in proof.md. No hidden full-key or discrete-log oracle is assumed.

## Exact bytes and real addresses

Metadata is48 bytes:the full183-bit fingerprint encoded in24 bytes0..23 with9 high zero padding bits,color24,pending occupancy25 initially zero,zero26..31,and ALL128 ordinal bits32..47. NH returns U32;its last8 bytes start zero and receive color/pending/padding. Group equality compares ONLY the first24 bytes. Sorting both96-bit halves includes only the192-bit encoding and never color24 or pending25. No probability is attributed to zero padding.

Baby ordinal=j-1 and giant ordinal=m+i tag the same actually emitted point/scalar as before. Payload scalar is canonical fullFr j or bi,not a truncated counter. Stable sorting preserves every raw48-byte row. Pending occupancy25 is restored to zero before every emission. Metadata ordinals retain both64-bit limbs. Confirmation gets color and ordinal from its full READ6 row;there is no absent payload color or stored S field.

Metadata ID=buffer*2^93+(i>>34),local48*low34(i). Each2^34-row bank occupies768GiB in a1TiB allocation. Payload ID=2^94+(i>>32),local32*low32(i). Each2^32 scalar bank uses128GiB. All software IDs fit96 data bits and actual pointers remain below2^47. Scalar payload address is a64-bit shift5;metadata address is a64-bit MADD with48. Retain eight extra scalar instructions per record beyond source main704,main712,for transient S/U scratch selection and new worker/Root address/constants. All old endpoint/resident allowances remain. The transient192+32 scratch per worker is inside retained private storage and replaces larger stored payload buffers.

## Complete framed source operations

Source batches contain512 scalar payloads=16384 bytes,EXACTLY four4096-byte/512-word frames. Their512 metadata rows total24576 bytes,EXACTLY six frames. Frame boundaries may cut48-byte rows;Root concatenates ALL six metadata frames before parsing,never dropping fragments. Pending prevents regenerating/advancing until all four payload,six metadata frames and ACKs complete. Ready sequences remain within the source1..32 bound. Source m is divisible2048,so512/2048/color boundaries align.

A scalar batch starts at a multiple512 and512 divides2^32. The whole4-frame payload is4096-aligned and cannot cross its scalar bank. Four metadata batches make2048 rows=98304 bytes=24 bulk frames;all fragments are copied and delivered once. All resident/source/Root/control copies and partial tails are charged.

Source F-frame edge=512F+264+64(F-1):payload4 costs2489,metadata6 costs3656,total6160/512=12.03125 per record. Generic3-edge F bulk=3*512F+792+192(F-1):payload4 costs7512/512=14.671875,metadata24 costs42072/2048=20.54296875. Generic output35.21484375 plus source12.03125,endpoint200,resident/control320 gives567.24609375,roundedUP568. Larger metadata and every extra frame are included;this is not just payload bytes divided by bandwidth.

## Confirmation and native I/O

Each pair can use two full READ6 metadata rows and two READ4 scalar payloads. Generic even READ6=3*6+792=810;READ4=804;total3228. Fetch exact rows/scalars,no unused bank tail. Phase formulas use metadata color/full ordinal and full scalars. No S comparison occurs. The complete expensiveGT verification bound is2^29 in C*H+R,including these RPCs and all rejected phases,so all confirmations together<=L under cap=floor(L/2^29). FinalG2/pairing occurs at most once,cold paid. This supersedes the old8192/full-S confirmation scheme.

Actual metadata_io.S now has copy4/copy6/copy384/copy386/copy512,with4/6/384/386/512 data instructions plusRET. All are unrolled LDP/STP post16,disjoint or exact alias,only x0..3 clobbered,noSP/x18/callee/LR changes,no overreads. Color LDRB offset24,pending get/set/clear offset25,ordinal LDP offset32 each costone data instruction plusRET. Baby scalar helper remains7 with complete32-byte output and exact128-bit increment. All new copy sizes and byte positions are interpreted by the shared-only audit.

A full64-row block is3072 bytes/384 words;reply adds destination128=3088 bytes/386 words. Header32 gives3120<6144-byte half-mailbox. Root/relay caches4096 hold64 records;pending63 holds3024 bytes. Singleton reply16+48=64 bytes/8 words. Full64 blocks are aligned and cannot cross2^34 metadata banks. Native copies have even sizes now;old odd-copy assumptions are not reused.

The bounded fixture assembles four payload/six metadata frames and24-frame bulk,checks full183-bit/192-encoded keys,full128 ordinals,padding,color,scalar-bank carry and all bytes. A separate TEST-ONLY full-r exponent oracle checks continuing after false fingerprints,BB/equalGG skips,mixed/GG formulas and cap-before-processing. It supplies no production oracle and does not establish productionGT costs or probability. Current analytical proof and exact references supply those bounds. No candidate is host compiled/executed or run with credentials.
