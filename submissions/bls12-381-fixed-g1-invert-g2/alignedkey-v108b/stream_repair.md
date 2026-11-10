# Current source framing and record preservation

Payload32 is the full canonicalFr scalar. Packedmetadata40 retains ALL183 hashbits,colorbit183,pendingbyte23 andALL128 ordinalbits24..39. Source512 batches have four payload and FIVE metadata frames. Concatenate every metadata frame before parsing40-byte rows;rows may straddle frames. No regeneration/advance until all data and ACKs finish. Four batches make2048 rows/twenty bulk frames. Stable sorting masks the high metadata word to55 bits before upper96 digit/equality,excluding flags. Pending23 is cleared before emission. Full byte/address/native/count/price contracts are in metadata40_certificate.md andradix80_certificate.md;protected cryptographic parameters/threshold/referenceprices are unchanged.


