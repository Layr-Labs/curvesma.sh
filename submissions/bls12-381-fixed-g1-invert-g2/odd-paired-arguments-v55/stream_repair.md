# Current byte and immediate contracts

Current payload is full canonical32-byte Fr scalar;metadata48 stores full183-bit fingerprint in24 bytes,color24,pending25,zero padding26..31,ALL128-bit ordinal32..47. Source512-record batch has four payload frames and SIX metadata frames;concatenate before parsing rows that straddle frames. metadata40_certificate.md,radix80_certificate.md andproof.md give the complete current routing/layout/ledger. Historical filenames do not imply40-byte rows or80-bit passes.

All36 NH key-pair loads are scalar postindexed LDP[x2],#16,so each offset is encodable and final x2 advances576 to affine coefficients. The actual shared interpreter checks initialized reads,exact full183-bit output,aliases/ABI and counts. This replaces the historical second-bank fixed-offset layout;no out-of-range LDP offset or omitted metadata frame is assumed. Native postindex does not supply bulk copying or arithmetic for free. Current exact cost is only the single proof.md ledger.
