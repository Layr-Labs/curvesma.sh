# Complete metadata reply and encodable NH key addressing

This revision corrects two concrete defects identified in the public review of v76 (c11d8371-e55b-41eb-b9a6-9e5d653d3835, recorded commit04126027c4dc71f7015d35c90ee1e146ad808a91). That review rejected the original full claim. Its scoped arithmetic checks did not establish native immediate encodability or complete metadata transport, and no prior rejection is treated as reversed.

## Every metadata record is transferred

A prepared batch has512 payload records of256 bytes and512 metadata records of64 bytes. Thus payload128KiB requires32 frames of512 words; metadata32KiB requires EIGHT such frames. The earlier four-frame specification was incorrect and omitted half the metadata bytes; no reconstruction or acceptance of that version is asserted.

The second source reply now transfers eight consecutive frames from the SAME initialized32KiB metadata buffer. It does not regenerate records, advance the sequence, discard another half-batch or reuse ordinals. The pending flag remains set after the payload reply and is cleared only after all eight metadata frames have completed. Both request counters remain full128-bit values. Root receives all512 metadata records, preserving their order, and aggregates four such replies for its2048-record bulk metadata output. Both color boundaries still align because m is a multiple2048. Every actual record therefore supplies its complete64-byte metadata to the stable sorting permutation.

With one-edge session price F*512+264+(F-1)*64, the eight-frame metadata reply costs4808, replacing2504. The32-frame payload still costs18632. Source total23440/512=45.78125 per record, exactly4.5 above the earlier41.28125. This counts all four extra frame handshakes, data words and pointer advances. Resident copies already charged40 words per full payload/metadata record (256+64 bytes), so this changes no existing resident-copy price. Main704 retains128 Root control instructions per item; the full session control prices remain charged separately.

Generation network becomes45.78125+136.46484375+200+320=702.24609375, rounded703. Sorting-network630 and scan36 remain unchanged, for total1369. Generation4696/base6441 includes the new quotient saving; complete hot7810 is below the retained7896 with86 spare. No clock-time substitution, missing transfer or omitted retry is used.

## Encodable scalar pair loads

The old NH listing used four LDP offsets above504. The actual corrected nh_fingerprint.S adds ONE ADD x16,x2,#288 before its NH rounds. For input pairs0..5, key LDP uses x2 with offsets0..272. For pairs6..11 it uses x16 with offsets0..272. Each key still has exactly the same effective address and independent row assignment. x16 was unused in the NH stage (input10/11,key12/13,product14/15); the later affine stage reuses it only after every NH key load has finished. There is no lost live register, added argument, callee change or hidden address operation.

The375 instructions are the earlier374 plus this ONE ADD. Caller4 gives379<=384, so the retained fingerprint budget is unchanged. Every64-bit LDP/STP offset in this source is a multiple8 in[-512,504]. ADD288 is encodable. Affine LDR offsets remain unsigned scaled multiples8 below32760. The source still has no18/SP/callee19..30 changes or BL, and exact output/input alias is safe because all S reads finish before output stores.

The shared-only audit now enforces signed scaled pair offsets[-512,504] with8-byte scaling on every LDP/STP, including writeback. Single64-bit load/store immediate bounds are also enforced (unsigned scaled0..32760, or signed unscaled[-256,255] with writeback). It checks the actual corrected375-instruction source against the independent NH128/affine oracle. This strengthens the interpreter's evidence scope; it is not a claim of a production-size run or a new pass before a genuine shared report.

## Complete unchanged threshold and honest attribution

N=7896L+2^124 still satisfies EXACT N^1000<20000^1000*2^124739, so the proposed124.739 rational bound remains supported by the repaired instruction accounting. The self-contained NH2^-63 proof, delta, exact group-event counts and final success0.3900268131920969 remain unchanged because every record is now explicitly transported and the corrected loads compute exactly the same fingerprint. Full protectedp381/r255, success0.39, C20000, six production tracks and one computation metric are preserved.

The four scoped reports described in v76/v73's public judgments are evidence for the interpreted arithmetic only. They did not prove the incorrect prior stream/encodability claims. This revision supplies actual repairs and stronger checks rather than asking for acceptance of the old versions. All candidate code remained inert on the host; execution/review/scoring belong only in shared credential-free llm-verifier. No new observed score, qualification, acceptance or promotion is asserted.
