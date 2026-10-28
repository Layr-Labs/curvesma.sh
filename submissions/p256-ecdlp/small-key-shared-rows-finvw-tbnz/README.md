# small-key-shared-rows-finvw-tbnz

P-256 sigma rho, declared **126.741620** (C=1000). On the shared-rows revision:
`finv` is 4-bit windowed (253 squarings, 39 multiplications), saving 11,814
instructions per round; the 2^60 anchor cap test is one `tbnz x12,#60`; every
unbatched rare-path addition shares `finv` (F_MAX 41,190 -> 29,376, wL
2,789,697 -> 2,021,787); and A=2^26 instead of 2^24 (same whole-binary mnemonic
sequence, 122,403 instructions). Step 594.552667, preprocessing 2^83.96119.
At 2^16 walks: plain 39,402,364 -> 39,390,550, anchor 39,665,371 -> 39,588,022,
single 40,884 -> 29,070; dumps identical. See proof.md.
