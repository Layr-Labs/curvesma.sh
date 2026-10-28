# small-key-shared-rows-finvw-tbnz

P-256 sigma rho, declared **126.741734** (C=1000). On the shared-rows revision:
`finv` is 4-bit windowed (253 squarings, 39 multiplications), saving 11,814
instructions per round; the 2^60 anchor cap test is one `tbnz x12,#60`; every
unbatched rare-path addition shares `finv`, so F_MAX 41,190 -> 29,376 and wL
2,789,697 -> 2,021,787. Step 594.599878. Organizer-image counts at 2^16 walks:
plain 39,402,364 -> 39,390,550, anchor 39,665,371 -> 39,588,022, single
40,884 -> 29,070; dumps identical. See proof.md.
