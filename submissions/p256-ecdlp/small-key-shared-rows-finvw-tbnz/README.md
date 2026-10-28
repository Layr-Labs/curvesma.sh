# small-key-shared-rows-finvw-tbnz

P-256 sigma rho, declared **126.741792** (C=1000). Two fixed-path changes
on the shared-rows revision (126.741857): `finv` walks p-2 in 4-bit windows
(253 squarings, 39 multiplications instead of 255 and 127), saving 11,814
instructions per round; the 2^60 anchor cap test is one `tbnz x12,#60`
instead of `lsr`+`cbnz`, saving one instruction per lane per anchor round.
Step 594.623592. Measured in the organizer image at 2^16 walks: plain
39,402,364 -> 39,390,550, anchor 39,665,371 -> 39,588,022; dumps identical.
See proof.md and `evidence/measurement/finvw-tbnz.json`.
