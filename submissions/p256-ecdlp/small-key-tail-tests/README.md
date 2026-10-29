# small-key-tail-tests (windowed inversion, one-instruction cap, rare-path repricing)

P-256 sigma rho, declared **126.735049** (C=1000), on the top-limb distinguished-test
package (126.735172). Three changes: `finv` is 4-bit windowed (253 squarings, 39
multiplications), saving 11,814 instructions per round; the 2^60 anchor cap is one
`tbnz x12,#60`; every unbatched rare-path addition shares `finv`, so F_MAX 41,190 ->
29,376 and wL 2,789,700 -> 2,021,790. PLAIN_FIXED 605,581,542, C_COMMON 577.577735,
total step 591.850860. Organizer-image paired counts at 2^16 walks: plain 39,222,158 ->
39,210,344, anchor 39,485,165 -> 39,407,816. The referenced branchy.json checkpoints are
in the base submissions; this package omits that file for the 256 KiB limit.
