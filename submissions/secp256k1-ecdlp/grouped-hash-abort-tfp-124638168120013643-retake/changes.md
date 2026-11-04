# Changes against the repair base grouped-hash-abort-tfp-124645118129451422-repair (read with proof.md 2.1, 7)

This package is the repair base 38df975c (124.645118129451422, the route-a sender-drop repair of LEAP's
fused hop-1 arrival) plus one change, DC ZVA head zeroing; everything else is the repair's, re-evaluated.

1. **DC ZVA head zeroing** (proof.md 2.1, 7; `kernels/route.s` zero128, `kernels/a64.py` gains `dc zva`,
   certificate `ZW`, `ZERO_PASS`, `ZERO_INIT`, `machine.py` grouped_paths): the 2^57 head words were zeroed
   by 2048 post-indexed STPs per pass (2050 per 4096 words, 0.5005/word). `DC ZVA` zeroes a 64-byte block
   per base-ISA step; a pass is 16 (DC ZVA, ADD #64) -- no offset form, one ADD per block -- plus
   SUBS/B.NE: 34 per 128 words (0.2656/word) plus RET. Machine: DCZID_EL0.BS = 4, DZP = 0,
   SCTLR_EL1.DZE = 1, heads in cacheable Normal memory, base/size 64B-aligned (emulator asserts
   alignment). HB stays 57. The grouped_paths toy re-verifies the same heads/zeros (garbage pre-fill,
   one full 128-word pass, count 35, every word zero). Success, caps, groups, every route/leaf row
   untouched. Delta: -0.006950009437779 bits, the zeroing alone.

Figures: E_total 2^134.603953, claim 124.638168120013643, per leaf 239.5989, grouped q1 145307571 (success
0.3900000031), worst case 2^188.505443, memory 2^135.9989, preprocessing 2^117.1226.

History (unchanged, from the repair base): vs rejected d0a0b7ae the repair restored the per-group
Chebyshev sender-drop bound (`for _q in qs`, drop-implies-count in success.md 6.6), re-derived PHI,
E-ledger, success calc, re-priced the claim (+0.000019745593090); vs ob2 the fusion alone is
-0.020232516156239. This package adds only item 1 above.
