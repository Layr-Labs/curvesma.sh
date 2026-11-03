# Changes against 46385d36 (read with proof.md 3, 4, 7 and success.md 6.6)

This package is the accepted 46385d36 (124.802364127839677) plus eleven changes; everything else is byte-identical.
Without the 16-byte next entries and the per-group exchange row it is our 26e42f93; without items 7-9 our b6b
(a9202261, 124.67474874163155), without items 8-9 our nr1 (763558b9, 124.67441353408947), without items 9-10 our forward
(37a3a6b6, 124.673783083764069), without item 10 our round 1 (124.670416020621392), without item 11 our round 2
(124.665350645607661).

1. **Truncated-fraction leaf kernels** (proof.md 3.1-3.4; `kernels/leaf.s`, `experiment.py`, `machine.py`, `sites.py`):
   bradcush's 4598c05a design with our 66043fba three-pair β cache (via bda80700), ported by the second climber
   (Claude Fable 5.1). `bwd` 1206.5 → 1285.5, `keys` 328.5 → 140.5 per incidence; D' 1305, widening 654,
   box 1 + 2^-17, 33 fixup sites; rows of β (72) replace the keys' block. Delta: -0.0868 on 46385d36.
2. **Hash join at 7 fixed and 5.5 per visit** (proof.md 4, 5; `kernels/route.s`; was 11 and 9): BFI into
   x10 = 3*2^61, one ADD for the 16-byte next entry at (A>>1) + C, walk unrolled by two (5 per visit, one B
   after an even walk: 5.5 expected, `VISIT`); K2 hits and listings stay in EK2/ELIST; the matched visit skips
   its CBZ/CBNZ (run counts 4, ledger 5). Next-array inside R1 (2^62.3357 of 2^62.9986). Delta: -0.0281.
3. **Key paths E-priced under PHI, meter removed** (success.md 6.6; `PJ`, `WGJ`, `DEFJ`, `JROW_W`; JH/JT):
   P(JH or JT, valid | ·) <= PJ = 2/3 + 2 d_wj + 1/(r - 1), group g costs at most PHI[g] 4 s_g PJ, defectives
   4 per leaf globally (first climber, Claude Opus 5.5); JPATH = 4, worst case charges every leaf 4, no d_j
   row. Delta: -0.0056.
4. **114 groups** (proof.md 2.3, 6.6, 7): PHI, drop rows and grouped q1 from the certificate; barriers 1943 of
   2048, phase control 2^16 (64, 96, 114, 120 screened). Exchange row under PHI at 3*2^39 (base keeps 2^41),
   worst case every group. Delta: -0.0052 (groups), -7e-5 (exchange).
5. **Y generators at itemization** (proof.md 7; `GEN_MULY` = 642,930): p256 1772672b idea. Delta: -1.5e-8.
6. **Toy witnesses** (`machine.py`): `two_group_finals` (group 2's first full block rewrites group 1's partial
   final byte-identically; counts exact); `grouped_paths` uses 64 buckets; part-B cases equal exact models.
7. **K2-hit d4 check at 4** (`kernels/route.s` paths 10f/11f; `HCHG_VER` 5 → 4): LDR, EOR of the two w3
   words, TST vs the 46-one mask, B.NE (was LDR, two UBFX, CMP, B.NE); same accept/reject, `grouped_paths`
   re-checks exact counts. Delta: -3.35e-4 on b6b.
8. **`fwd` unrolled 8x** (`kernels/leaf.s` fwd; `FWD_INC` 98 → 97.25, `FWD_FIX` 5 → 9; captain's
   bundle 7, pre-measured -6.3e-4): 8 FWDI per main pass, 2-per-pass tail, 9 fixed;
   `fwd_count(q2) = 9 + 97 q2 + 2 floor(q2/8) + (q2 mod 8)`. Delta: -6.30e-4 on nextround.
9. **Y words in 512-byte batch entries, β pair 3 in registers** (`kernels/leaf.s`; proof.md 3.1, 3.2, 5; `BWD_INC`
   1285.5 → 1281.5): x3 addresses the Y words, so x1 is free and βu²_3, βv²_3 stay in x1, x20 (was one STP, three
   LDPs); same products and sites. Delta: -3.367e-3 on 37a3a6b6.
10. **hchg head by SWP** (`kernels/route.s`, `kernels/a64.py` gains `swp`; `HCHG_NEW` 7 → 6): one SWP replaces the
   head's LDR and final STR; the walk reads no head. Delta: -5.065e-3 on item 9.
11. **DC ZVA head zeroing** (`kernels/route.s` zero128, `kernels/a64.py` gains `dc zva`; proof.md 2.1, 7):
   16 (DC ZVA, ADD #64) per pass, 34 per 128 words (was 2048 STPs, 2050 per 4096); DCZID BS = 4, DZE = 1,
   heads 64-byte aligned; HB stays 57 (58 gains 1e-4, not taken). Delta: -6.853e-3 on item 10.

Figures: E_total 2^134.624282, claim 124.658497653855674, per leaf 242.9991, grouped q1 145307571 (success
0.3900000031), worst case 2^188.505443, memory 2^135.9989, preprocessing 2^117.1226.
