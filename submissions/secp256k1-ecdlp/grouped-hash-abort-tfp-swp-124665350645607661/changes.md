# Changes against 46385d36 (read with proof.md 3, 4, 7 and success.md 6.6)

This package is the accepted 46385d36 (124.802364127839677) plus nine changes; everything else is byte-identical.
Without the 16-byte next entries and the per-group exchange row it is our 26e42f93; without items 7-9 our b6b
(a9202261, 124.67474874163155), without items 8-9 our nr1 (763558b9, 124.67441353408947), without item 9 our forward
(37a3a6b6, 124.673783083764069).

1. **Truncated-fraction leaf kernels** (proof.md 3.1-3.4; `kernels/leaf.s` bwd and keys, the fraction model in
   `experiment.py`, `machine.py` and `sites.py`): bradcush's 4598c05a design with the three-pair β cache of our
   66043fba (via bda80700), ported onto the grouped line by the crew's second climber (Claude Fable 5.1). `bwd`
   1206.5 → 1285.5 and `keys` 328.5 → 140.5 per incidence; D' 1305 residues, widening 654, box factor 1 + 2^-17, 33
   fixup sites; the rows of β (72, once per node) replace the keys' block. Delta: -0.0868 bits on 46385d36 alone.
2. **Hash join at 7 fixed per record and 5.5 per expected visit** (proof.md 4, 5; `kernels/route.s` hchg; was 11
   and 9): BFI of the bucket into x10 = 3*2^61 (check "R1 layout"), one shifted-operand ADD for the record's 16-byte
   next entry (next address, K2) at (A>>1) + C, written before the walk, the walk unrolled by two with one LDP per
   visit (5 per visit, one B after an even walk without a match: 5.5 per expected visit, `VISIT`); K2 hits (4; item 7) and listings (4, or 3 from
   the second half) stay in EK2/ELIST; the matched visit skips its CBZ or CBNZ (the run counts 4, the ledger 5).
   The next-array (2^62 bytes per node) stays inside R1 with the slots and heads (2^62.3357 of 2^62.9986).
   `grouped_paths` checks heads, next entries, candidates, cursors and every call's count against the chain model
   (toy heads at 3*2^61, 64 buckets). Delta: -0.0281 bits.
3. **Key paths E-priced under PHI, meter removed** (success.md 6.6 "Key-path means"; certificate `PJ`, `WGJ`,
   `DEFJ`, `JROW_W`; `kernels/leaf.s` JH/JT): a valid batch's leaf is uniform on G minus a point given the history,
   so P(JH or JT, valid | ·) <= PJ = 2/3 + 2 d_wj + 1/(r - 1), group g's valid key paths cost at most
   PHI[g] 4 s_g PJ, defective batches 4 per leaf globally (the crew's first climber's argument, Claude Opus 5.5);
   the two `add x18,x18,#5` are gone, so JPATH = 4, the worst case charges every leaf 4 and 6.4 has no d_j
   row. Delta: -0.0056 bits.
4. **114 groups** (proof.md 2.3, 6.6, 7): PHI, the per-group drop rows and the grouped q1 follow from the
   certificate; barriers 1943 of 2048, phase control 2^16 per node (64, 96, 114, 120 screened). The phase-4 count-exchange allowance is a per-group row under PHI at 3*2^39 per node and sub-phase
   (itemized 2^40.13; the base row keeps 2^41), the worst case charging every group. Delta: -0.0052 bits (groups),
   -7e-5 (the exchange row).
5. **Y generators at their itemization** (proof.md 7; `GEN_MULY` = 642,930, checked): no mixed addition, no X
   header; pepedesigner's p256 1772672b idea. Delta: -1.5e-8 bits.
6. **Toy witnesses** (`machine.py`): `two_group_finals` runs `h2s` and `fin` for two groups of 300 records to
   one destination (group 2's first full block rewrites group 1's partial final byte-identically; records at the
   one-record rule's addresses; counts exact); `grouped_paths` uses 64 buckets. The base-path cases of part B
   (`block_paths`, `unrolled_paths`) equal their exact models.
7. **K2-hit d4 check at 4** (`kernels/route.s` paths 10f/11f; `HCHG_VER` 5 → 4): LDR, EOR of the two w3
   words, TST vs the 46-one mask, B.NE (was LDR, two UBFX, CMP, B.NE); same accept/reject, `grouped_paths`
   re-checks exact counts. Delta: -3.35e-4 on b6b.
8. **`fwd` unrolled 8x** (`kernels/leaf.s` fwd; `FWD_INC` 98 → 97.25, `FWD_FIX` 5 → 9; captain's
   bundle 7, pre-measured -6.3e-4): 8 FWDI per main pass, 2-per-pass tail, 9 fixed;
   `fwd_count(q2) = 9 + 97 q2 + 2 floor(q2/8) + (q2 mod 8)`. Delta: -6.30e-4 on nextround.
9. **Y words in 512-byte batch entries, β pair 3 in registers** (`kernels/leaf.s`; proof.md 3.1, 3.2, 5; `BWD_INC`
   1285.5 → 1281.5): x3 addresses the Y words, so x1 is free and βu²_3, βv²_3 stay in x1, x20 (was one STP, three
   LDPs); same products and sites. Delta: -3.367e-3 on 37a3a6b6.

10. **Hash-head replace by one LSE SWP** (`kernels/route.s` hchg; `HCHG_NEW` 7 -> 6): one `swp x11,x1,[x10]` replaces the head LDR and final STR. hchg is single-threaded per
    home, so the atomic early publication is safe.

Figures: E_total 2^134.631135, claim 124.665350645607661, per leaf 244.1561, grouped q1 145307571 (success
0.3900000031), worst case 2^188.505443, memory 2^135.9989, preprocessing 2^118.0361.
