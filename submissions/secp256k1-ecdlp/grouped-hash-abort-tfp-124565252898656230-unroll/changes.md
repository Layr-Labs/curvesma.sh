# Changes against the parent 38df975c (grouped-hash-abort-tfp-124645118129451422-repair)

kaankolcu's 38df975c (124.645118129451422, declared dependency) plus seventeen changes; prose condensed for the cap.

1. **Hoisted emit base** (3.2; its x14 is superseded by 16; `EMIT` 7 → 6, `KEY_INC` 140.5 → 138.5, `KEY_FIX` 5 → 4, `KEY_SIGN` 2): x18 = x3 << 30
   and x14 = x18 | (2^36 - 1) once per sign; w3 by one ADD or EOR (x3 a multiple of 64: records bit-identical).
   -0.001716334357277 bits.
2. **Key-path tree** (3.2, 7; success.md 6.6; `JPATH` 4 → 3, `SJOIN` 1): a j = 1 or 2 block runs on into its own
   copy of the sign's later leaves, one B per sign, charged on every sign. -0.001718378639432 bits.
3. **Fused hop-2 join** (2.3, 4, 7; success.md 6.6; `dlvh`, HCHGI, STP of count and p in `h2s`/`fin`; `HCHG_NEW`
   6 → 5, `DLVC`, `DLVH_FAST`, `DLVH_SLOW`, `BLOCKH`, `FINALH`, program 2^17 bytes; `grouped_paths`): `dlvh` copies
   each row block and joins its records at or past F, the slot's previous final cursor (E2 by group parity), from
   the copy's registers; no `hchg` pass. -0.005979373771831 bits.

4. **Beta by table lookup** (3.1, 3.4, 5, 7; `leaf.s` BLK0/PAIRB, `route.s` tbld; `BWD_INC` 2563/2 -> 2493/2, site
   BLK with FIX_BLK 30 C0, `TBLD`, tables 2^51 B per node): beta u^2, beta v^2 as sums of six window entries
   beta*t*2^(43i) mod p (117 vs 152 per pair). The base path (bkt, prefix4, pass4, scan_ranges, ledger R) is dropped
   (the grouped claim never used it); grouped q1 145307202. -0.030329340838130 bits.

5. **DC ZVA head zeroing** (2.1, 4, 7; `route.s` zero4096, `a64.py` DC ZVA; `ZERO_PASS` 2050 -> 1026, program then 20,262 (19,873 after 6-8)
   kernel instructions): a pass zeroes 4096 head words by 512 (DC ZVA, ADD #64), 64-byte blocks (DCZID_EL0.BS = 4,
   DZP = 0, SCTLR_EL1.DZE = 1), instead of 2048 STPs; the technique is kaankolcu's 8c3852be (declared). The zeroing
   row goes 2^128.0014 -> 2^127.0028; no other row changes. -0.007606291059698 bits.

6. **Entry STP by writeback** (4, 7; `route.s` HCHGI, dlvh; `HCHG_NEW` 5 -> 4, `DLVH_FAST` 14 -> 15, `DLVH_SLOW`
   37 -> 38): consecutive joined slots are 32 bytes apart, so their next entries (A>>1)+C are 16 apart; x13 steps by
   the STP's writeback instead of an ADD per record (one ADD per segment). -0.005305169302925 bits.
7. **Entry-address chains** (4, 5, 7; success.md 6.6; `HCHG_VISIT` 5 -> 4, `HCHG_VER` 4 -> 5, join state C + 16 and
   -2C): heads and next fields hold entry addresses, so a visit loads at the pointer itself (no ADD); a K2 hit forms
   the slot address A = 2E - 2C for its d4 check and listing. -0.002642020086674 bits.
8. **13 key sums** (3.2; `KEY_INC` 277/2 -> 269/2; `a64.py` SUB/SUBS of a shifted register): i1 = b0 - g0,
   i2 = b1 - g1, i3 = a0 - b0, i5 = a1 - b1, i6 = i1 + i2 give W(1) = (i3 + i6, i5 - i1) and W(2) = (i3 - i2,
   i5 + i6): 13 two-word ops per sign instead of 14, same values mod 2^128. -0.003569075326975 bits.

9. **One-ORR window addresses** (3.1, 7; `leaf.s` BLK0; `BWD_INC` 2493/2 -> 2477/2, site BLK 1 + 54): the table base
   0xfff8<<48 has bits 51-63 set, so `ORR x13, base, window, LSL #8` absorbs the bits above a 43-bit window and the
   AND and ADD of each window become one ORR (-4 per product). The certificate's informational ratio of analytic
   allowances to E-work is now stated below 2^-12.4 (2^-12.498; E-work fell, the allowances did not).
   -0.007164748920173 bits.
10. **Window 3 by one ORR** (3.1, 5, 7; BLK0, `route.s` tbld 7 slots, `BWD_INC` 2477/2 -> 2473/2, `TBLD` 4 + 84 +
   20*(6*2^43 + 2^41)): window 3 (bits 129..171) is limb 2 shifted left 7; its bit 0 lands on address bit 7, so window
   3's entries sit at slot 2 and again at slot 6 (+128) of each 256-byte row, and slots 2 and 3 swap (windows 3, 2).
   -0.001758990708814 bits.

11. **Out-of-line hash hits** (4, 7; success.md 6.6; `route.s` HCHGI, `a64.py` `\@` and `.pushsection`; `HCHG_END`
   1 -> 0, listing 4 in both halves, worst case 2^188.505439 -> 2^188.398524): the K2-hit blocks move after all other
   code, so a chain end falls through to the join's end after either half of the unrolled walk (4 per visit, not
   4.5). -0.001513804220626 bits.

12. **E/F constants after the first row** (3.1, 3.4; `leaf.s` CPK; -2 per pair): b⁴ or a⁴ joins the first row's five
   words (b_0*c_0 + K < 2^320) by ADDS, 3 ADCS and ADC instead of AD4's 6 after the last row; same integer to the fold.
13. **γ's V0 in (s, hi) form** (3.1, 3.2; TFINS for both γ words; -2 per pair): its low word is 8s, so keys takes it by
   the shifted operands it already uses (ADDS/SUBS ..., LSL #3); same values.
14. **BLK0's fold branches at limb 1** (3.1, 3.4; -4 per pair; site BLK 1 + 51, extra 5, bad set 30*C0*2^128): the
   common path is MUL, ADDS, ADCS, B.LO; the rare fixup finishes limbs 2 and 3 and adds C0 on their carry by CSEL.
   12-14 together: `BWD_INC` 2473/2 -> 2457/2, -0.007224502188604 bits.

15. **fwd unrolled by 32** (3.1, 3.2, 7; `leaf.s` fwd, `FWD_U`, `FWD_INC` 389/4 -> 1553/16; program 19,833 -> 22,545
   kernel instructions alone; 27,151 with 17, within 2^18 bytes): SUBS/B.HS once per 32 incidences instead of per 8 (-0.1875 per pair);
   the run adds q2 = 34 so the main loop executes. Claim 124.568590100030263 -> 124.568420519247833.
16. **EON for the j = 1 record word** (2.2, 3.2; `keys` KEY_SIGN 2 -> 1, `KEY_INC` 269/2 -> 265/2; `a64.py` EON): JH's w3 is
   EON of x18 with S >> 28, so the per-sign ORR forming x14 is gone; c for the - sign is re-coded 5 -> 4 (the - emit
   pointer moves to entry - 256) so that c in {2, 4} and the complemented {5, 3} stay distinct.
17. **h2s fully unrolled on full blocks** (4, 7; `route.s` h2s, `H2S_FULL` 6; program 2^18 bytes): a full arriving
   block runs 256 unrolled H2R with 6 fixed and no loop; other blocks pay 12 fixed plus 2 per 8-record pass.
   16-17 together: -0.003167247868309 bits; 15-17 together: -0.003337201374033 bits.

18. **Delivery processors stated (no count, row or figure changes):** 2.1 restores the parent's delivery processors (one per kind
   of message), completion and the acyclic wait graph; `h2s0` sets h2s's persistent registers once, each `h2s` run resets x1 to
   the window and loads p from a node word at 2^48 - 168 written per group by phase control (+4 within its allowance); h2s keeps
   6/12 fixed. 7 discloses the no-stop cost (about 2^124.785) and DC ZVA's fallback.

Figures: E_total 2^134.531038, claim 124.565252898656230, per leaf 227.7874; the rest as proof.md 7.
