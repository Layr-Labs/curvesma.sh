# Changes against the parent 38df975c (grouped-hash-abort-tfp-124645118129451422-repair)

kaankolcu's 38df975c (124.645118129451422, declared dependency) plus four changes; prose condensed for the cap.

1. **Hoisted emit base** (3.2; `EMIT` 7 → 6, `KEY_INC` 140.5 → 138.5, `KEY_FIX` 5 → 4, `KEY_SIGN` 2): x18 = x3 << 30
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

Figures: E_total 2^134.571159, claim 124.605374701844752, per leaf 234.2111; the rest as proof.md 7.
