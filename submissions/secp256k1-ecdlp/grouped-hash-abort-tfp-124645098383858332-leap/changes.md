# Changes against ob2 grouped-hash-abort-tfp-124665350645607661 (read with proof.md 2.1, 4, 5, 7, 8 and success.md 6.6)

This package is ob2 (124.665350645607661) plus one change, fused hop-1 arrival with a header byte count; everything
else follows from it. Without it this is ob2.

1. **Fused hop-1 arrival** (proof.md 2.1, 4, 5, 7, 8; success.md 6.6; `kernels/route.s` h2s, H2R, fin;
   `kernels/leaf.s` FLUSH; `machine.py` part B; certificate `H2S_ARR`, `H2S_ARR_FIX`, `_route`, `d_p1`): each
   arriving column block runs `h2s` directly from the link window (the same 7-instruction H2R body: 2 LDP, BFI,
   LDADD, 2 STP, TBZ) instead of being stored to R1 slots by `dlvk` (4 per record) and re-read later, so the
   R1 slot loop, its 2 STP store, its 2 LDP load and its per-slot set-up are gone (saving ~4.25 per record, partly
   offset by 11 per arrival block and the count words below). Every block header gains a byte-count third word
   (32 per record, 8192 when full), stored just before the header STP in the keys FLUSH, H2R and `fin` send paths
   (one STR each: `FLUSH` 8 → 9, `FIN_SENT` 3 → 4); `h2s` loads it (LDR) and derives its end (ADD), 11 per block
   (`H2S_ARR_FIX`: 3 MOV, LDR, ADD, 3 TBZ, CMP, B.HS, RET); `dlvk` never reads it. R1 slots, E1 finals and hop-1
   RC cursors are retained (same layout, same `fin` writes, still charged) as a conservative control; only their
   reads are gone. With no R1 writes the hop-1 slot-overflow event is empty (`d_p1 = 0`, empty loop), and the
   grouped q1 falls to base (145307571 → 145305376). Part B holds arriving column blocks and runs `h2s` on each
   (byte count from the header, row blocks byte-exact, counts exact, including 1..19/24/44/256-record and empty
   prologue cases). Delta: -0.020252261749329 bits (-0.02023 the fusion, -0.00002 the same lever's q1 side-effect).

Figures: E_total 2^134.610883, claim 124.645098383858332, per leaf 240.7526, grouped q1 145305376 (success
0.3900009466), worst case 2^188.505421, memory 2^135.9989, preprocessing 2^118.0361.
