# Changes vs bank aa8e1f24 grouped-hash-abort-tfp-124638168120013643-retake (read with proof.md 2.1, 3.2, 7)

This package is the bank aa8e1f24 (124.638168120013643) plus two validated ports, jointly measured; all
else is the bank's, re-evaluated. q1, success, worst case, memory unchanged.

1. **Fourfold zeroing** (pepedesigner 6c296971's validated change; proof.md 2.1, 7; `kernels/route.s`
   zero512; certificate `ZW`, `ZERO_PASS`; `machine.py` label + two-pass canary witness): four zero128
   bodies per loop (64 DC ZVA + ADD #64, SUBS, B.NE: 130 per 512 words, plus RET), same byte sequence.
2. **Hoisted emit base + key-path tree** (soispoke 45552a5c's validated changes; proof.md 3.2; `kernels/leaf.s`
   S1/S2/S3, JHB/JTB, C2/C3, TJX; certificate `EMIT`, `KEY_INC`, `KEY_FIX`, `KEY_SIGN`, `JPATH`, `SJOIN`;
   `machine.py` sign-B meter): x18 = x3<<30, x14 = x18|(2^36-1) once per sign; w3 in one ADD/EOR (emit
   7->6, keys 140.5->138.5); j = 1, 2 blocks run into per-path leaf copies ending in one B per sign
   (JPATH 4->3 + SJOIN), records bit-identical, program 11,102 instructions (44,408 bytes < 2^16).

Figures: E_total 2^134.600153, claim 124.634368311615235, per leaf 238.9687, grouped q1 145307571
(success 0.3900000031), worst case 2^188.505443, memory 2^135.9989, preprocessing 2^117.0575.
Delta vs bank -0.003799808398408; vs live lead c3e0530e 124.637604945834149: -0.003236634218914.

History (from the bank): DC ZVA zeroing replaced zero4096; the repair base restored the per-group
Chebyshev sender-drop bound (success.md 6.6) over rejected d0a0b7ae. This package adds only 1-2 above.
