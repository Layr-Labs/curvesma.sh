# Changes against the parent 38df975c (grouped-hash-abort-tfp-124645118129451422-repair)

kaankolcu's 38df975c (124.645118129451422, declared dependency) plus two changes in `keys`; the code is otherwise the
parent's; prose and comments are condensed for the archive cap.

1. **Hoisted emit base** (proof.md 3.2; certificate `EMIT` 7 → 6, `KEY_INC` 140.5 → 138.5, `KEY_FIX` 5 → 4,
   `KEY_SIGN` 2): w3 = d2' + (x3 << 30), with x3 fixed for a sign's three leaves. Each sign sets x18 = x3 << 30 and
   x14 = x18 | (2^36 - 1) once; each emit forms w3 by one ADD x18, W, LSR #28 (j = 0, 2) or EOR x14, S, LSR #28
   (j = 1). x3 is a multiple of 64, so the EOR equals the parent's EOR-then-ADD: records bit-identical.
   Delta -0.001716334357277 bits.
2. **Key-path tree** (proof.md 3.2, 7; success.md 6.6; `leaf.s` S1-S3, TSUM, JHB, JTB, C2, C3, TJX; certificate
   `JPATH` 4 → 3, `SJOIN` 1; `machine.py` meter; experiment TBZ count 38 → 158): a j = 1 or 2 block runs on into its
   own copy of the sign's later leaves, each path ending in one B to the sign's end; the B is charged on every sign.
   Program 13,022 instructions, within 2^16 bytes. Delta -0.001718378639432 bits.

Figures: E_total 2^134.607468, claim 124.641683416454713, per leaf 240.1834; grouped q1 145307571 (success
0.3900000031), worst case 2^188.505443, memory 2^135.9989, preprocessing 2^118.0361 unchanged. Base claim
124.894599358028. The parent's own change (fused hop-1 arrival vs ob2) is in its `changes.md`.
