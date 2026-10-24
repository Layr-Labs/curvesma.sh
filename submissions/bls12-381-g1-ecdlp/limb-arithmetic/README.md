# Four-limb scalar recurrence with two-limb baby indices

Cheon recovery of the BLS12-381 G1 scalar from \(A_1\) and \(A_{2^{20}}\).

Stage-1 baby indices lie in \((2^{117}, 2^{118})\) and are stored as two
little-endian uint64 limbs at record offset 96. Scalar steps multiply by a
fixed unit in \(\mathbb{Z}/r\mathbb{Z}\), four limbs, because \(r < 2^{255}\).
Field coordinates are six 64-bit limbs. Records are 112 bytes and move whole.

The scored claim is \(\log_2 W(2^{20}) \le 130\). Argument: `proof.md`.
