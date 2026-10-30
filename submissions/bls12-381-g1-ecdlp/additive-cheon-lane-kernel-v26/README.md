# additive-cheon-lane-kernel-v26 (selected complete candidate)

Cheon's attack on bls12-381-g1-ecdlp with planned pair walks, a measured side-output kernel and
synchronous wave transport. This revision replaces the min-rule orbit key with a **parity key**,
which cuts every key extraction from 17 instructions to 10.

Claimed bound: log2 W(2^20) <= 116.322746 (exact rational 116322746/1000000,
exact total 207809621468482542913424300296884246308 instructions), heuristic success at least
0.3900 (K = 0.4952, floor 0.390050), peak memory 2^121.5854 bytes, preprocessing at most
2^115.3262 units included.

## The change (v25 → v26)

The six members of a U-orbit have y-coordinates ±y. Before, the key was the low 192 bits of
min(Y, p − Y) on Montgomery residues: a full 6-limb p − Y, a 6-limb compare, 3 `csel` and a
3-word store, 17 instructions in all.

p is odd, so exactly one of Y and p − Y is even, and that even member is also an exact orbit key.
The record keeps its **limbs 1..3** (bits 64..255). Bit 0 is always 0, so the low limb is dropped
and the key stays 192 uniform bits; routing, head index, pool fields and H3 are unchanged.

Computing it takes limbs 0..3 of p − Y in one borrow chain (borrows only move upward), one `tst`
of bit 0, 3 `csel` and a 3-word store: **10 instructions**. The same sequence is used in three
places:

- the kernel's P − T key: `side_lanes` goes from 1061 to **1054 per lane**, measured in the
  organizers' pipeline, with self-check 1 comparing every lane's key against an independent C
  parity-key computation;
- the P+T key in the pair listing: 63 → 56;
- the root key in the root listing: 60 → 53.

The toy end-to-end experiment (`evidence.py` mode `ec`) now uses the same parity rule on affine
y. Its key is the even member shifted right by one, the toy analogue of dropping the zero bit.
It still recovers every scalar correctly, with coverage events matching the exact sets.

## Ledger

Pair = 2000 + 1054 + 56; per point **1605.46** (v25 1612.46). The exact total falls by
899641898646385640234978834250172876 instructions, to 2^116.32274577, claimed at the six-decimal
ceiling 116.322746. Everything else is inherited from v25 (ed50caea), which builds on v24
(63cce2c0), v23 (34705e18) and the accepted d421a23d. No full-size execution is claimed.

Files: `proof.md`, `evidence.py`, and `measure/` (`side_lanes.S`, generator, driver, harness,
`run_measure_lanes.sh`, `reproduce_lanes.sh`, `measurement.log`).
