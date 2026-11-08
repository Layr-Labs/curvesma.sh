# Routed exact-key BSGS in GT

Track `bls12-381-fixed-g1-invert-g2`. The input is a nonidentity point `P` in G1 and an element `z` in GT. The output is `Q` in G2 with `e(P,Q)=z`. The unit is one dense Fp12 multiplication, and `C = 20000` scalar instructions. This argument has no heuristics, no randomness and no retries. Success is exactly 1. The bound is analytic. No production-size table is built.

The previous STNFS total `122.745034` is not reused. That claim had no target-relation construction, an incorrect reduction modulo `r`, and no derivation of its matrix dimensions. This submission does not claim a factor base or a sparse matrix.

## Reduction

`Q0` is the pinned nonidentity G2 generator. The pinned optimal ate pairing is nondegenerate and G1, G2 and GT are cyclic of prime order `r`. For every nonidentity `P`, `D = e(P,Q0)` therefore has order `r` and generates GT. The challenge sampling makes `z = D^a` for a unique `a` in `0..r-1`. The algorithm returns `Q = [a] Q0`. Bilinearity gives `e(P,Q) = D^a = z`. The identity `z = 1` has `a = 0` and the output is the point at infinity, which the track accepts for that target.

## Search

Let `m = h = 2^128`. Then `m*h = 2^256 > r`, since `r < 2^255`. Every `a` in `0..r-1` is uniquely `a = i*m + j` with `0 <= j < m` and `0 <= i < 2^127 < h`.

1. Reject an input that fails the track's canonical checks. Charge those checks in the fixed reserve.
2. If the twelve coefficients of `z` are the identity, return infinity.
3. Compute `D = e(P,Q0)` by one Miller loop and final exponentiation.
4. Baby list: `D^0, ..., D^{m-1}` by successive multiplication. Store each element as its twelve canonical residues.
5. Compute `D^m` by 128 squarings and `D^{-m} = (D^m)^{r-1}` by square-and-multiply.
6. Merge-sort the baby records on the exact 72-limb key.
7. Giant list: start at `z` and multiply by `D^{-m}`, once per `i = 0..h-1`. Binary-search the sorted baby keys. A hit `z * D^{-i*m} = D^j` returns `a = i*m + j`.
8. Set `Q = [a] Q0`, check `D^a = z`, and emit the canonical encoding.

Equal 72-limb keys are equal GT elements because each coefficient is the unique residue in `[0, p)`. The baby keys are distinct, so the binary search finds the unique match. Every exponent below `r` occurs, so the match exists. There is no fingerprint and no hash table.

## Field prices

An Fp element is six little-endian limbs, the unique residue in `[0, p)`. One limb product step is 4 address updates, 4 loads, `mul`, `umulh`, an 8-instruction carry chain, 2 stores and 2 branches: 22 instructions, charged 64. A reduction row is a 14-instruction quotient digit, six 16-instruction limb mul-subs and 16 instructions of carry cleanup: 126 instructions, charged 256. The final conditional subtraction is charged 256. A multiply is therefore at most

```
36*64 + 6*256 + 256 = 4096 <= 8192.
```

An addition is six limb adds (30 instructions) and a conditional subtraction (32 instructions), charged 256.

Fp2 multiplication is 4 Fp multiplications and 2 additions: `4*8192 + 2*256 = 33280`, charged `2^16`. Fp2 addition is charged `2^10`. Inversion of a nonzero Fp element by `x^(p-2)` is at most `761` multiplications (`381` squares and `380` extra multiplies), at most `761*8192 = 6230016 < 2^23`. An Fp2 inverse (norm, one Fp inverse, two multiplications) is charged `2^24`.

One complete affine G2 addition or doubling uses one Fp2 inverse, at most 4 Fp2 multiplications, 10 Fp2 additions and a `2^12` control allowance:

```
2^24 + 4*2^16 + 10*2^10 + 2^12 = 17053696 < 2^25.
```

A scalar of at most 256 bits uses at most 512 such steps, which is `2^34` instructions. The same cap covers one G1 membership test. Two membership tests are `2^35`.

One Fp12 multiplication is schoolbook: 144 Fp multiplications and 165 Fp additions for the product and the reduction `w^12 = 2*w^6 - 2`, plus stores into the twelve canonical residues.

```
144*8192 + 165*256 = 1221888 < 2^21.
```

Squaring uses the same cap. A `k`-bit exponentiation is at most `2k` Fp12 multiplications.

The Miller loop follows the pinned 64-bit loop count, so it has at most 64 iterations. One iteration is at most two G2 steps (`2^26`), a line evaluation (`2^20`) and five Fp12 multiplications (`5*2^21 < 2^24`). That is less than `2^27`, and 64 iterations are less than `2^33`. The final exponentiation is charged 1024 Fp12 multiplications (`2^31`), one Fp12 inverse (`2^30`) and 16 Frobenius maps (`2^24`). One pairing is therefore less than `2^34` and is charged `2^36`. The algorithm computes at most two pairings (the base `D`, and a second one only if the implementation checks `e(P,Q)` directly). Both are inside the fixed reserve.

## Addressing

A table of `2^128` records is not a 64-bit array. Memory is a binary tree of depth 128, with fewer than `2^129` nodes. One request is outstanding at a time. A record is 1024 bytes, which is 128 words: 72 limbs, a 256-bit index and padding. Node state is 256 bytes.

A routed transaction walks at most 256 hops and moves at most 128 words. One word on one hop is a load, a store, two address updates, a hop-counter update, a checksum add and a branch: 8 instructions, charged 32.

```
256 * 128 * 32 = 2^20.
```

The routing key is a 4-by-4 limb product (16 limb multiplications at 8 instructions, plus 64 instructions of carry) and the key comparison is 72 limb comparisons at 4 instructions, plus 64 instructions of loop control:

```
16*8 + 64 + 72*4 + 64 = 544 <= 4096.
```

Word motion plus that allowance is less than `2^21`, the transaction cap. Placement counters and record indices are four-limb integers, so a count of `2^128` does not overflow a two-word register.

Installing the node program and clearing node state is one transaction per node:

```
2^129 * 2^21 = 2^150.
```

## Work

Per baby record:

- one Fp12 multiplication, `2^21`;
- the initial routed write, `2^21`;
- 128 merge-sort passes, each one read, one 1024-instruction comparison and one write: `128 * (2*2^21 + 1024) = 2^29 + 2^17`.

The sum is `2^29 + 2^22 + 2^17 = 541196288 < 2^30`.

Per giant iteration:

- 129 binary-search probes (128 descent steps and the equality probe), each a read plus a 1024-instruction comparison;
- one Fp12 multiplication.

```
130 * 2^21 + 129 * 1024 = 272761856 < 2^29.
```

All `h` iterations are charged, including those after a hit would have returned. The fixed reserve is the sum of

| item | instructions |
| --- | ---: |
| parsing and canonical checks | `2^30` |
| two membership scalars | `2^35` |
| `z^r = 1` | `2^30` |
| two pairing caps | `2^37` |
| `D^m` | `2^28` |
| `(D^m)^(r-1)` | `2^30` |
| `[a] Q0` | `2^34` |
| `D^a = z` | `2^30` |
| output encoding | `2^24` |
| controller setup | `2^30` |

That sum is `194632482816 < 2^48`. The reserve is charged in full.

```
R = m*2^30 + h*2^29 + 2^150 + 2^48
  = 2^158 + 2^157 + 2^150 + 2^48
W = R / 20000
log2(W) = 144.301... < 145
```

`experiment.py` checks these inequalities on the integers above, including `R < 20000 * 2^145`. Preprocessing is the tree installation, the baby table and the fixed reserve. It is included in `R`, and it is also below `2^145` units. Nothing is amortized over a second target. Advice is empty. One processor is charged.

## Memory

Two record buffers are `2 * 2^128 * 1024 = 2^139` bytes. The tree is `2^129 * 256 = 2^137` bytes. The peak is `5 * 2^137`, and `log2` of that quantity is `139.321... < 150`.

## What the experiment checks

`experiment.py` recovers every exponent in the multiplicative group of the field of 101 elements, and every exponent in five additive groups of order at most 200, using the same exact-key sort and binary search. It also checks 40 products in the real base field against integer multiplication, and it checks the instruction inequalities. Those runs do not measure a production table. The production bound is the integer ledger.
