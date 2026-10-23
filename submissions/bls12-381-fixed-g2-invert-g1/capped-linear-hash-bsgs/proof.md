# Capped linear-hash baby-step giant-step

Analytic classical attack for `bls12-381-fixed-g2-invert-g1`. No production-size run is claimed. The norm is `classical-group-work-v2`: one dense Fp12 multiplication costs one reference unit, and every other scalar AArch64 instruction-equivalent costs `1/C` with the track constant `C = 20000`. A reference multiplication is not also expanded into its internal instructions.

## Parameters

Let `p` and `r` be the field modulus and subgroup order in the track profile. Let

```text
N = ceil(391 * r / 1000)
m = 143186686509166558328661515423460682036
```

`m` is an integer satisfying `(m-1)^2 < N ≤ m^2`, `m^2 ≤ r`, and `m < 2^127`. The search uses `h = m` baby exponents and `m` giant exponents. Let `B = 2^84` and `T = B * 2^193 = 2^277`.

## Reduction

The challenge gives a nonidentity `Q` in the order-`r` subgroup of G2 and `z` in GT. `r` is prime, so the order of `Q` is `r`. Fix the profile generator `P0` of G1 and the pinned pairing, which the profile realises as `py_ecc` `pairing(Q, P)` with final exponentiation. Nondegeneracy on generators gives `e(P0, Q0) ≠ 1`. GT has prime order `r`, so `g0 = e(P0, Q0)` generates GT. The challenge samples `Q = [b]Q0` and `z = g0^{a b}` with `b ∈ {1,…,r-1}` and `a ∈ {0,…,r-1}` uniform and independent, which is the same law as `Q` uniform on nonidentity r-torsion and `z` uniform on GT.

Set `g = e(P0, Q) = g0^b`. Then `g` generates GT, and `z = g^a` for that same `a`. Returning `P = [a]P0` yields `e(P, Q) = g^a = z` by bilinearity. `a = 0` returns the point at infinity and is the unique solution of `z = 1`. No other `Q` is substituted.

## Canonical key

A GT element is twelve coefficients in `[0, p)`, degree order. Each coefficient is the unique six-limb little-endian base-`2^64` integer (`6 * 64 = 384 ≥ 381`). The key is these 72 limbs. Two GT elements are equal if and only if the 72 limbs are equal. Byte length of one key is `72 * 8 = 576`.

## Linear hash

Draw 72 seeds `s_0,…,s_71` independently and uniformly from `{0,1,…,2^256-1}`. These are fresh algorithm coins; generating and storing them is charged below. For a key `(k_0,…,k_71)` define

```text
F(k) = (sum_{i=0}^{71} k_i * s_i) mod 2^256
```

with integer products. `F` is deterministic given the seeds. Equal keys have equal fingerprints.

**Lemma.** If two keys differ, then over the seeds `P(F(k) = F(k')) ≤ 2^{-193}`.

**Proof.** Let `d_i = k_i - k'_i ∈ Z`, so `|d_i| < 2^64` and some `d_j ≠ 0`. Equality of fingerprints is `sum_i d_i s_i ≡ 0 (mod 2^256)`. Let `t = v_2(d_j)`. Then `0 ≤ t ≤ 63` and `d_j = 2^t * u` with `u` odd, hence invertible modulo every power of two. Fix every seed except `s_j` and set `C` to the corresponding partial sum. The congruence becomes `2^t * u * s_j ≡ -C (mod 2^256)`. If `2^t` does not divide `C` there is no solution. If it does, there are exactly `2^t` solutions modulo `2^256`. In every case the number of bad `s_j` is at most `2^t ≤ 2^63`. Therefore the conditional probability is at most `2^63 / 2^256 = 2^{-193}`. The bound does not depend on the other seeds or on which keys were fixed, so it is worst-case over distinct 72-limb strings.

## Algorithm

All Fp12 products below are reducing dense multiplications and are charged one reference unit each. Coefficients written into tables are in `[0, p)`.

1. Decode and check the profile acceptance conditions: `Q` on the G2 curve, `Q ≠ ∞`, `[r]Q = ∞`, `z ≠ 0`, `z^r = 1`. Invalid input is rejected. The challenge distribution never feeds that branch; the work is still charged.
2. Compute `g = e(P0, Q)`.
3. Baby list, `j = 0..m-1`. Start at `γ = 1`. For each `j`, append the 576-byte key of `γ` at slot `j` and append an index entry `(F(γ), j)` with `j` a 128-bit integer. Then replace `γ` by `γ * g`, except after the last slot. One extra product produces `g^m`.
4. Compute `D = (g^m)^{r-1}`, so `D = g^{-m}`.
5. Sort the `m` index entries by `F` with LSD radix sort, digit width 16 bits, 16 passes, low digit first. Each pass is a counting sort on `2^16` counters. Each counter and each prefix is a 128-bit integer. Scatter is stable. By induction, after pass `d` the index is sorted on the low `16(d+1)` fingerprint bits, and after 16 passes it is sorted on all 256 bits. Equal fingerprints occupy one contiguous run.
6. Giant list, `i = 0..m-1`. Start at `δ = z`. For each `i`, append the key of `δ` at slot `i` and the index entry `(F(δ), i)`, then replace `δ` by `δ * D` except after the last slot. Sort these index entries by the same radix procedure.
7. Merge the two sorted index lists. A merge step compares fingerprints only. When the current fingerprints are equal, compare the two full 72-limb keys. On equality of keys with indices `(j, i)`, set `a = i*m + j`, compute `P = [a]P0`, check `[r]P = ∞` and `e(P, Q) = z`, and return `P`. The merge stops after `B` full-key comparisons if none of them was a key equality.
8. If no key equality was accepted, fail.

Slot `j` is the baby exponent and slot `i` is the giant index, so a key equality says `z * g^{-i m} = g^j`, hence `z = g^{i m + j}`. The exponents `i*m+j` for `0 ≤ i,j < m` are exactly the `m^2` consecutive integers `{0,…,m^2-1}`. They are pairwise distinct and all lie in `{0,…,r-1}` because `m^2 ≤ r`. There is therefore at most one true pair.

Addresses of key slots are 192-bit integers. The byte offset of baby slot `j` is `j * 576`. Since `j < 2^127`, this offset is `< 2^136` and fits in three 64-bit limbs (`(m * 576)` has bit length 136). An index-entry offset `j * 48` has bit length at most 133. Counter tables occupy `2^16 * 16` bytes each and sit in the same address space; a 16-bit digit shifted by 4 bits is a 64-bit offset added into a 192-bit base. Counter values never exceed `m < 2^127`, so a 128-bit counter does not overflow.

## Success

Let `S = {0,…,m^2-1}` and let `A` be the event `a ∈ S`. Uniformity of `a` gives `P(A) = m^2 / r`.

Let `X` be the number of index pairs whose fingerprints agree and whose 72-limb keys differ. For each of the at most `m^2` pairs, the lemma bounds the coin-probability of a fingerprint collision by `2^{-193}` whenever the keys differ, and a differing-key pair is the only kind that contributes to `X`. Linearity, conditional on the challenge, yields `E[X | challenge] ≤ m^2 * 2^{-193}`. Taking expectation, `E[X] ≤ m^2 / 2^{193}`, and Markov's inequality yields `P(X ≥ B) ≤ m^2 / T`.

If `A` holds and `X < B`, the true pair has equal keys, hence equal fingerprints, and it is the only such exponent pair. The merge meets at most `X` unequal fingerprint matches before or beside that true pair, then one full-key comparison that accepts it. That is at most `X+1 ≤ B` full-key comparisons, so the cap does not fire first. The returned point passes the pairing check by the reduction.

Thus

```text
P(success) ≥ P(A) - P(X ≥ B) ≥ m^2 / r - m^2 / T
           = m^2 * (T - r) / (r * T).
```

The integer comparison `100 * m^2 * (T - r) ≥ 39 * r * T` holds, so the right-hand side is at least `39/100`. The algorithm therefore succeeds with probability at least `0.39`. No heuristic distribution on GT elements is used: the collision bound is uniform over coins for every fixed pair of distinct keys, and Markov applies to the nonnegative integer `X`.

## Instruction ledger

`C = 20000`. Each dense Fp12 multiplication contributes `C` instruction-equivalents and nothing else. Word transfers below are scalar instruction-equivalents, not a hardware page-table walk: every index that does not fit in 64 bits is a three-limb add, and the per-limb charges include those adds. Upper bounds are what the total uses.

**Hash, 2304 instructions, once per stored element.** For each of the 72 limbs charge 32 instructions: load the limb, load four 64-bit seed limbs, four `mul` and three `umulh` forming the low 256 bits of a 64-by-256 product (`umulh` of the top seed limb falls at bit 256 and is omitted), three adds with carry to assemble that product, four adds with carry into a 256-bit accumulator, and the rest for the loop pointer and carry slack. `32 * 72 = 2304`.

**Key store, 576 instructions, once per stored element.** Eight instructions per limb for the 72 stores and the running 192-bit cursor (`cursor += 576` is one three-limb add). `8 * 72 = 576`.

**Index append, 128 instructions, once per stored element.** Six words (four fingerprint limbs and two index limbs) plus the 192-bit slot address.

**Radix sort, 3072 instructions per index entry, plus a fixed term.** Sixteen passes. Each pass, per entry, charge 64 instructions to extract the 16-bit digit and increment its 128-bit counter, and 128 instructions to read the prefix, advance it, form the 192-bit destination, and copy the six-word entry. `16 * (64+128) = 3072`. Separately, each pass zeros `65536` counters and scans them for prefix sums: 64 instructions to zero a counter and 32 to fold it into the 128-bit prefix. Per table that fixed term is `16 * (64+32) * 65536 = 100663296`. Two tables cost `201326592`.

**Merge fingerprint step, 64 instructions.** At most `2m` steps (two lists of length `m`), charged as `64 * 2m`.

**Full-key comparison under the cap, 2048 instructions.** At most `B` of them: 72 limbs at 16 instructions, plus address formation. `B * 2048 = 2^84 * 2^11 = 39614081257132168796771975168`.

**Pairing, at most 16384 reference units.** The ate loop count has bit length 64. Charge 64 squarings and 64 line multiplications as dense multiplications. The final exponent `(p^12-1)/r` satisfies `log2((p^12-1)/r) < 12*381 - 254 = 4318`, so binary exponentiation costs at most `2*4318 = 8636` dense multiplications. Frobenius maps and the easy part of the final exponent are charged inside the slack up to 16384. One pairing builds `g`. One pairing verifies the output.

**Subgroup and power checks.** `[r]Q`, `[r]P`, and `[a]P0` are binary scalar multiplications of at most `2*255 = 510` point operations. Each point operation is charged as one dense Fp12 multiplication, which is loose for both G1 and G2. `z^r` is at most 510 dense multiplications. The budget uses 3072 reference units on the input side and 3072 on the output side to cover these checks, canonical decoding, and the final 256-bit product `i*m+j`. Inversion `(g^m)^{r-1}` is at most 510 dense multiplications and is charged as 512.

Seed materialisation and fixed control are `2^16` instructions on top of the input-side reference units.

Per stored element the variable instruction count is

```text
side = C + 2304 + 576 + 128 + 3072 = 26080.
```

Both lists together, plus the merge steps:

```text
(26080 + 64) * (2m) = 7486945464191301001889053318461912142298368.
```

Setup instructions, input checks plus pairing plus inversion plus `2^16`:

```text
(16384 + 3072 + 512) * C + 2^16 = 399425536.
```

Verification instructions:

```text
(16384 + 3072) * C = 389120000.
```

Sum, which is the worst-case budget even on the paths that return early or abort at the cap:

```text
total = 7486945464191340615970310450630709904145664.
```

Reference-unit total is `total / C`. Its base-2 logarithm is

```text
log2(total / C) = 128.137644851442…
```

and the value rounded upward to six decimal places is `128.137645 = 128137645/1000000`. Rounding upward preserves a valid upper bound.

Preprocessing is the work that does not depend on `z`: setup, the baby list, and the baby radix fixed term.

```text
pre = 3734308784159063841211492322243855087587712
log2(pre / C) = 127.134108832415…
```

rounded upward `127.134109 = 127134109/1000000`. It is included in `total`. The giant list, merge, cap, and output check depend on `z` and are online. Nothing is amortised across challenges. Advice is empty: the seeds are coins of this run, and both tables are built from `Q`, `z`, and public profile constants.

## Memory

Peak storage, in bytes:

- two key tables: `2 * m * 576`
- two index tables and one radix scratch index: `3 * m * 48`
- one counter table: `65536 * 16`
- code, field scratch, seeds, and the live Fp12 value: `2^20`

```text
1296 * m + 2^21 = 185569945715879859593945323988805046015808.
```

`log2` of that integer is `137.091014965730…`, rounded upward `137.091015 = 137091015/1000000`. Advice bytes are 0, and `log2(max(1, 0)) = 0`.

## Accounting notes

The cap, not an expectation, bounds the full-key work. A run that sees more than `B` unequal fingerprint matches fails and still stops inside `total`. The success probability, not the instruction total, is what absorbs that tail. The two quantities are not the same claim: every execution performs at most `total` instruction-equivalents, and the fraction of challenge/coin pairs that return the correct `P` is at least `0.39`.

No finite-field NFS/TNFS cost is asserted. No Pollard walk, distinguished-point predicate, or Frobenius-orbit representative is used. The join is a sorted merge of a proved universal hash, followed by exact key equality.
