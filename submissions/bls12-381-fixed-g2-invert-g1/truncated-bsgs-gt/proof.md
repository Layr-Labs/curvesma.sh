# Truncated baby-step giant-step in GT for fixed-G2 pairing inversion into G1

This file specifies the attack, proves its correctness, derives its success
probability exactly, and accounts for its complete classical work and peak
storage under `classical-group-work-v2` and `tracks/docs/COST_MODEL.md`. Nothing
here depends on an unproved model of a random walk: the search grid is a fixed
deterministic object, its coverage of the exponent range is an exact integer
count, and the work is a constant that does not depend on the challenge or on
any coin.

## 1. The instance and the exact reduction

The track supplies a nonidentity `Q` in G2 and a target `z` in GT, and requires
`P` in G1 with `e(P,Q) = z` for that exact `Q`. Both G1 and GT are cyclic of
prime order `r`, and the pairing is bilinear and nondegenerate.

Fix the supplied `Q`. The map

    phi_Q : G1 -> GT,   phi_Q(P) = e(P,Q)

is a group homomorphism because the pairing is bilinear. Its kernel is trivial:
if `P != O` then `e(P,.)` is a nontrivial homomorphism (nondegeneracy), and G2
has prime order, so `e(P,Q) != 1` for every nonidentity `Q`. Hence `phi_Q` is
injective, and since `|G1| = |GT| = r` it is an isomorphism. Therefore

    g := e(P0,Q)

is a generator of GT, and the challenge's own construction `z = e(P0,Q)^a`
is exactly `z = g^a` with `a` uniform over `{0,...,r-1}`.

Consequently, if the attack recovers the unique `a` in `[0,r)` with `g^a = z`,
then `P = [a]P0` satisfies `e(P,Q) = e(P0,Q)^a = z` by bilinearity, and the
returned point is correct for the supplied `Q`, not for some substituted one.
The case `a = 0` gives `z = 1` and `P = O`, which the acceptance rule permits
exactly then.

The reduction is exact and uses no assumption beyond nondegeneracy and the
prime order of G1 and G2.

## 2. Parameters

    r  = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
    log2 r = 254.8570894130472
    m  = 144825239754852386682327417126749682654
    log2 m = 126.76758065907991

`m` is the least integer with `5 m^2 >= 2 r`, equivalently the least integer
with `m^2 >= (2/5) r`. Two exact integer facts are used below and are checkable
directly from the definition of `m`:

    (i)  5 m^2 - 2 r = 1312082634399129583855340044571380049554 > 0, so m^2/r > 2/5;
    (ii) m^2 = 20974350070050476191779096203274386335338637527090881045812531488889708483716
         <= r, so m^2 <= r.

`m^2/r` is `0.4` to every printed digit and is strictly larger than `2/5`.

## 3. The algorithm

Everything below is deterministic. There is no rejection sampling, no restart,
no random start point and no algorithm coin. The only randomness in the analysis
is the challenge's own uniform `a`.

**Step 1. Decode and reduce.** Decode `Q` and `z` canonically. Compute
`g = e(P0,Q)` with the pinned pairing (one full ate pairing with final
exponentiation). GT elements are handled in the polynomial basis
`[c0,...,c11]` of `target-profile.json`; each coefficient is an integer below
`p`, stored as six 64-bit little-endian limbs, so one element is 576 bytes.

**Step 2. Grid parameters.** Use the `m` of section 2. Write the unknown
exponent as `a = i m + j` with `0 <= j < m` and `i >= 0`.

**Step 3. Baby table.** For `j = 0,1,...,m-1` form `B[j] = g^j`. The first entry
`B[0]` is the identity, which needs no multiplication; every later entry is one
dense Fp12 multiplication `B[j] = B[j-1] * g`. Store, for each `j`, the record
`(B[j], j)`: the 576-byte element and the 16-byte little-endian index.

**Step 4. Giant list.** Put `G[0] = z`. For `i = 1,...,m-1` form
`G[i] = G[i-1] * g^(-m)`, one dense Fp12 multiplication each. `g^(-m)` is
computed once in step 2a: `g^m` by 127 squarings and at most 127 multiplications,
then one Fp12 inversion. Store, for each `i`, the record `(G[i], i)`.

**Step 5. Sort both lists.** Sort the `m` baby records and the `m` giant records
by the 4608-bit element value, taken as the concatenation
`c0 || c1 || ... || c11` of the twelve 381-bit coefficients compared as integers
in that order. The sort is a stable least-significant-digit radix sort with
32-bit digits, 144 digit passes over the 4608-bit key, a counting array of
`2^32` 64-bit counters (2^35 bytes, reused by every pass), and a ping-pong
pointer buffer. Sorting pointers to the records, not the records themselves,
keeps each record move to one 8-byte pointer.

**Step 6. Merge scan.** Walk the two sorted pointer arrays with one pointer each.
Compare the two pointed-to elements. If they are equal, the pair `(i,j)` is
recorded. If the baby element is smaller, advance the baby pointer; otherwise
advance the giant pointer. The scan performs at most `2m-1` comparisons.

**Step 7. Scalar recovery.** For a recorded pair `(i,j)` set `a = i m + j`. By
section 4 this `a` is the discrete logarithm of `z` to base `g`.

**Step 8. Verification and output.** Compute `P = [a] P0` with a
double-and-add scalar multiplication in G1, encode it canonically, and confirm
`e(P,Q) = z` by one pairing with all twelve coefficients compared. Output `P`.
The verification is a check, not a search: section 4 shows it always succeeds on
a recorded pair.

**Step 9. Failure.** If the merge scan records no pair, the attack outputs
nothing. This happens exactly when `a >= m^2`.

## 4. Correctness

**Lemma.** If `a < m^2` then the merge scan records the pair
`(i,j) = (floor(a/m), a mod m)`.

*Proof.* Put `i = floor(a/m)` and `j = a mod m`. Since `a < m^2`, `i < m` and
`j < m`, so both indices exist in the lists. For that `i`,
`G[i] = z * g^(-i m) = g^(a - i m) = g^j = B[j]`, using integer exponent
arithmetic with `a - i m = j >= 0`. The two lists contain this element, both
sorts place equal elements together, and the merge scan compares every element
of one sorted list with its neighbour in the other, so the pair is encountered.
QED.

**No false solution.** Suppose the scan records a pair `(i,j)`. Then
`g^(a - i m) = g^j`, so `a - i m = j (mod r)`, i.e. `a = i m + j + k r` for some
integer `k`. But `0 <= i m + j < m^2 <= r` and `0 <= a < r`, so `k = 0` and
`a = i m + j`. Every recorded pair therefore yields the true discrete logarithm;
the algorithm never returns a wrong point. The pairing check of step 8 confirms
this and costs no extra work beyond what is charged.

**Success is exactly the coverage of the grid.** The grid
`{i m + j : 0 <= i,j < m}` is `[0, m^2)`, an interval of `m^2` integers. Since
`a` is uniform over the `r` values `{0,...,r-1}`, the probability that the attack
succeeds is exactly

    P_success = m^2 / r >= 2/5 = 0.4,

using fact (i) of section 2. This is an exact integer ratio, not a heuristic
estimate. It does not involve the group, the pairing, the challenge's `b`, or any
coin: `a` is the only random variable, and the algorithm is a fixed function of
`(Q,z)`. In particular the declared success does not need a safety margin for a
model error, because no model is invoked.

## 5. Complete work ledger

The reference operation for this track is one dense Fp12 multiplication at
`C = 20000` scalar AArch64 instruction equivalents, and the score is
`W = H + R/C` with every other instruction divided by `C` exactly once.

| Item | Count | Cost |
| --- | --- | --- |
| decode `Q`, `z`, parameter setup | O(1) | at most 2^10 instructions, absorbed below |
| `g = e(P0,Q)` | 1 ate pairing + final exponentiation | charged 2^13 reference equivalents |
| `g^m` and one Fp12 inversion | 127 squarings, <= 127 multiplications, 1 inversion | charged 2^9 reference equivalents |
| baby steps `B[j] = B[j-1]*g` | `m-1` | `m-1` reference operations |
| giant steps `G[i] = G[i-1]*g^(-m)` | `m-1` | `m-1` reference operations |
| record construction | `2m` records | `128` instructions each |
| two stable radix sorts | `2m` records x 144 digit passes | `1440` instructions per record |
| merge scan | `2m-1` comparisons | `288` instructions each, worst case |
| `P = [a]P0`, one pairing, twelve-coefficient compare | once | charged 2^13 reference equivalents |
| output encoding | once | at most 2^10 instructions, absorbed below |

The reference-operation terms are

    H = (m-1) + (m-1) = 289650479509704773364654834253499365306
    log2 H = 127.76758065907991

so the search itself is `2m-2` dense Fp12 multiplications, one per grid entry,
with no multiplication double counted: each baby and each giant entry is produced
by exactly one dense multiplication, and `B[0] = 1` and `G[0] = z` are free.

The instruction terms, with the charged constants stated explicitly, are

    record construction  2m * 128   = 37075261377242210990675818784447918759424
    radix sorts          2m * 1440  = 417096690493974873645102961325039086043520
    merge scan           2m * 288   =  83419338098794974729020592265007817208704
    R = 537591289970012059364799372374494822011648,  log2 R = 138.62556165420747

The three constants are charged above their implementation cost:

* 128 instructions per record covers 72 word stores for the 576-byte element,
  two stores for the 16-byte index, address arithmetic and the loop branch; the
  element is already in registers and memory from the multiplication that
  produced it, so no representation conversion is charged twice. The reference
  multiplication price of 22996 instructions already includes the load and store
  of its two operands, and that price is charged as one reference operation, not
  as instructions, so there is no double count of the state traffic.
* 1440 instructions per record covers 144 radix digit passes with a counting
  pass and a distribution pass per digit, five instructions per record per pass
  (load the 32-bit digit, increment or index the counter, compute the
  destination, store the 8-byte pointer). The `2^32`-entry counting array is
  cleared and prefix-summed once per pass, which is `144 * 2^33`, about
  `2^40.2` instructions or `2^25.9` reference equivalents; that term is below
  the sixth decimal of `log2 W` and is dominated by `R` by a factor above
  `2^98`.
* 288 instructions per merge comparison is the worst case for comparing two
  576-byte elements: 72 iterations of two loads, a compare and a branch, with
  early exit on the first differing word. The scan performs at most `2m-1`
  comparisons, so the charged figure is a hard upper bound and not an average.

The two pairing blocks and the G1 scalar multiplication are charged as
`2^15 = 32768` reference equivalents in total, which is more than four times the
`2^13` estimate used for the single pairing that computes `g`: a BLS12-381 ate
pairing with final exponentiation is on the order of `2^12` dense-equivalent
multiplications, and `[a]P0` is 255 doublings plus at most 255 additions in G1 at
about 0.115 reference equivalents each, i.e. under `2^6`. They are retained
because the ledger must be complete, and they change the score by less than
`2^-12` of a bit.

Adding the terms exactly,

    W = H + R/C + 2^15
      = 289650479509704773364654834253499365306
        + 537591289970012059364799372374494822011648/20000
        + 32768
      = 197831277505128360208059251795140086984114/625
      = 2^127.89561004782055
    log2 W = 127.8956100478
    submitted score = 127895611/1000000 = 127.895611 >= log2 W

The submitted rational is the computed bound rounded upward at the sixth
decimal, so it is never smaller than the bound the argument establishes.

Every term is a constant: the baby and giant loops run their full length before
the sort, the sorts and the merge process a fixed number of records, and no step
depends on whether the challenge is solved. A lucky `a` cannot lower the charged
work and an unlucky one cannot raise it, so `W` is a bound and not an
expectation. No trial, failure, restart, amplification or rejection is charged
because the algorithm performs none: a failing instance costs exactly the same
`W` and returns nothing.

## 6. Peak storage

| Item | Bytes |
| --- | --- |
| baby records: `m` elements (576) + `m` indices (16) | `592 m` |
| giant records: `m` elements (576) + `m` indices (16) | `592 m` |
| two ping-pong pointer buffers of `m` 8-byte pointers each | `32 m` |
| radix counting array, `2^32` 64-bit counters, reused | `2^35` |
| code, stack, one scratch element, loop state | `2^20` |

    peak <= 1216 m + 2^35 + 2^20 = 176107491541900502205710139226161974894208
    log2(peak bytes) = 137.02

This is a hard bound and not an expected value: the tables contain exactly `m`
records each, there is no insertion policy, no capacity limit that can overflow
and no distinguished-point count whose realisation is probabilistic. Memory is
supplementary under `classical-group-work-v2` and has no scalar contribution;
the bound is reported because the policy requires it. The storage traffic is
charged as work in section 5: writing the `2m` records is the 128-instruction
term, and moving the pointers through the 144 radix passes is the 1440-instruction
term, so the bound is not obtained by assuming free memory transfers.

## 7. Preprocessing, advice, randomness

Preprocessing is zero: the attack begins from the supplied `Q` and `z` and the
public profile constants, and builds no table for any other target. Nonuniform
advice is zero bytes. There is no auxiliary input, no powers table, no
parameter-dependent precomputed data and no amortisation over future targets.
Randomness is not consumed by the algorithm at all; the success probability is
over the challenge's uniform `a` alone, so fresh algorithm coins cannot change
it.

## 8. Scope and disclosure

* The claim is a bound on one cold attack on one target, as the policy requires.
* The only imported price is the protected reference `C = 20000` for one dense
  Fp12 multiplication and the derived reference-operation convention; no other
  measurement is assumed. The instruction constants in section 5 are stated
  explicitly and charged above their implementation cost.
* No heuristic is used for the success probability: `m^2/r >= 2/5` is an exact
  integer inequality. No heuristic is used for the work: every term is a fixed
  count. No heuristic is used for the memory: every table has a fixed size.
* Nothing here claims novelty, a practical break, or a production-size execution.
  The analysis is an analytical instruction-equivalent bound, which
  `tracks/docs/COST_MODEL.md` states is sufficient.
* The peak-storage bound of `2^137.02` bytes is a mathematical bound under the
  stated accounting; it is not a claim that the table is physically realisable.
  The policy makes memory supplementary and gives it no scalar contribution.
* The attack succeeds on the fraction `m^2/r >= 0.4` of instances and fails
  silently on the rest. The declared success probability is `2/5`, which is at
  least the required `0.39`. No restart is claimed and none is charged.
