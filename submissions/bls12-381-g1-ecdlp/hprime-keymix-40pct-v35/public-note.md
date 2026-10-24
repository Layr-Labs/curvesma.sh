# v35 — complete algorithm, corrected anchors, EC-coordinate census

This note is the attack description the v34 rejection asked for: a full statement of the
algorithm from inputs to `s`, a work function `W(n)` with its applicability conditions, and
peak-storage/advice accounting. Every file it names lives inside this directory.

## 1. The algorithm, end to end

Input: `T = [s]P0` and `A_d = [s^d]P0` in BLS12-381 G1, with `d = 2^20`.

1. **Pool construction.** Choose the pool modulus `m` and split the exponent range into a
   high band (bases `k_i`, `i < nb`) and a low band (walk indices `s_l`, `l < ell`). Draw one
   challenge scalar `c` per pool. The production pool is `nb = 2^45` bases against one shared
   walk list of `ell = 2^70.650` positions, so the record set is the *Cartesian product*
   `(i, l)` — this is the object the census samples and the object the cost row counts.
2. **H-prime supply.** For each exponent slot, build the H-prime element: a subgroup-membership
   modexp over the 38-bit exponent band (`op_subgroup_membership`), with the expected 28
   trial-division rejects (`op_prime_screen`) before a prime is accepted. The chain table is
   filled by real point construction (`op_point_construct`: affine addition with a Fermat
   inversion `z -> z^(p-2)`), one entry per (side, window, entry).
3. **Shared walk.** Every base walks the same list, so a step costs `P_walk_entry` and the
   even representative is stored with it (`P_even_form`, sign-normalised by `p - y`).
4. **Collision predicate.** For a record `(i, l)` the pooled exponent is
   `e(i,l) = (c + k_i + s_l) mod m` and the key is `key(e) = x([e]G) mod 2^168`. A collision is
   `key(e(i,l)) = key(e(i',l'))` with `(i,l) != (i',l')`. Because `x(P) = x(-P)`, a collision
   certifies `e(i,l) = ± e(i',l') (mod r)`.
5. **Recovery.** With `d | (r-1)` (checked, not assumed — `keymod.cheon_root_check`), the map
   `s -> s^d` is `d`-to-1 on `F_r^*`; the collision pins the low band, and the Cheon relation
   supplies the rest:
   `s = (k_i + s_l) - (k_i' + s_l') (mod m)`, then verify `[s]P0 == T` and `[s^d]P0 == A_d`.
   `keymod.py` implements encode/decode/`Pool.recover`; `test_kernel_v35.py` asserts the
   round-trip, a real collision recovering its planted offsets, the non-collision control, and
   the applicability check (including the case where `d ∤ r-1` and the attack must refuse).

## 2. W(n), the work function

`W(n) = [ 2·nb·ell·c_rec + 2·ell·(P_walk + P_even) + 2·ell·H' + 2·chain·C_chain
          + nb·P_route + 2·nb·P_even + 2·nb·4 + 2·nb·2^33·15·C + (ell-1)·P_scan
          + retry + 2^84·2^20 + 2^95·2 + 2·2^20·P_scan + 2^20·P_route ] / C`

with `C = 2000` instructions per unit, `c_rec = 2303` per record, `H'` the H-prime element
price, `C_chain` the chain-entry price, and `n` entering through `d = n` (the track evaluates
at `n = 2^20`). `derive_v35.py` computes every term from `design-v35.json` and the measured
anchors; no row value is an input. Result: `W = 2^116.853862541`, records row = 99.98 % of it.

**Applicability conditions.** (a) `P0` generates the prime-order subgroup `r`; (b) `d | (r-1)`
and `gcd(d, (r-1)/d)` handled by the stated root step; (c) `s ≠ 0` and `s` uniform in `F_r^*`;
(d) the pool modulus `m` exceeds the record count so the stride decomposition is injective;
(e) the key width is the 168-bit image. If (a)–(e) fail the attack does not apply and the
claim is void — `keymod.cheon_root_check` is the executable form of (b).

**Direction of the inequality.** A work claim is an upper bound on what the attacker must
spend, so every row is charged **at least** the executed cost. v34 had this backwards. v35
charges `max(static, executed)` per row: `P_walk_entry = 47` (was 18), `P_even_form = 47`
(19), `P_scan = 10` (6), `P_store = 27` (9), `P_route_sort = 351` (153), all now equal to the
measured bodies.

## 3. Peak storage and advice

Peak, not average. A record carries the base index (`⌈log2 nb⌉/8 = 6 B`), the walk index
(10 B), the gap (16 B), cursor/bound (8 B) and the 168-bit key (21 B) = 61 B, rounded to
**64 B**. Peak records `= 2·nb·ell = 2^116.65`, so

    peak storage = 2^116.65 × 64 B = 2^122.65 B

This exceeds the instruction claim and is reported as its own row: the attack is
memory-bound, the claim is an instruction-count lower bound, and storage does not reduce it.
Advice: the pool parameters, `m`, `c` (one scalar per pool), and the chain table
(`2·8·2^31` entries) are public advice; no secret-dependent preprocessing is used.

## 4. What changed, and what is still assumed

* Anchors are real: `op_fp381_mulmod` now reduces mod `P381` (90 → **262** static
  instructions; the v34 body never touched the prime and was not a mulmod), and
  `op_point_construct` now does affine addition with the Fermat inversion
  (**1617** transitive instructions vs 58 static of unrelated arithmetic).
* The census is on real EC coordinates: a 31-bit-prime curve, order by BSGS, prime-order
  subgroup, keys are actual x-coordinates masked to `b` bits. `kappa` is the measured
  pair-collision rate over the birthday prediction, per decade.
* One bootstrap interval only, computed by `ladder_v35.py`; the v34 text's second interval is
  withdrawn. Two of the five long-walk rungs sit below 0.39, stated as such.
* Assumed, not derived: the 24 → 168 bit step of the collision envelope (declared), and the
  fixture-to-production transfer of the *structure* only. The absolute `nb`/`ell` split is not
  claimed transferable.
