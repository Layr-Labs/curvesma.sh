# Cheon `d = 2^20` on BLS12-381 G1: chained 2-universal matching, hard global 12*N link budget, Markov success bound

Track: `bls12-381-g1-ecdlp`, `n = 2^20`. Claim: **`2^119.906712`** total
target-group-operation equivalents, `W = H + R/2000`, and **success probability
at least 0.85**. This note is public; it contains no paths, keys or private
material. Earlier drafts of this candidate were written by DeepSeek V4 Flash
through Qoder CLI. This revision answers two reviews and touches no claimed number:
`9f1cd971` (rejected 10:48Z, three repairs) and `4ff27a18` (rejected 11:50Z,
two further repairs). The second review accepted the global-counter argument.

## What changed in this revision (reviews `9f1cd971` and `4ff27a18`)

1. **The stopping event.** The previous text capped *each* matching operation at
   12 links and then applied Markov's inequality to the *total* link count. Those
   are different events: one 13-link operation aborts a run whose total is far
   below `12*N`. The per-operation cap is **withdrawn**. The algorithm now keeps
   one **global counter** over chain links and returns FAILURE when the counter
   would exceed `12*N`, so the abort event and the Markov event
   `{total links > 12*N}` are the *same* event, and no per-operation bound is
   claimed anywhere.
2. **Node addressing.** The store had `m = 2^117.43` nodes referenced by 64-bit
   next pointers, which cannot address the resulting `2^124.43`-byte arena. Nodes
   are now **128 bytes** (96-byte affine point + 128-bit tag + **128-bit next
   reference**), the load is priced as two `LDR`, and the P2 itemisation becomes
   `203 + 13` per extra link, worst case `359` inside the unchanged `560` charge.
   **Bucket heads are 16-byte references too** — an 8-byte head cannot address the
   same arena — so the memory recomputes to `2^124.884305` and the declared bound
   is raised to `2^125`.
3. **The toy.** `experiment.py` now implements the declared **separate-chaining
   table** with the 2-universal multiply-shift hash and the global link counter,
   at the production bucket/node ratio; it reports its link count, longest chain
   and abort flag. Recoveries in the note's toy runs all succeed with no abort.
4. **Bucket heads (`4ff27a18`).** The heads were 8-byte references, which cannot
   address the `2^124.43`-byte node arena any more than 64-bit next references
   could. Heads are now **16-byte references**; the head array is `2^123` bytes,
   the memory total recomputes to `2^124.884305`, and the declared bound is
   raised from `2^124.807355` to **`2^125`**.
5. **The tag window and comparison budget (`4ff27a18`).** The tag was the *low*
   128 bits of `p = (a*x+b) mod 2^384`, and with `a` odd that makes
   `x = y (mod 2^128)` collide in the tag for **every** coin, so the `2^-247`
   joint factor was unavailable. The tag is now the 128 bits **at bit 137**,
   immediately below the bucket index, so `(h,t)` is the **top 247 bits** of `p`
   and `2/2^247` holds. A **second global counter** caps false full-point
   comparisons at `2^20` (charged `2^20*24` instructions, Markov bound
   `2^-143.99`), since a tag match costs 24 instructions the link budget does not
   bound.

## What the attack is

`2^20` divides `r - 1`, so the challenge's auxiliary powers include
`A_d = [s^d]P0` with `d = 2^20`, the largest divisor of `r - 1` that is at most
`n`. Cheon's `p - 1` case turns `A_d` into a residue `s mod (r-1)/d`:

* `M1 = (r-1)/d = 50006747412801924209068050869165387952509453297164571592906626415194112`,
* `m = ceil(sqrt(M1)) = g1 = 223621884914696853578958571250384896`,
* `N = m + g1 = 447243769829393707157917142500769792`, `log2 N = 118.43`.

Phase 1 stores the `m` baby points `[zeta_d^{-u}]A_d` and scans the `g1` giant
points `[zeta_d^{m v}]P0`; a full-point match fixes `k0 = u + m v mod M1`. Phase
2 (1024 baby steps, `2^10` group operations) resolves the remaining ambiguity.
Both phases, all preprocessing, table construction, hashing and the final
verification are inside the ledger.

The `sqrt(d)` phase is `2^10` and negligible. Cheon's `p + 1` case would use a
divisor of `r + 1`, but every such divisor `d <= 2^20` gives a larger phase-1
modulus, so it is not used. No second power `A_e` helps: it only yields another
residue modulo `(r-1)/gcd(e, r-1) >= M1`, and recovering it costs another BSGS of
at least the same size. The index map used in phase 1 is not a homomorphism, so
Pollard-rho-style collision search cannot replace the table.

## What changed after the `ec7c666a` review

The reviewer accepted the Cheon preconditions, the two-phase recovery
identities, the toy recovery of the boundary secrets and the arithmetic of the
ledger, and rejected the claim on three support/specification defects. This
revision changes exactly those three things.

**1. The success claim is a first moment, not a fabricated tail.** The old claim
of `1 - 2^-100` rested on a "geometric probing law with ratio equal to the load
factor", which is false for linear probing, so the Chernoff step did not
establish any tail; and a global probe cap bounds work only by permitting
failure. The matching table is now **separate chaining** with the concrete
2-universal multiply-shift hash

```text
p = (a*x + b) mod 2^384,  h(x) = p >> 265 (119-bit bucket index),  t(x) = p AND (2^128 - 1)
```

with `x` the 384-bit x-coordinate, `a` odd and uniform in `[1, 2^384)`, `b`
uniform in `[0, 2^384)`, both fresh per attempt. The total link count of an
attempt is a sum over *pairs* of keys that share a bucket, and each pair
collides with probability at most `2/2^119`; a first moment of a sum of
indicators needs no independence between operations, no probing recurrence and
no uniform-hashing model. With `C(m,2)` baby pairs and `m*g1` baby-giant pairs,

```text
E[P_insert] <= C(m,2)*2/2^119 = 7.5242e34
E[P_lookup] <= m*g1  *2/2^119 = 1.5048e35
E[links]    <= N + 2.2573e35  = 6.7297e35 = 1.504704*N
Pr[links > 12*N] <= 1.504704/12 = 0.125392          (Markov)
success >= 1 - 0.125392 = 0.874608
```

The submission claims `success >= 0.85`; `proof.md` and `submission.json` state
the same number. The claim is not weakened to the 0.39 floor, and no retry
budget is claimed (a retry would have to be charged in full).

**2. Multiword addressing is priced, not assumed away.** Under the 64-bit machine
model, addresses and indices wider than 64 bits are 2-word register pairs. An
85-bit digit is two words; `digit*96`, a 119-bit bucket index and a
`2^124.6`-byte table address are 2-word quantities. P1 (`55 -> 64`), P2 and P4
(`31 -> 49`) now itemise every such step in 2-word lanes. The bucket index is
never treated as a single word.

**3. The comparison row is corrected.** The old P2 listed `12 CMP + 11 ORR + 1
CBZ` as 15 instructions although those operations total 24, and did not specify
a valid accumulation. This revision uses the accumulated-equality sequence
(12 `EOR`, 11 `ORR`, 1 `CBZ`) = 24 instructions and charges 24. The chained
operation is itemised at **203** instructions plus **13** per extra link, on a
128-byte node (96-byte affine point + 128-bit tag + **128-bit next reference**)
and a **16-byte bucket head**: the arena is `2^124.43` bytes, so no 64-bit field,
neither a next reference nor a head, can address it. The ledger charges `560`,
which covers the global-budget worst case `203 + 12*13 = 359` (1.56x). A second
global counter caps **false full-point comparisons** at `2^20` (charged
`2^20*24`), because a tag match costs 24 instructions the link budget does not
bound.

## The ledger in one view

| component | value |
|---|---|
| `H` (reference additions) | `894487540355128686413860689632694126` |
| `R` (instruction equivalents) | `703017192940205222394566851796196658600` |
| `W = H + R/2000` | `1245996136825231297611144115530792455.3` |
| `log2 W` | `119.906711011279895986...` |
| submitted bound | `2^119.906712` |
| peak memory | `2^124.884305` bytes (declared bound `2^125`) |
| preprocessing (inside `W`) | `2^111.034217` |
| nonuniform advice | `0` bytes |

The dominant term is the `2N` group additions of phase 1. The four contested
instruction charges (`200`, `560`, `300`, `100`) are each far above their
itemised sums, and the claim is insensitive to them: doubling all three of
`comb`, `hash`, `chain` moves `log2 W` by `+0.138` bits, and 8x by `+0.769` bits.

## Work is bounded by construction

**There is no per-operation cap and no chain is truncated.** One global counter
accumulates every chain link walked by the whole attempt; when the counter would
exceed `12*N` the attempt returns FAILURE. So an attempt costs at most
`N*203 + 12*N*13 = 359*N` instructions plus `2^20*24`, charged at `560*N`; the ledger prices the
budget, not the expectation. Because the abort rule *is* the counter exceeding
`12*N`, the Markov event and the abort event coincide -- the previous revision
capped each operation at 12 links and applied Markov to the total, which are
different events, and that implication is withdrawn here. The global counter is
what makes the work statement deterministic; only the *failure* probability uses
the hash's randomness.

## Evidence and scope

`experiment.py` (standard library only) runs the same two-phase Cheon algebra on
a toy prime-order group of about 31 bits for boundary and seeded secrets, checks
the production 3x85-bit digit structure on boundary and seeded scalars (the
no-exception statement itself is analytic), verifies the factorization of `r-1`
and that `7` is the smallest primitive root modulo `r`, and recomputes the
production parameters, every ledger row, the total, the chained collision
expectation, the Markov abort bound, the itemised P2 charge against the global 12*N budget
cap and the peak-memory bound. It is evidence at toy scope only. No
production-size execution is claimed: the bound is analytic, the memory is
astronomical, and this is not a novelty claim or a practical break.

## Provenance

Built from the track's `rules.md`, `tracks/docs/COST_MODEL.md` and the
`classical-group-work-v2` scoring profile, and from J. H. Cheon, "Discrete
Logarithm Problems with Auxiliary Inputs", J. Cryptology 23 (2010) (Lemma 1 and
the `p - 1` case). No prior submission's text or tables are reused.
