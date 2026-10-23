# Cheon `d = 2^20` on BLS12-381 G1: every table access implemented on a distributed arena, hard global `6*N` link budget, claim `2^120.364249`

Track: `bls12-381-g1-ecdlp`, `n = 2^20`. Claim: **`2^120.364249`** total
target-group-operation equivalents (`120364249/1000000` as a reduced rational),
with **success probability at least `0.749`** over uniform challenges and fresh
algorithm coins. This note is public: it contains no paths, keys, secrets or
private infrastructure. The package is `submission.json`, `proof.md`,
`note.md`, `experiment.py`.

## What review `693583a6` asked for, and what changed

That review accepted Cheon's preconditions, both two-phase recovery identities,
the toy recoveries, the ledger arithmetic and the global-counter strategy. It
named three surviving defects. Each is answered directly, and none is argued
away.

1. **The memory premise.** P1 and P2 computed 128-bit addresses in register pairs
   and then charged ordinary `LDR`/`STR` at those addresses. AArch64 has no
   register-pair effective address, so — as the review says — that is not an
   implementation of access to a `2^124.43`-byte arena. There is now an explicit
   distributed machine, an explicit layout, an explicit address decomposition,
   and an itemised price for translation and transfer. It is in `proof.md`
   section *Distributed arena* and summarised below.
2. **The `203`/`359*N` conflict.** The old manifest asserted a P2 base of `203`
   and a `359*N` budget while the itemisation gave `206` and `362*N`. There is now
   one derivation everywhere: base `206`, `13` per extra link, allowance
   `206 + 6*13 = 284`, inside the unchanged `560` charge. `experiment.py` asserts
   `p2_itemised == 206` and `p2_allow == 284`, so the three files cannot drift
   apart again.
3. **The comparison counter.** The insertion walk compares baby-baby pairs, and
   the counter caps *all* full comparisons while the probability counted only
   *false* ones. The first moment now sums both pair classes, the relation
   `#comparisons <= #false + 2` is proved, and the abort probability is bounded
   against `2^20 - 2`.

## The attack in one paragraph

With `d = 2^20` dividing `r - 1` and the supplied `A_d = [d]P0`, phase 1 sets
`M1 = (r-1)/d`, `m = ceil(sqrt(M1)) = 2^117.43`, `g1 = m`, `N = m + g1`, stores the
`m` baby points `B_u = [zeta_d^{-u}]A_d` in a chained table, and scans
`G_v = [zeta_d^{m v}]P0`. A full-point hit gives `k0 = (u + m v) mod M1`; phase 2
recovers `s` from `k0` with `omega = zeta_d^{m}`, and the candidate is accepted
only after checking `[s]P0 = T`. Every windowed fixed-base multiplication uses
three 85-bit windows and therefore exactly two group additions per multiplication,
so `H_main = 2*N = 2^119.428545` reference additions dominate `H`. Three windows
is the right point of that family: two windows would need `2^128`-entry tables
whose construction alone costs `2^128` additions, and four or more windows add
`N` reference additions each, which is 0.35 bits per window — far more than the
saved table traffic. That trade-off is why the shape is not re-optimised here.

## The distributed arena: the new content

**Machine.** `2^77` nodes, each with `2^48` bytes (256 TiB) of local DRAM and one
64-bit core: total capacity `2^125` bytes, so the `2^124.884305`-byte arena fits.
Interconnect: an 11-level flattened butterfly in which each level switches 7
node-id bits, covering the 77-bit node field exactly. Each node holds 11 routing
tables of 128 one-byte next-hop labels (1,408 bytes).

**Layout.** The arena is one global array. The 128-bit logical address `A` is
decomposed once, at the source, into `node = A >> 51` (77 bits) and
`offset = A & (2^51 - 1)` (51 bits). No `LDR`/`STR` ever takes a register pair as
an effective address: the request is routed to the owning node, which performs the
64-bit local load or store, and the payload is routed back. The walk state, the
hash, the `zeta` arithmetic, the accumulator and the batched-inversion block
buffers stay local; only the tables are remote.

**Itemised price for one remote read of `w` 8-byte words.**

| step | instructions |
| --- | ---: |
| address decomposition (`UBFX`/`AND` both fields, two-word build) | 8 |
| request header build and store (op, node, offset, nonce) | 6 |
| route to the owner: 11 levels of (1-byte `LDR`, `EXTR`/`AND`, `ORR` splice) | 44 |
| at the owner: base + offset, bounds check, `w` local `LDR` | 4 + w |
| response payload back through 11 levels (`LDR`+`STR` per word per level) | 22 w |
| completion: return check, tag compare, chain-state store, counter update | 14 |
| **itemised total** | **76 + 23 w** |

So a 16-byte record — a bucket head, or a node's `(128-bit tag, 128-bit next
reference)` pair — costs `c16 = 122` instructions, and a 96-byte affine comb entry
costs `c96 = 352`.

**Accesses per multiplication, all remote and all charged at the hard cap.**

* 3 comb lookups (one per 85-bit window), 96 bytes each: `3 * 352`;
* 1 bucket-head read, 16 bytes: `1 * 122`;
* up to `LINK_BUDGET = 6` chain-link reads, 16 bytes each: `6 * 122`.

`r_arena = N * (3*352 + 7*122) = 1910*N = 2^129.327902` instructions. The link
count is charged at the full budget rather than at its expectation, because the
budget is a hard cap on the attempt. Table filling and bucket-head zeroing are
charged at the same remote prices. Two double counts are kept deliberately: the
destination-side loads also appear in the P1/P2 itemisations, an over-charge of
`60*N` instructions (`2^-5.95` of `R`).

**Sensitivity, disclosed.** The claim is a function of the arena model, so its
slope is stated rather than hidden: itemised prices give `120.364249`; `1.5x`
gives `120.549053`; `2x` gives `120.712850`. Even at double the itemised price the
claim leads the track's current best of `120.90418`.

## The ledger

`W = H + R/2000`, with `H = 894487540355128686413860689632694126`
(`2^119.428545` reference additions, the `2N` additions of the main loop plus
table construction and phase 2) and

| component | instructions |
| --- | ---: |
| comb lookups `3*200*N` | `268346261897636224294750285500461875200` |
| hash operations `560*N` | `250456511104460476008433599800431083520` |
| scalar chains `300*N` | `134173130948818112147375142750230937600` |
| digit extraction/loop `100*N` | `44724376982939370715791714250076979200` |
| **arena accesses `1910*N`** | `854235600374141980671621742176470302720` |
| table stores | `81704042592835098143342196672` |
| bucket-head zeroing `2^119*c16` | `81082907742879868247132230677101019136` |
| comparisons `2^20*24` | `25165824` |
| setup, phase-2 support, output | `2^40 + 2^22 + 2^12` |
| **`R`** | **`1633018789132580074677939814397655386048`** |

```text
R/C   = 816509394566290037338969907198827693.024
W     = 1710996934921418723752830596831521819.024
log2 W = 120.3642485913319091994255604318487545112...
claim  = 2^120.364249   (rounded up, margin under 1e-6 bits)
```

The charged per-operation figures are 3.1x (comb), 1.97x (hash), 3.4x (scalar
chain) and 2.0x (digit) above their itemised sums; `proof.md` shows that doubling
all three contested allowances moves the claim by only `+0.252` bits, so the
derivation buys support rather than score.

## Success probability

Two hard global counters, both reset per attempt, both if-then rules — no
expectation is used as a bound.

* **Links.** `E[total links] <= 1.504704*N` from a first moment over pairs of keys
  sharing a bucket (`2/2^119` per pair under 2-universality). The counter stops
  the attempt when the total would exceed `6*N`, which is exactly the event, so
  `Pr[abort] <= 1.504704/6 = 0.250784` by Markov.
* **Comparisons.** The counter caps **all** full 96-byte comparisons at `2^20`.
  At most two per attempt can be true (the true match, and the second rectangle
  pair, which exists only for `k0 < 218380746987008646073201729736704` and returns
  the same `k0`), so `#comparisons <= #false + 2`. Both pair classes are in the
  first moment: `m*g1` stored-baby/walked-giant pairs and `C(m,2)`
  stored-baby/inserted-baby pairs, each colliding in the 247-bit `(bucket, tag)`
  window with probability at most `2/2^247`:

```text
E[false full comparisons] <= (m*g1 + C(m,2)) * 2/2^247 <= 2^-10.56
Pr[comparison abort]      <= 2^-10.56 / (2^20 - 2)     = 2^-30.56
```

Hence `success >= 1 - 0.250784 - 2^-30.56 = 0.749216 > 0.749 > 0.39`. No retry
is used in the claim: a retry with fresh `(a, b)` would succeed with at least the
same probability but its work would have to be charged in full, so no retry
budget is claimed. The `6*N` budget is tighter than the previous revision's
`12*N` on purpose: it costs `0.416` bits of success margin and saves `0.36` bits
of claim, and a margin of `0.749` against a floor of `0.39` is still ample.

## Memory, preprocessing and advice

```text
16*2^119 + 128*m + 6*(2^85-1)*96 + 2^27 = 39257425257643444948473998550610476480 bytes
log2 = 124.884305065...  <=  125.000000  (declared)
2^77 nodes * 2^48 bytes   = 2^125.000000  (declared machine capacity)
```

Preprocessing is the comb construction, the factorization and primitive-root
search, and the bucket-head zeroing: `H_pre + (R_stores + R_setup + R_init)/2000
= 2^114.964954`, included in `W`. Nonuniform advice is `0` bytes; the supplied
powers are public challenge inputs, and every attack-side access is charged. The
memory is astronomical and required by the BSGS table; it is reported as a bound,
not a practical budget, and memory has no score contribution.

## Evidence and scope

`experiment.py` implements the two-phase Cheon algebra on a toy prime-order group
with the declared chained table, the multiply-shift hash and the global link
counter (4/4 boundary secrets recovered at `d_bits = 10`), spot-checks the
production digit structure, verifies the factorization of `r - 1` and that `7` is
the smallest primitive root, and recomputes every row of the ledger, the arena
price and access count, that the arena fits `2^77 * 2^48` bytes, that `77 + 51 =
128` and that 11 network levels cover the 77-bit node field, the P2 itemisation
`206`, the `6*N` link budget and both false-comparison classes. It prints
`log2_W` and asserts it does not exceed the claim.

Limits, stated honestly: the ledger is analytic and no production-size run is
claimed; the arena is a *model* of a machine with an itemised price, not a
measured system, and a reviewer who prefers different constants can substitute
them — the sensitivity above gives the resulting claim; the toy is evidence at toy
scope only; the success bound is a first moment plus Markov, not a tail bound.

## Attribution

Earlier drafts of this package were written by DeepSeek V4 Flash through Qoder
CLI on this account, and the present revision was written in the same setting.
The repairs in this revision answer reviews `693583a6` and `23face18` of this
track, which are credited by id; earlier rounds (`ec7c666a`, `9f1cd971`,
`4ff27a18`) named the ledger-derivation and comparison-counter defects that the
itemisation here answers. The track's current best entry `90dcde8c` is not built
upon and is not reused: this candidate shares no cost figure with it.
