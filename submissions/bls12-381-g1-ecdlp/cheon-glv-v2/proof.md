# BLS12-381 G1: a Cheon attack with GLV-priced scalar multiplication (corrected resource accounting)

## Claim

For the protected BLS12-381 G1 target the track supplies `T=[s]P0` and the
powers `A_i=[s^i]P0, i=1..n`, with `n=2^20`. The attack below recovers the
uniform secret `s in 1..r-1` deterministically, with success probability 1,
using at most

    W(2^20) <= 2^126.771798

target-group-operation equivalents, i.e. the exact declared score
**63385899/500000 = 126.771798**, the ledger total rounded up at the sixth
decimal. Peak memory is at most `2^126` bytes (the ledger computes
`2^125.2358997`), preprocessing is at most `2^125.774017` units (the ledger
computes `2^125.774016486...`) and is included in the total, and nonuniform
advice is zero bytes. Every stated bound is at least the ledger's exact
value; the score rational is the rounded-up logarithm.

The attack uses only `T=A_1` and `A_d` for `d=2^20`, both supplied at the
scoring `n`. The rules state that `d=2^20` divides `r-1` and that Cheon's
attack using these inputs is admissible. No prior submission is used.

## Inputs and the Cheon precondition

The profile fixes

    r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

(255 bits) and the G1 generator `P0` of order `r`. The factorisation

    r-1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2
              * 2508409 * 2529403 * 52437899 * 254760293^2

is verified by multiplying the factors and testing primality in the bundled
experiment, so `d=2^20` divides `r-1`. Set

    M = (r-1)/d = 500067474128019242090680508691653879
                 52509453297164571592906626415194112   (235 bits).

`zeta=7` is a generator of `Z_r^*` (the experiment checks
`zeta^((r-1)/q) != 1` for every prime `q | r-1`). Define

    zeta_hat   = zeta^d mod r      (order M),
    zeta_check = zeta^M mod r      (order d).

Every `s in 1..r-1` has the unique representation

    s = zeta^(k0 + M*k1) mod r,   0 <= k0 < M, 0 <= k1 < d.          (1)

Raising to `d` removes `k1`:

    s^d = zeta_hat^k0 mod r.                                          (2)

Dividing `s` by `zeta^k0` isolates `k1`:

    s * zeta^(-k0) = zeta_check^k1 mod r.                            (3)

Both equations are exact; there is no heuristic.

## Phases

**Phase 1.** With `m = ceil(sqrt(M))` and `m_hat = ceil(M/m)`,

    m = m_hat = 223621884914696853578958571250384896,
    m + m_hat = 447243769829393707157917142500769792,

and `m*m_hat >= M`. Write `k0 = u + m*v`, `0 <= u < m`, `0 <= v < m_hat`.
Equation (2) becomes

    [zeta_hat^(-u)] A_d = [zeta_hat^(m*v)] P0.

The attack builds the `m`-entry baby table, scans the `m_hat` giant points,
and recovers `u+mv` from a full-affine collision. The index is reduced modulo
`M` if the scan first meets an alias in `[M, m*m_hat)`; because `zeta_hat` has
order `M` and the solution `k0 < M` is unique, this recovers `k0`.

**Phase 2.** Compute `S=[zeta^(-k0)] T`, which by (3) equals
`[zeta_check^k1]P0`. With `m'=m_hat'=1024=sqrt(d)`, the same two-list search
over `{[zeta_check^(-u')]S}` and `{[zeta_check^(m'v')]P0}` recovers
`k1 in [0,d)`. Then `s = zeta^k0 * zeta_check^k1 mod r`, and the attack
verifies `[s]P0=T`. Both phases are exhaustive and deterministic, so success
probability is 1.

## Scalar multiplication: production GLV

The curve `y^2=x^3+4` has `j=0` and the endomorphism `phi(x,y)=(beta*x,y)`
with `beta^3=1 mod p`, `beta!=1`. Let

    lambda = 7^((r-1)/3) mod r,   lambda^2 + lambda + 1 = 0 mod r.

`lambda` is 128 bits. The matching cube root `beta` is selected and the
production pair is verified on the protected generator: the experiment
computes `[lambda]P0` and `phi(P0)` in the production group, checks equality,
and checks `[r]P0=O`. The short lattice
`{(a,b): a+lambda*b=0 mod r}` has the reduced basis

    (-228988810152649578064853576960394133503, 1),
    (1,   228988810152649578064853576960394133504).

The nearest-plane decomposition gives `k=k1+lambda*k2 mod r` with
`|k1|,|k2|<2^127` (verified on random scalars; 130 bits charged). Then
`[k]P=[k1]P+[k2]phi(P)` via interleaved double-and-add. Non-adjacent form
has no two consecutive nonzeros, so one multiplication uses at most 130
doublings and at most 130 additions, plus one endomorphism application.

## Batching schedule and affine prices

All `m` baby multiplications are independent (fixed base `A_d`, known scalars
`zeta_hat^(-u)`); so are the `m_hat` giant ones (fixed base `P0`). They are
processed in batches of `B=2^64` live chains. Within a batch, the 130 rounds
are run in lockstep and the `B` affine denominators of each round are
inverted together by Montgomery's trick: the prefix/back-substitution costs
`3M` per element and one inversion per round. `B=2^64` makes that single
inversion negligible per element. This is exactly the protected reference
operation, a batched incomplete affine G1 addition, `6M` with `S=M`, `A=0`,
charged **1 unit**. The workspace for a batch is

    B * (112 bytes current point + 48 bytes denominator
           + 48 bytes prefix + 32 bytes control) = 240*B = 2^71.9 bytes,

negligible against the retained table. The table itself is built batch by
batch and sorted once; the giant side is streamed in batches against the
sorted table. No giant-chain state beyond one batch is live at once.

A doubling is not the reference addition: the affine `a=0` formula is
`5M+2S = 7M`; with the protected `M=383` and `C=2000` it costs
`7*383/2000 = 2681/2000 = 1.3405` units.

## Table: merge sort and binary search, no hashing

The baby table is an array of `m` records, each a 96-byte affine point
(`x,y`, 48 bytes each) and a **16-byte index** (an index below `m` needs 118
bits, so 8 bytes would be insufficient; 16 bytes is charged). It is sorted by
merge sort with two record buffers and `L=ceil(log2(m))=118` passes.

Per entry per pass charge at most: one key comparison (8 limb comparisons and
branches, 32 instructions) plus one record move (14 load/store instructions)
plus 16 instructions of index/loop control and a 12-instruction margin, i.e.
74 instructions; over 118 passes this is 8,732 instructions = 4.37 units per
entry. The ledger charges **5 units per baby entry for the sort**.

Each giant lookup is a binary search: `L+1=119` comparisons, each 32
instructions for the key plus 16 for addressing/control, i.e. 5,712
instructions = 2.86 units; the ledger charges **4 units per giant lookup**.
There is no dictionary, hash function, load factor, collision handling or
expected-time assumption; the table is a deterministic sorted array.

## Per-multiplication overhead

The 16-unit overhead per scalar multiplication is itemised and is an upper
bound:

| item | instructions |
| --- | ---: |
| GLV decomposition (4 128-bit products, rounding, signs) | <= 1,024 |
| endomorphism `phi(P)` (one `F_p` multiplication and a store) | <= 512 |
| scalar update for the next chain (`F_r` multiplication, loads/stores) | <= 1,024 |
| loop, point selection, negation, branch and index control | <= 1,024 |
| **total** | **<= 3,584 = 1.792 units** |

The ledger charges 16 units, about nine times the itemised bound.

## Complete cost

One scalar multiplication costs

    E = 130*(2681/2000) + 130*1 + 16 = 64053/200 = 320.265 units.

Phase 1 uses `m+m_hat` scalar multiplications, `m` sort entries and
`m_hat` binary-search lookups. Phase 2 uses 2048 scalar multiplications, a
tiny sort of 1024 records (`1024*12*48` instructions) and 8 units. The fixed
allowance is `2^40` scalar instructions (`2^40/2000` units) for reading and
validating `T` and `A_d`, the factorisation, generator test, powers of
`zeta`, GLV basis, reconstruction and final verification. Hence

    W = (m+m_hat)*320.265 + m*5 + m_hat*4 + 2048*320.265
        + 1024*12*48/2000 + 8 + 2^40/2000
      = 2^126.7717970576168947035681899827184488498828103844885138...,

rounded up at the sixth decimal to `63385899/500000 = 126.771798`. Component
logarithms: phase 1 (scalar multiplications, sort and search) `126.771797`,
phase 2 `19.324`, fixed `29.034`. The bundled experiment recomputes every term
with exact integers and rationals, checks the high-precision logarithm, and
asserts the rounded score.

## Memory, preprocessing and advice

The table uses two merge-sort buffers of `m` records of `96+16=112` bytes,
i.e. `2*m*112 = 2^125.2358996` bytes. The batch workspace is `240*2^64 =
2^71.9` bytes; the fixed reserve is `2^22`. Peak memory is therefore at most
`2^125.2358997` bytes, and the claim states the rounded-up bound `2^126`.
The baby table is the preprocessing: `m*(320.265+5) = 2^125.7740165` units,
included in `W`, with no reuse across targets. Only `A_1` and `A_d` are read.
Nonuniform advice is zero bytes.

## W(n), validity and sensitivity

The attack needs `A_d` for `d=2^20`, so it applies to every `n >= 2^20`, with
`W(n)=W(2^20)` (extra powers are ignored); it is not claimed for `n < 2^20`.
The protected score is evaluated at `n=2^20`.

A conservative sensitivity check: pricing doublings at `8M/2000` and the
per-multiplication overhead at 24 units gives `E=353.16` and a total of about
`2^126.909`, still below `2^127`. The declared score uses the protected
reference price for the reference operation, the `7M` doubling count, and the
itemised 16-unit overhead.

## Evidence, checks and limitations

`experiment.py` (stdlib only) has three labelled parts:

1. **Production ledger.** Verifies the factorisation and primality of the
   factors of `r-1`, `zeta=7` as a generator, the order of `zeta_hat`, the
   production GLV pair on `P0`, `[r]P0=O`, the reduced lattice, the 130-bit
   decomposition bound on random scalars, and all cost terms, including the
   sort/search schedule; asserts the rounded score.
2. **Production-curve checks**, all with the protected `p`, `r`, `P0`.
3. **Toy end-to-end Cheon attacks**: a 29-bit curve (order 4759, `d=2`, six
   fresh secrets) and a 20-bit prime-order group (`d=1024`, four fresh
   secrets), including the geometric-progression tables and reconstruction.

The toy runs are evidence at their scope only; they are not a production-size
run or a production success measurement. The bound is analytic and the
success probability exact; the remaining limitations are the cost-model
pricing (the `6M` batched-affine reference, the `7M` doubling, `M=383`,
`C=2000`), the enormous `2^125`-byte table, and the absence of a full-size
execution. Cheon's auxiliary-input attack and the GLV method are known; no
novelty is claimed.

## Provenance

No Yukon submission is reused. The method is Cheon's auxiliary-input attack
(J. H. Cheon, EUROCRYPT 2006) and GLV (Gallant, Lambert, Vanstone, CRYPTO
2001). The target, powers, `n=2^20` and the admissibility of Cheon's attack
are those of the protected `bls12-381-g1-ecdlp` profile and rules; the cost
model is protected `classical-group-work-v2` with `C=2000` and `M=383`.
