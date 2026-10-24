# argument.md — p256 v10g: kernel repaired and exercised, theorem restated with its sign condition, price re-derived

Revision for verdict `b8b00688`. Three items are *repaired artefacts* (kernel, kernel test,
theorem), two are *withdrawals of unsupported scope* (success transfer, budget label and
anchor), and one is a *re-measurement* (instruction count with and without the defect).

## 1. The kernel was wrong; it is now right and it is now tested against bigint

`p` = 2^256 − 2^224 + 2^192 + 2^96 − 1 has little-endian 64-bit limbs

    P0 = 0xffffffff_ffffffff   P1 = 0x00000000_ffffffff
    P2 = 0x00000000_00000000   P3 = 0xffffffff_00000001

v10f shipped `P1 = 0`, `P2 = 0x00000000ffffffff` — limbs 1 and 2 transposed. The review is
right that this was not a cosmetic slip: the CIOS reduction uses `p2 = 0` (it adds the carry
into limb 2 instead of a product `q*p2`), and the final conditional subtraction omits limb 2
as well, so the shipped header could not have computed in the field it claimed. `fe_inv`'s
exponent table already carried the correct limbs, i.e. the artefact contradicted itself.
v10g sets `P1 = 0x00000000ffffffff`, `P2 = 0`, masks `P2` in `fe_sub`'s add-back (it was
unmasked), and replaces the hard-coded `0` in `fe_mul`'s final subtraction with `P2`.

`fe_inv` was *also* wrong, for an independent reason found while writing the test: binary
exponentiation in a Montgomery product accumulates one factor `R^-1` per call and squaring
**doubles** that factor, so the v10f routine returned `a^(p-2) * R^-j` with

    j = 0x1ffffffff00000001000000000000000000000000fffffffffffffffffffffffc

which is not an inverse in this representation. v10g appends one `fe_mul` by
`R^(j+3) mod (p-1)` and states the law it now satisfies:
`fe_mul(a, fe_inv(a)) = R` (Montgomery one) for every `a != 0`. `fe_inv` is not on the
counted walk path (the step uses the shared batch inversion over `PROD`/`SUF`), so this
repair moves no count.

**The test that was missing.** `evidence/p256_kernel_test.c` includes the *shipped*
`evidence/p256_fe_cc.h` and the *shipped* counted step source (for `canon_y`), emits
10,240 vector records (10,000 random + a 15×15 edge grid + 15 self-pairs), 215 inversion
records and 315 canonicalisation records, and `evidence/p256_kernel_test.py` recomputes every
field in Python big integers at the real P-256 modulus. Edge values exercised: 0, 1, 2,
p−1, p−2, (p−1)/2, 2^63, 2^64, 2^128, 2^192, 2^224, 2^255, 2^255+1, 2^256−2^224, p−1−2^96.

| build | header | mismatches |
| --- | --- | ---: |
| fixed | `evidence/p256_fe_cc.h` (v10g) | **0 / 10,770 records** |
| mutation control | `evidence/p256_fe_cc_broken.h` (v10f, transposed) | **25,772** |

The control is the v10f header verbatim (only an include guard added), compiled from the same
driver by swapping the include directory; it fails on mul, add, sub and inv, so the test
demonstrably detects the defect it was written for. `evidence/kernel_test.json` holds the
counts, the two compile command lines and both header digests.

## 2. Involution: the theorem as stated needed a sign condition, and it no longer claims one

The v10f claim was "if `idx(S_{n+1}) = idx(S_n)` then `S_{n+2} = S_n`". That is false as a
universal statement, exactly as the review says: a repeated index makes the second step add
the same `T_j`, so the second state is `canon(S + 2T_j)`, which equals `S` only when the first
canonicalisation produced the **negated** representative.

Restated, and this is now the only form the package uses:

> **Lemma (conditional 2-cycle).** Let `S_{n+1} = canon(S_n + T_j)` with `idx(S_{n+1}) = idx(S_n) = j`.
> If `canon(S_n + T_j) = -(S_n + T_j)`, then `S_{n+2} = canon(-(S_n + T_j) + T_j) = canon(-S_n) = S_n`,
> and the pair `{S_n, S_{n+1}}` is a 2-cycle.
> **The sign condition `canon(S+T_j) = -(S+T_j)` is required and is not implied** by
> `x(-S) = x(S)`, by the index rule, or by a birthday repetition of an index.

What the package therefore claims about cycles is *measured on the fixture*, not universal:
on the 2^40-class fixture, a single shared table puts 62.0% of 500 starts into a cycle with
median length 2 (index rule `x & 2047`; `mix64`, per-walker coset and xor-fold read
59.2 / 60.0 / 63.4%, all median 2), 40.2% of a 5,000-step trail repeats a class, and the
independent bigint re-derivation (`evidence/diag_cycle.py`, 6/6 starts) reproduces the
2-cycle and the closure relation. The universal "every walker is trapped" sentence is
**withdrawn**: it is not proved, and the negating branch is measured behaviour, not an
identity. The v10f prose's consequence — "no index rule can fix it" — is likewise downgraded
to "every index rule we measured behaves the same way on this fixture".

The two-table construction is unchanged, and its justification is also the honest one:
alternating two independent 2048-row tables means a 2-cycle through the sign-condition route
needs `T^B_j' = T^A_j` as **group elements**, not an 11-bit index collision. Measured:
0/500 starts cycle within 5,000 steps and 0 of 1.5M trail steps repeat
(`evidence/basin.json`, mode 2), against 62% / 40% for the single-table walk.

## 3. Census: what it certifies, and what it does not

| Λ | T (steps) | per walker | trials | recovery | CP-95 | merges = verified | self-hits | dc0 | check_fail |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.60 (design) | 812,224 | 50,764 | **8,000** | **0.430875** | **0.421725** | 37,839,675 | 0.126% | 0 | 0 |
| 0.53 (sensitivity) | 763,376 | 47,711 | 4,000 | 0.382750 | 0.370050 | 14,864,870 | 0.123% | 0 | 0 |
| 0.53, DP-2 control | 763,376 | 47,711 | 2,000 | 0.263500 | 0.247311 | 1,936,218 | 0.297% | 0 | 0 |

Conversion is 1.000 in every cell: each merge is re-verified by recomputing `[x]P` and
comparing the full point, so a merge that does not solve is not counted.

Three scope statements the review asked for, in the place the next reader will look:

* **The certified row is the every-step-record configuration (`dp_bits = 0`).** The
  production predicate is a sparse filter (`mix64(x ^ DPSALT) mod 2^61 == 0`); the only
  sparsification measured is the DP-2 control, and it *loses* a large fraction of the floor
  (0.247311 against 0.421725). The package therefore claims its 0.39 floor for the certified
  every-step configuration on the fixture, and **does not claim** 0.39 for the sparse-DP
  production predicate: the transfer from dense records to a 2^-61 filter is unproven, and
  the DP-2 row is the measured direction of that loss.
* **The budget label was wrong.** 763,376 steps is *94.0%* of the design step budget
  (763,376/812,224 = 0.939859), not "40% of the design budget". The "40%" phrasing came from
  Λ being at 88.3% of the design Λ (0.53/0.60); in steps the two cells are 6% apart. The row
  is now labelled by its actual fraction, and it remains a disclosed sensitivity.
* **The fixture is a 2^40-class prime-order curve** against P-256's 2^254.999999999664 classes.
  Both regimes are pre-capture (first revisit ≈ √(classes) ≫ budget), which is the premise
  under which the measured coverage statistic is transferred; it is a declared premise, not a
  proved transfer.

## 4. Price: the anchor is gone, the count is re-measured, and the claim is the review's own reading

v10f priced with an anchor `A = 127.370918 − log2(887.904982)` and a fixture-sized budget
`T = 812,224`, and charged 1,053 instructions obtained by differencing a re-built anchor
against its shipped 1,027. The review rejected both moves, correctly: the anchor's derivation
was never supplied, and the submission's own `T` is the fixture's.

v10g prices the way the review priced it — no anchor, no fixture budget:

    T_prod = sqrt(2 · Λ · N'),  Λ = 0.60,  N' = (r+1)/2,  r = P-256 group order
    log2 T_prod = 127.631517202749
    charge = 1,101 counted + (1 + 2^-10)·0.639140 reserve = 1,101.639764
    W = log2 T_prod + log2(charge/1000) = 127.771170
    claim = 127.78  (= 6389/50, ≥ W + 2×10^-5, headroom 0.0088 bits)

This is deliberately *above* the review's own recomputation of the v10f package
(127.7069 with the 1,053.64 charge), so the claim no longer contradicts the reader's reading
of our own numbers; it is also above every count variant (1,053 → 127.62; 1,055 → 127.68;
1,101 → 127.77). The 0.639140 reserve and its `(1+2^-10)` multiplier are carried from the
v10f ledger unchanged; the repair adds no rare path.

**The count, and the delta proof.** `evidence/stepcount.json` now carries two builds of the
*same* source `evidence/step_v10f_o3.c` from the same compiler invocation, differing only in
which header is on the include path:

| build | header | counted instructions | umulh | mul |
| --- | --- | ---: | ---: | ---: |
| v10g | `p256_fe_cc.h` (repaired) | **1,101** | 114 | 94 |
| control | `p256_fe_cc_broken.h` (v10f) | 1,054 | 114 | 94 |

The delta is **+47 instructions, and it is not zero**: the transposed header let the compiler
fold `q*p1 = q*0` and the `p2` add out of every reduction round, so the v10f tally priced a
kernel that did less work than P-256 needs. The 47 extra instructions are all constant
materialisation and carry-chain work (`mix_key` 17 → 65); the multiply counts `umulh = 114`,
`mul = 94` are identical, so the arithmetic schedule is unchanged and only the modulus
constants now cost what they should. The earlier 1,053 figure (anchor-differenced) is
superseded by this measurement and is not used.

## 5. Kernel evidence vs mechanism evidence

The P-256 kernel supplies the *price*: it is the arithmetic the counted step executes. It does
not supply the *success* evidence — the census runs on a 2^40-class fixture with its own
40-bit field path (`evidence/census_v10f41.c`, gated 2,000/2,000 against bigint and 6 steps ×
16 slots bit-identical between the shared batch inversion and one inversion per slot). The
review's complaint that the bigint comparisons "concern separate 40-bit fixture arithmetic"
is answered for the kernel by §1 and is accepted for the census by §3: the two are now
explicitly separated, and neither is allowed to certify the other.

## 6. What this package does not claim

* It does not claim the universal involution theorem; it claims the conditional lemma of §2
  plus the measured cycle census.
* It does not claim 0.39 for the sparse-DP production predicate; it claims 0.39 for the
  certified every-step configuration and discloses the DP-2 penalty.
* It does not claim a production-length transfer of the recovery statistic; the pre-capture
  premise is declared, not proved.
* It does not claim the v10e floor 0.401692, whose kernel header is the one repaired here.
* It does not claim any authenticated prior acceptance as support; `dependencies.json` is `[]`
  and every premise above is derived or measured inside this package.
