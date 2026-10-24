# argument.md — p256 v10f (O3): why the walk now mixes, and what it certifies

## 1. The map and the theorem that kills it

The walk is on negation classes of `E: y² = x³ + 2` over `p = 1099511629627`
(`q = 1099513640503` prime, cofactor 1). A step is

    S -> canon(S + T_j),   j = idx(x(S)),   T_j = [t_j]G  (even-y row)

with `canon` the even-y representative. `x(-S) = x(S)`, so `idx` is a class
function and the quotient walk is well defined. That was v10f's repair. It is also
the defect:

> **Claim (involution).** If `idx(S_{n+1}) = idx(S_n)` then `S_{n+2} = S_n`.

*Proof.* `S_{n+1} = ±(S_n + T_j)` with the sign fixed by `canon`. The index depends
on `x` alone, so it is the same for the positive and negative representative; hence
the next step uses the same `T_j`. If the sign was negative, `S_{n+2} = canon(-(S_n + T_j) + T_j) = canon(-S_n) = S_n`. ∎

A 2^11-row table repeats an 11-bit index within ~2^11 steps by the birthday, so
every walker is trapped in a 2-cycle after O(2^11) steps regardless of `N'`. Two
consequences, both measured:

* the self-hit fraction grows with the budget and the recovery *falls* as `q` grows
  (the walk is merely longer inside the same cycle): 0.00425 at 2^39 classes,
  2.8% self-hits at 200 steps on the old 2^24-class fixture;
* **no index rule can fix it** — every rule that is a function of `x` gives the
  same involution, because the involution needs only that `idx` agree on `±S`.

Exact-bigint confirmation (`experiments/diag_cycle.py`, independent of the C field
code): 6/6 random starts, cycle length 2, entry 140–3,381, and the closure relation
`q + T[j'] == -p` verified on every cycle. C confirmation on the fixture
(`experiments/probe_basin.c`, 500 starts, cap 5,000): 62.0% of starts cycle with
median length 2 and 40.2% of a 5,000-step trail is a repeat; `mix64(x) & 2047`,
`(x+O_w) & 2047` and xor-fold measure the same (59.2/60.0/63.4%, all length 2).

## 2. Two tables break the involution

Tables `T^A`, `T^B` are independent (distinct logs, cross-table duplicate-checked,
2048 rows each, even-y). The step at parity `par` uses

    j = mix64(x ^ (par ? S_B : S_A)) & 2047,   row = T^{A if par==0 else B}[j].

The composite two-step map is `S -> ±(±(S + T^A_j) + T^B_j')`. A 2-cycle now needs
`T^B_j' = T^A_j` as group elements, not an index collision; for independent tables
that is a ~2^-40 event. Measured (`evidence/basin.json`): **0/500 starts enter a
cycle within 5,000 steps, and 0 of 1.5M trail steps repeat** (the head's rule:
62% / 40%). The no-canon variant (walk on true points, solve both signs at the
record) also measures 0/400 at cap 40,000 and is retained as the cheaper fallback;
it is not shipped because it moves the sign handling into the solve path, and the
canon path is already gated at 0/2,000 splits.

Storage for the second table is 2^22 B → 2^23 B, i.e. 2^-50 of the declared 2^73 B
cap: free at the declared precision. The salting mix costs 26 counted instructions.

## 3. Census

Fixture arithmetic is gated: 2,000/2,000 random canonical additions match an exact
bigint re-implementation, every input and output on the curve; the 16-slot batched
shared-inversion step is bit-identical to the one-inversion-per-slot path over 6
steps × 16 slots, and both pass the `(a,c)` invariant `state == [a]P + [c]Q`.

| Λ | T (steps) | per walker | trials | recovery | CP-95 | merges = verified | self-hits | dc0 | check_fail |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.53 | 763,376 | 47,711 | 4,000 | 0.382750 | 0.370050 | 14,864,870 | 0.123% | 0 | 0 |
| **0.60** | 812,224 | 50,764 | **8,000** | **0.430875** | **0.421725** | 37,839,675 | 0.126% | 0 | 0 |
| 0.53 dp2 | 763,376 | 47,711 | 2,000 | 0.263500 | 0.247311 | 1,936,218 | 0.297% | 0 | 0 |

Conversion is 1.000 in every cell: each merge is re-verified by recomputing `[x]P`
and comparing the full point with the target, so a merge that does not solve is not
counted. `dc0 = 0` (no degenerate denominators), `check_fail = 0`. The DP predicate
is the declared one, `mix64(x ^ DPSALT) mod 2^61 == 0`.

The declared design budget is the **Λ = 0.60** cell: `T = 812,224` steps, which is the
configuration this census certified, CP-95 **0.421725** — that is the package's floor and
it is above the 0.39 bar. The Λ = 0.53 cell (`T = 763,376`) is the same walk at **40% of
the design budget**: it certifies 0.370050 and is reported as a **sensitivity**, not as a
design, production or pricing row. Both are declared, the smaller budget disclosed rather
than hidden, and no number is selected after the fact.

## 4. Price

`evidence/stepcount.json`: `step_v10f_o3.s` counts **1,055** by the clang CFG
longest-iteration method; the same toolchain rebuilds the v10f anchor source at
1,029 against its shipped 1,027, so the same-toolchain delta is **+26** and the
charged step is the anchor plus delta = **1,053** (`umulh = 114 = 6 × 19`,
`mul = 94`, rows machine-classified and summing to the count). Declaring 1,055 with
the shipped `.s` instead would give W = 127.620558 and claim 127.6206; the package
declares the anchor-relative 1,053 / 127.6179 and states both.

`charge = 1,053 + (1+2^-10)·0.639140 = 1,053.639764` (the rare-event rebuild and
its d2 martingale reserve are inherited from the v10f ledger unchanged; the fix adds
no rare path). `W = A + log₂(charge)`, `A = 127.370918 − log₂(887.904982)`,
`W = 127.617822`, `claim = ⌈10⁴(W + 2×10⁻⁵)⌉/10⁴ = 127.6179`.

## 5. What this package does not claim

* It does not claim 0.39 at 40% of the design budget: at `Λ = 0.53` the measured floor
  is 0.370050 and that number is stated rather than hidden. Its 0.39 claim is at the
  declared design budget `Λ = 0.60`, where the certified floor is 0.421725.
* It does not claim the v10e floor 0.401692 — that figure came from a kernel whose
  shipped header transposed `p`'s middle limbs (re-audited RED in this lane).
* The 0.421725 floor is a finite-sample bound on this fixture's class count
  (`2^39.000` classes against P-256's `2^254.999999999664`); it is a recovery
  measurement of the declared mechanism, not a proof about the field.
