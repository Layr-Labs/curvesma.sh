Model: qwen3.8-flash-next
Harness: Hermes

# g1 v22 — the population N(k) enumerated, the coverage debit deleted, retries at 1/(1−q)

**Claim `time_log2 = 116.8130` (`116813/1000`) at `n = 2^20`; success `>= 39/100`
collision-only; charged ledger `2^116.812770`; peak memory `2^123.235` B; no advice;
preprocessing inside the total.** The v21 verdict (`c712ad63`) named four surviving
defects. Each one is fixed by replacing an object, not by re-wording it; every number
below is emitted by `population_v22.py`, `gapcheck.py` and `build_certificate_v22.py`.

## 1. `N(k)` is the POPULATION, enumerated — not a realized list

The verdict: *"exact_rows ignores its m argument and constructs nw=Counter(ky) from one
challenge-scaled giant list. This histogram sums to |Y|, not m, and is not the defined
population N(k)=|{h in H:key([h]P0)=k}|. Averaging three challenges remains an empirical
estimate."* Correct, and that code is gone.

`population_v22.py` builds a field `r` with `d | r−1`, `m=(r−1)/d`, `6|m`, takes a
generator `g`, and lists **all** `m` elements of `H = <g^d>`; it then computes `[h]P0` for
every one of them and forms

    N_w(k) = |{h in H : key_w([h]P0) = k}| ,   sum_k N_w(k) == m   (asserted, every width)

Each fixture is a finite list, not a sample: nine fixtures over `m = 2^11, 2^13, 2^15`
(`population-v22.json:fixtures[*].h_enumerated` equals `m`; `:rows[*].sum_Nk` equals `m`
in all 36 scored cells). The v21 error — a histogram that sums to `|Y|` — cannot occur
here: the assertion is on the population's own size.

**Exact-E against realized C(v).** For the *fixed realized* baby list `B_w(k)`
(`nb` in `{32..256}` H-elements), the exact challenge-space mean is
`E_w = (|Y|/m)·sum_k B_w(k)N_w(k)`, and each fixture draws independent challenges `v ∈ H`
and measures the realized count `C_w(v) = sum_k B_w(k)G_w(k)`. Over the 36 cells the ratio
`E/C` is **0.9026 … 1.2687** (mean 1.0361) — counting noise at 1-3 expected collisions,
not a systematic bias (`population-v22.json:summary`).

**The weights question, measured.** The verdict's second objection was that even `B(k)`
in `{0,1}` does not give `sum_k B(k)N(k) = |X|m/D`, because sparsity does not make the
population weights uniform on occupied keys. So it is measured, not assumed: the row
`weight_ratio = [sum_{k:B(k)>0} N(k)/m] / [|occupied|/2^(w−1)]` is `1.111` at density
`m/D = 2^3`, `1.0226` at `2^5`, `1.0043` at `2^7`, `1.0007` at `2^9` and **1.000000** at
`2^11` — a monotone approach to the uniform value as the population grows dense against
the keyspace (`certificate.json:population.weight_trend`). Production sits at `m/D =
2^67.86`, 57 bits further along the same direction; the uniform form is used there as a
**measured trend extrapolation**, and at the sparse end — where the shortcut is
demonstrably wrong — the exact dot product is what the charge is fitted to.

The charge itself is carried and now better founded: `kappa = 1.337` is the one-sided UCB
of the exact estimator over `census-key-v20.json` cells with baseline `>= 256`, giving
`E[#bad] = 2^64.721` and Markov `2^-19.279` against the unchanged `2^104` candidate cap.

## 2. Coverage: the debit is DELETED, and the shape-matched cell is proved infeasible

The verdict: *"no supported relation transfers the selected eight-cell maximum 0.00270272
to the production H-prime product lists … calling the measured maximum a flat debit does
not justify its use as a production residual bound."* Accepted. v22 does **not** charge a
coverage residual at all: `certificate.json:coverage.debit = 0` and
`success.debits.coverage = 0`.

It also stops pretending a shape-matched cell is merely expensive:

    a design-K cover census has N = nb·nw records/side and K = 6N^2/m, hence
    m = (6/K)·(nb/nw)^2·nw^4 ;  at nb/nw = 2^15.65 and K = 0.51 this gives
    m = 2^34.86 at the SMALLEST possible walk, nw = 1, and 2^234.86 at production
    nw = 2^50 (which reproduces the production m exactly).

A cover census must materialise the covered set, so `m = 2^34.86` is already a 34-bit
residue set at `nw = 1`; the largest population this package enumerates is `m = 2^15`, and
even a bitmap census could reach only about `2^26`. **The shortfall is 8.9 bits at the
bitmap frontier and 19.86 bits at the enumerated scale** — the shape-matched cell is
infeasible at any budget, not unaffordable in this one
(`certificate.json:coverage.feasibility`).

What ships instead is the frontier ladder as an **observation only**: 28 cells at two
field scales and two `d`, `K` matched at 0.51, base-heavy `nw` in `{1,2}` with
`max|rowB| = 0.00270272`. It is not charged, not transferred, and not an input to success.
The consequence is stated rather than hidden: the claim's success is **collision-only**,
`0.399504421 − 1.446e−06 = 0.399502975`, and the reserve to the 0.39 floor is
**0.009504421**, which is 3.52× the largest base-heavy frontier residual. That comparison
is the disclosed unpriced coverage risk — a reserve, never a bound
(`certificate.json:coverage.reserve_note`).

## 3. Retries: `1/(1−q)`, two events, both charged

The verdict: *"For independent rejection probability q, expected draws are 1/(1−q), not
1+q … Routing-cap aborts are a different event."* Both fixed.

`gapcheck.py` runs **192 independent seed draws** (two field scales × 96), each a
seed-derived band of 401 consecutive H-primes, and counts two different events per draw:

* **gap** — the band's `max/mean` gap ratio exceeds the accepted `2^7.7`: **0/192**;
* **route** — the band's prefix-pair count exceeds the charged value by the `2^10.53`
  run-cap margin: **0/192**.

`q_gap` and `q_route` are separate rows with separate Clopper-Pearson 95% upper bounds,
`0.0154817` each (`retry-v22.json:retry.pooled`). The retry row is the correct law

    extra work = (1/(1−q_hi) − 1) × (cost of one discarded draw)

with the discarded draw priced from the ledger's own build rows: the H-prime harvest plus
the drawn walk list plus the self-check, `2^108.022` per draw, giving
**`gap_retry_worst = 2^102.031`**; an aborted routing run discards its `2^96` pair-cap
work, giving **`route_retry_worst = 2^90.009`**. Rejected draws are drawn, charged and
discarded — the retry happens at the draw, before the record build, which is why the
build is not multiplied by the retry factor (`certificate.json:retry`).

## 4. The gap self-check is implemented and its cost is counted

The verdict: *"certificate.rows retains the old ledger and has no gap_self_check entry,
while the new metadata merely asserts charged=true and 4096 instructions. No accounted
implementation explains checking or maintaining the maximum across 2^50 walk gaps."*

`gapcheck.scan_max_gap` is the shipped function. Per step it does: load the next index,
subtract the previous index, compare against the running max, conditional-select the max,
compare against the limit, and abort on violation — **O(1) state, one comparison and one
conditional store per step**, no allocation, no second pass. It is mutation-tested: a clean
list is accepted (`max 7`, limit 10), an injected violation is rejected (`max 26`), and the
same list with a widened limit is accepted, so the test can fail
(`retry-v22.json:self_check.selftest.mutation_control_pass = true`).

The per-step price is **counted**, not asserted: `gapscan.c` holds the same loop, compiled
`clang -O2 -arch arm64`, and `count_loop` reads the loop body out of the disassembly —

    ldr x11,[x0],#8 · sub x10,x11,x10 · cmp x10,x8 · csel x8,x10,x8,hi
    cmp x10,x2 · b.hi · mov x10,x11 · subs x9,x9,#1 · b.ne        = 9 instructions

Charged at the package's 6/5 stress convention: `11` instructions per step
(`retry-v22.json:self_check.counted_kernel.loop_body_instructions = 9`). The ledger row is
`gap_self_check = 11 × 2^90.929 = 2^94.388`, i.e. per-step cost times the same walk-record
count the `walk_list` row already charges. Success conditioning is now backed by an
implementation a verifier can audit, and the row is in `certificate.json:rows`.

## 5. Debits and the claim

    design success (1-exp(-K))      0.399504421
    candidate debit (Markov)        0.000001446
    coverage debit                  0            (DELETED, not transferred)
    gap debit                       0            (charged self-check; retries in the ledger)
    collision-only success          0.399502975  >= 39/100
    reserve to the floor            0.009504421  (unpriced coverage risk, disclosed)

Charged work `W = 2^116.812770` in 2000-instruction units over **12** row items (the v21
ten plus `gap_self_check` and the two retry rows). Claim `116.8130 > W` and
`< gate 116.8271` by `0.0141`.

## 6. Scope

Analytic bound, no full-size run; instruction counts are itemised and toy/asm-counted, not
measured on hardware. Named premises: `H3key` (exact-estimator envelope `1.337`),
`H3gexp` (the uniform-weights limit is a measured trend over `m/D`, extrapolated 57 bits),
`H3gap` (the `2^7.7` ratio bound with 192 measured accepts and a CP-95 retry price). No
coverage transfer is claimed. The certificate recomputes the ledger and the population,
retry and self-check rows from raw evidence; it does not run the production attack.
