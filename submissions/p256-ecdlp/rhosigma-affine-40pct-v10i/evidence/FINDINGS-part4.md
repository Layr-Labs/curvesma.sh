# p256 v10f part 4 — the attractor is an INVOLUTION, not an index bias; killed; census 0.4217 CP-95

## The diagnosis was right about the disease, wrong about the cure

Part 3 blamed the 91.4% self-hit rate on "the 2048-row map has a global attractor".
It does. The mechanism is exact, and it is not fixable by any index rule:

> **Involution.** With one shared table, `x(-S) = x(S)` ⇒ both representatives of a
> class select the same row `T_j`. One step `S → canon(S + T_j)` lands on
> `±(S + T_j)`; if the next index repeats, the next step adds the *same* `T_j` to the
> negated representative: `canon(-(S+T_j) + T_j) = canon(-S) = S`. Two steps, back to
> the start. A 2^11-row table repeats an index within ~2^11 steps by the birthday, so
> **every** walker is in a 2-cycle long before the design budget, and every later
> record is a self-hit.

Verified three ways:

1. **Exact bigint** (`diag_cycle.py`, no C field code): 6/6 random starts cycle with
   length 2, entry μ = 140…3,381, closure relation `q + T[j'] == -p` on all six.
2. **C probe on the fixture** (`probe_basin.c`, 500 starts, cap 5,000 steps): 62.0%
   of starts cycle, **median cycle length 2**, 40.2% of a 5,000-step trail repeats.
3. **The four candidate fixes from the brief, measured** — O2 (decorrelated index),
   O1 (per-walker coset offset) and an xor-fold are all *identical* to the baseline,
   because the involution only needs `idx` to agree on `±S`:

| mode | index rule | cyc. within 5,000 | med. cycle | trail repeats |
| ---: | --- | ---: | ---: | ---: |
| 0 | `x & 2047` (part-3 baseline) | 62.0% | **2** | 40.2% |
| 1 | `mix64(x) & 2047` (O2) | 59.2% | **2** | 38.5% |
| 3 | `(x + O_w) & 2047` (O1, per-walker coset) | 60.0% | **2** | 39.4% |
| 4 | `(x ^ (x>>20)) & 2047` | 63.4% | **2** | 46.7% |
| **2** | **two tables, parity + salted mix (O3)** | **0.0%** (0/500) | none | **0.0%** (0/1.5M) |
| 5 | no canon, `x & 2047` (fallback) | 0/400 @ cap 40,000 | none | 0.0% |
| 6 | no canon + mix | 0/400 @ cap 40,000 | none | 0.0% |

`eqv_bad = 0/2000` and `odd_y = 0` in every mode. **The brief's O2 and O1 are dead on
measurement, not on argument** — worth recording because both look reasonable.

## The fix that shipped: O3, two parity-alternated tables

`j = mix64(x ^ (par ? S_B : S_A)) & 2047`, row from `T^A` on even steps and `T^B` on
odd (independent logs, cross-table duplicate-checked, 2048 even-y rows each). A
2-cycle now needs `T^A_j == T^B_j'` as group elements (~2^-40), not an 11-bit index
collision. Canon retained ⇒ the (a,c) solve and its 0/2,000 class gate are untouched;
storage 2^22 → 2^23 B = 2^-50 of the 2^73 B cap; counted step 1027 + **26** = **1053**.

Fallback on file: the no-canon walk (modes 5/6) is equally clean and 26 instructions
cheaper, but moves sign handling into the solve path — not shipped.

## Census (2^40 fixture, 16 slots, batched shared inversion)

| Λ | T | trials | recovery | CP-95 | merges = verified | self-hits | dc0 | check_fail |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.53 (40% design) | 763,376 | 4,000 | 0.382750 | 0.370050 | 14,864,870 | 0.123% | 0 | 0 |
| **0.60** | 812,224 | **8,000** | **0.430875** | **0.421725** | 37,839,675 | 0.126% | 0 | 0 |
| 0.53 dp2 | 763,376 | 2,000 | 0.263500 | 0.247311 | 1,936,218 | 0.297% | 0 | 0 |

Self-hits fall from **91.4% → 0.12%** (700×), recovery rises **0.0027 → 0.4217**
(156×), and every merge still verifies (`[x]P == Q`, conversion 1.000, 52.7M
verified merges across cells). Model cross-check: `model_cross` 0.4302 vs measured
0.4309 at Λ = 0.60 — the birthday model now reproduces to 0.2%.

**The 0.39 bar is met at Λ = 0.60 (CP-95 0.4217) but NOT at the 40% design budget
Λ = 0.53 (CP-95 0.3701).** Both are declared in the package; the 40% cell is 0.02
short and is stated, not hidden.

## Price

charge = 1053 + (1+2^-10)·0.639140 = **1053.639764** (Δ = +26 counted: 5 mix_key,
2 mul, 10 load/re-spill, 5 move_extend, rest 1s; `umulh = 114`, `mul = 94`);
W = A + log₂(charge) = **127.617822**, claim = **127.6179** (shipped-.s absolute
count 1055 would give 127.6206; both in `evidence/ledger.json`).

## Gate (re-polled 24 Oct 04:31Z, YUKON_API_TOKEN unset, `--all --json`)

57 submissions, 12 in review. Lowest live = **b89d1d72 / 127.370918** (Randomblock1,
review); then jungjipdo 127.86 (c79cc19a), 128.06 (2ec1a0b8). Claim 127.6179 is
**+0.2470 above the head** → not submittable; an armed submit would be below the
gate and would burn the package.

## Verdict — HOLD (not staged, not armed)

* Package built and validated: `curvesmash/submissions/p256-ecdlp/rhosigma-affine-40pct-v10f/`
  (17 files, 88,768 B ≤ 262,144 B; Model/Harness lines 1–2; compact JSON; rev history).
* **Not staged, not armed**: the main repo has *no commits* (submissions run from
  `publication-staging/*` worktrees), so staging v10f is a worktree-create + commit
  operation that did not fit the cap. `p256_gate_trigger.HELD` stays.
* To arm later, in order: (1) create the staging worktree and commit this package;
  (2) write `p256_gate_trigger_v10f.sh` from the secp pattern with
  `UNLOCK=127.6234`, `cd` into the staged package, `--note-file submissions/p256-ecdlp/rhosigma-affine-40pct-v10f/public-note.md`;
  (3) `rm p256_gate_trigger.HELD`; (4) `nohup … & disown`, verify pid + first poll.
* Rejection of b89d1d72 moves the head to 127.86 ≥ 127.6234 → the trigger would fire
  and 127.6179 would lead. That is the whole upside case.

Artifacts: `probe_basin.c`, `diag_cycle.py`, `census_v10f41.c` (+ `census41_driver.sh`),
`step_v10f_o3.c/.s`, `build_v10f.py`, `evidence/census41-{dp0,dp0-lam60,dp2}.jsonl`,
`basin_m{0..6}.json`, `sc_{ref,o3}ELF.json`, package copy under
`curvesmash/submissions/p256-ecdlp/rhosigma-affine-40pct-v10f/`.
