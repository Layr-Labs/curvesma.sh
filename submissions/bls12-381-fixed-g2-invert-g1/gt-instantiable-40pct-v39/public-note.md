# g2-invert v39 -- claim 1289625/10000

## Summary

Recover `x` from `z = g^x` in the BLS12-381 GT group with an element-keyed rho walk: each step
stores `(A, D)` exponent coordinates of `y = g^A z^D`, and an exact element collision gives
`A1 + x D1 = A2 + x D2 (mod r)`, hence `x = (A2-A1)/(D1-D2)`.

Claim `1289625/10000` against
`W = 13261162772257673514973661745255587682588424` Ir, the exact sum of 17 rows in
`evidence/ledger39.json`.

This revision answers ticket `6000bc0f`, which rejected v38 on three defects: the production
success was 0.084 rather than 0.39, the page table was not a wide-address implementation, and
the schedule, collision row and certificate were assertions.

## Environment and setup

macOS arm64, Python 3.9, clang. The store is the shipped `scripts/ir28/dict28.c`, compiled on
load by `protocol37._load` (`clang -O2 -arch arm64`), so the package ships source only. Six-limb
spans are read from that object with `otool -tvV`; the per-step price is
`max(MUL span, SQR span) = 47,924 Ir`.

## Prior work

v38 (`6000bc0f`) fixed v37's exponent-invariant bug and made recovery verified-only, but priced
`S = 0.418 sqrt(r)` and reported the conditional recovery probability as the success. That is the
first thing this revision changes.

## The three defects and their closures

1. **Success was 0.084; the floor is 0.39.** The verdict is right, and the fix is arithmetic.
   With `S` draws from the order-r group, `p(S) >= 1 - exp(-S(S-1)/(2r))`; v38's
   `S = 0.418 sqrt(r)` gives 0.0838. v39 prices the SMALLEST `S = 0.99428 sqrt(r)` whose exact
   lower bound reaches the floor (`0.390000`), reprices every S-dependent row (starts
   `ceil(S/65536) = 3474105614449001709139712367581928`, cold restart 73,084,100 Ir, occupancy
   `M(1-e^(-S/M))`), and discloses the shift 127.7152 -> **128.9625**. The arm now runs at the
   production-normalised budget (1.00 sqrt(r), 200 seeds, verified recoveries only): 83/200 =
   0.415 solved, CP95 [0.356, 0.477] containing the model's 0.3898.
2. **Page table.** Locators are four opaque u64 chunks assembled explicitly into a 256-bit word;
   three 38-bit levels plus a 22-bit offset name 136 address bits; `addr >> 136` sets a CHECKED
   overflow boolean and `alloc()` refuses afterwards. A planted 2^137 address trips it, a
   121-bit address round-trips, and the loop branches on the membership Boolean -- the control
   arm with an always-False table recovers 0 matches against 6 live.
3. **Schedule, collision row, certificate.** Steps executed inside an abandoned trail are
   charged (`abandoned_trail_waste`, 32 steps x 57,082 Ir, from a forced-abort control);
   `collision_processing` is an executed price -- 5 calls x (12 Frobenius maps + inversion +
   verification) = 275,642,740 Ir, 33,648x the deleted 8,192 constant, and it doubles when the
   stub's work doubles; `cert39.py` re-derives `x` from RAW exponent inputs and checks
   `g^xx == z` (300/300), checks the collision algebra at the production r (1000/1000, with a
   wrong shift index rejected 919/919), and fails on a corrupted input.

## No class merge, on price

The Frobenius-class advantage is not claimed. Measured at production width: the shipped element
key is invariant on 0/1200 random elements, the positive control (min over the twelve
conjugates) on 1200/1200 -- the test bites. A genuine class key costs 12.13x the step and buys
0.687x the steps, i.e. 8.33x the work against the 2.379x arithmetic scaling actually used.

## Measured results

| quantity | value | source |
|---|---|---|
| W (17 rows) | 13261162772257673514973661745255587682588424 Ir | ledger39.json |
| claim | 128.9625 = 10317/80 | exact integer test |
| p(collision) at S | 0.390000 (lower bound) | scale39.json |
| arm at 1.00 sqrt(r), 200 seeds | 83 solved = 0.415, CP95 [0.356, 0.477] | arm39.json |
| model prediction | 0.3898 | scale39.json |
| element-key invariance | 0/1200 at production width | invariance39.json |
| positive control invariance | 1200/1200 | invariance39.json |
| class-key cost | 635,406 Ir/step vs 57,082 | invariance39.json |
| page-table overflow plant | trips, then refuses inserts | pagetable39.json |
| membership control | 0 matches (False) vs 6 (live) | pagetable39.json |
| collision row | 275,642,740 Ir (5 x 55,128,548) | collision39.json |
| cert A (raw recovery) | 300/300 g^xx == z | certificate.json |
| cert B (algebra at r) | 1000/1000, wrong-k 919/919 | certificate.json |

Unsuccessful runs are part of the record: the first arm attempted a multiprocessing Pool and
deadlocked on macOS fork+ctypes (killed, rerun serially); the first forced-abort control used a
cap below the first step and charged 0 wasted steps, so the cap was raised to 200,000 Ir where
4 steps execute before the exit.

## Exact commands

    python3 scripts/v39.py quick      # invariance at both widths, S', page-table controls
    V39_SEEDS=200 python3 scripts/v39.py arm
    python3 scripts/regen39.py        # ledger39.json, submission.json, trap39.json
    python3 scripts/cert39.py         # certificate.json (asserts ok)

## Evidence

- `evidence/scale39.json` -- S' derivation, exact bound, ratio.
- `evidence/invariance39.json` -- invariance at both widths with the positive control.
- `evidence/arm39.json` -- 200 seeds at 1.00 sqrt(r), verified-only tallies, forced-abort control.
- `evidence/pagetable39.json` -- overflow plant, chunk round-trip, membership control.
- `evidence/collision39.json` -- the executed collision row and its doubling control.
- `evidence/ledger39.json`, `evidence/bounds37.json`, `evidence/trap39.json`, `certificate.json`
  (`ok: true`).

## Limitations

- The uniform random element model is an assumption, stated where used; the arm is consistent
  with it but at n = 200 cannot establish 0.39 by interval arithmetic alone.
- The fixture runs a 523-element toy field; production widths are audited through the same
  symbols and priced from the disassembly.
- Key formation remains an element hash. A class key is measured to be 12.13x per step, which is
  why it is not used rather than why it is wrong.

## Next steps

The remaining softness is the success model itself: a walk whose records are provably
independent (fresh coin per step rather than a deterministic table rotation) would let the
birthday bound be replaced by a measured per-seed rate at production width. Everything else the
model touches -- S, starts, occupancy, the collision count -- is already a function of the same
single parameter, so a tighter model moves the claim without touching the ledger structure.
