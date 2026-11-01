# Unit-entry linear-algebra kernel

Base charge (STORAGE.md): 88 instructions per sparse entry = 33 (`la_mac.s` multiply-accumulate)
+ block selection 4 + sign 3 + transfers 32 + finishing/row headers 16.

**Change.** After filtering almost all coefficients are +-1. A +-1 entry needs no multiplication: the
4-limb element is added into the positive or negative five-limb partial sum of its (block, sign) group,
which the base layout already keeps separately. Each group stores a unit list (index only) and a general
list (index and 8-bit magnitude 2..255). `group_run.s` is clang 18 output (`--target=arm64-apple-macos14
-O3 -fno-vectorize -fno-slp-vectorize -mgeneral-regs-only`) for one group: unit loop `LBB0_1` is
**11** instructions per entry (ldr, add, 2 ldp, adds, 3 adcs, cinc, subs, b.ne); general loop
`LBB0_3` is **33**, identical to `la_mac.s`. The partial sums stay in registers across both loops.
Outside the loops the function has 13 instructions against 9 for `mac_run`; with two more header words
per group, the extra fixed cost is at most 6 per group. Columns occupy 14 blocks (STORAGE.md), so a row
has at most 28 groups: at most 168 extra instructions per 101-entry row, below 1.67 per entry, charged 2.

Per-entry charge at charged unit fraction f: `4 + 3 + 32 + 16 + 11 f + 33 (1 - f) + 2 = 90 - 22 f`;
at **f = 0.88** this is 70.64, charged **71** (`C_ENTRY`). With the cheaper product the Block Wiedemann
cost of proof.md section 4 is minimised at **(m, n) = (22, 4)** instead of (26, 5); the failure bound
4N/r < 2^-192 does not depend on m, n. Linear algebra 2^119.36845 -> **2^119.1608**, generator
2^106.4172 -> **2^106.1537**.

**H4u (heuristic): at least 88% of the production filtered matrix row weight is +-1 entries.**
`unit_entries.py` reruns the package's filtering simulation (`filtsim.py`: signed sigma-quotient system,
mean norms 748.39/498.93 bits, B = 2^69, excess 1.3, purge, Markowitz merge to mean weight 100,
|coefficient| <= 255) with an unchanged copy of `filtsim.merge` plus counters, and measures the unit share
of row weight when the mean weight reaches 100. Explicit entries use exact coefficients; filtsim's
implicit mid-range weight (`r[-1]`) counts as unit unless it entered a row from a pivot with |m| > 1.

Measured unit fractions: 2^12 seeds 1-3: 0.9330, 0.9350, 0.9173; 2^13 seeds 1-2: 0.9338, 0.9405; 2^14 seeds 1-3: 0.9452, 0.9376, 0.9452; 2^15 seed 1: 0.9461.

The fraction rises with scale; the charge 0.88 is at least 0.037 below every value. Like N/U (H3) this is
an extrapolation to production. A shortfall changes cost only (22 instructions per missing unit entry),
never correctness. `evidence16.py` checks both loop lengths and `C_ENTRY >= 88 - 33 + 0.88 x 11 +
0.12 x 33 + 2`.
