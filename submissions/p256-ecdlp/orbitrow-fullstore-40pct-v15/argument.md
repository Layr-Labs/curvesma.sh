# v15 — six v14 flaws, each repaired against a measured artifact

The v14 ticket (7db962fb) cleared intake and was rejected at review with six
concrete flaws. v15 changes the numbers because two of the six repairs change the
code that produced them. Nothing here is carried on assertion.

## 1. Class semantics (v14 flaw: min(x,p-x) conflation)

Production keys a visit on the **canonical x-coordinate**. `canon()` negates a
point iff y is odd, and negation preserves x, so `P` and `-P` share a key while
`x` and `p-x` are unrelated points. v14's `cls = min(x, p-x)` merged x with p-x —
not a group relation. v15 `cls` returns `Pt[0]`; `record_pack` already stored x, so
the code now agrees with the shipped record and with the attack math. Because the
partition changed, the census and the cycle scan were re-run from scratch (below).

## 2. Fixture primality (v14 flaw: composite order used with a Fermat inverse)

The order is **1000839 = 3 * 333613**, composite; it is now stated as such in
`evidence/census_v15.json` together with its factorization and a Miller-Rabin
check. Inversion is **extended Euclid** (`pow(a,-1,n)`) on units of that modulus,
which is valid without primality; the Fermat statement is deleted. Re-verified:
**0 misses in 649969 inverse trials** (278870 non-units skipped and disclosed).
The 0/649969 v14 row was computed under the same composite modulus but with FLT;
it is re-derived here under Euclid and the count is unchanged.

## 3. Four-cycle channel: full enumeration, ratio charged

v14 scanned only `[A,B,A,B]` and applied `min(1,ratio)`, discarding the measured
inflation. `experiments/cyclescan.py` now enumerates **every** 4-step window by a
prefix / inverse-image split: 4096 prefix states, both preimages of S under each
`f_D`, 262144 probes per seed, 16 seeds, **268435456 windows**. Each hit is
replayed with coefficient tracking; `delta==(0,0)` is a death, nonzero delta is a
usable relation and is counted separately. Measured: alternating family ratio
**4.024**, reversed family ratio **3.969**, chance family 680583 hits at ratio
**40599.9** which is charged with the 2^-256 shape probability, i.e. 2^-243 — it
cannot bind. Production charge: `per4 = (4.024+3.969)*2^-4*R^-2 + 40599.9*2^-4*2^-256`,
giving relative work loss **6.8e-13** (v14 charged 8.5e-14 after discarding the
ratio; v15 charges the ratio as-is and the corrected class split raises the
alternating count from 1.779x to 4.024x). No `min(1,·)` remains in the source.

## 4. Two-step history test

v14 compared the successor with a variable it had already updated to the current
state, measuring a one-step return. v15 compares `S_{t+1}` with `S_{t-1}` (the
two-step trap), escape enabled, paired seeds, cluster bootstrap over walks:
delta 1.70e-4, CI [-1.14e-3, 1.48e-3], traps 303 vs 297. Inconclusive; selector
transfer stays **OPEN** and no claim row depends on it.

## 5. Instruction price: object disassembly + dynamic trace

v14's 147 pointed at an absent object. v15 ships the build: `evidence/p256_fe_v15.c`
(v10m field kernels verbatim + counted wrappers), `evidence/insn_disasm.txt`
(objdump), `evidence/insn_nm.txt`, `evidence/trace_counts.json`.
`fe_mul_obj` spans **148** instructions (0x0..0x250, 4-byte arm64), `fe_add_obj` and
`fe_sub_obj` **27** each; the affine-addition row is 3*148 + 6*27 = **606** Ir.
Dynamic trace: one affine addition executes **388** `fe_mul` with its own inversion
(3 + 385 chain = 256 squarings + 128 multiplies + 1 fixup) and **3** when the inverse
is supplied by the production batch. v10's 161 lesson applied: the executed number
is shipped, so the per-visit total rises from 1100 to 1220.054 Ir.

## 6. Ledger regenerated from the 48-byte config

Rows: 606 + 444 (3 muls batch share) + 0.054 (385-mul chain / 2^20) + 10 + 20 + 12
+ 128 = **1220.054 Ir** = 1.220054 units at C=1000. Work = 127.67 + log2(1.220054)
= 127.9569, plus restart 6.8e-13 and replay 2^-19, ceilinged to **127.9570**.
Storage is recomputed from 48 B: log2 = 127.67 + log2(48) = **133.2549625**, so the
v14 row's 133.25 was a rounded-down lie and the honest bound (133.2550) still does
not bind: **claim = min(bounds) = 127.9570**. The in-package assertion
`claim >= log2(total)` holds, and the claim ships as the reduced rational
**127957/1000** (gcd 10).

## Success semantics (tracks/p256-ecdlp/rules.md, quoted)

"Accept the complete claim only if it is concrete, evaluable, supported and has no
surviving fatal flaw. Supported plausible heuristics may qualify." … "Algorithmic
success must be at least 0.39 over uniform full-production challenges and fresh
algorithm coins." … "An analytic argument need not run an exponentially large
attack. Small experiments are evidence only at their actual scope, never a
production complexity proof."

v15 therefore states the success law as the birthday heuristic it is: lambda =
V^2/r = 0.6329 with V=2^127.67 and r~2^256, success 0.4689 >= 0.39. The toy census
is disclosed at its own scope and does **not** reach that regime: the toy walk
recurs in about sqrt(1000839)/7 ~ 140 steps, so a lambda=0.63 toy store is already
cycle-dominated (measured 0.297, CP-95 0.260). That number is reported, not hidden,
and it is not used as the production value. Record capacity is demonstrated by
replay, not by identity serialization: 38400/38400 reconstructions equal truth and
the replayed coefficient pair recovers the secret.

## Negative results

* The generalized scan's chance family (680583 hits) does not bind at production:
  its shape probability is 2^-256, not R^-2.
* The toy census cannot validate the birthday law at lambda=0.63 at any feasible
  toy size; the recurrence scale and the birthday scale are both sqrt(group size).
* Selector transfer remains unproven in both v14 and v15.

## Files and commands

`python3 experiments/insn_price.py`, `fullstore_census.py`, `cyclescan.py`,
`record_recovery.py`, `selectors.py`, `ledger_v15.py`, `literal_sweep.py`,
`audit.py`. Evidence JSONs are written to `evidence/`. The literal sweep requires
every integer of 5+ digits in this package to reconcile to an evidence value.
