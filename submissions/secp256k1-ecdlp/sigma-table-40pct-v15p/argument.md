# Argument v15p

1. Layout. `scripts/rec_layout.h` defines `rec_t` once: x[4], y[4], idx, len, esc
   = 88 bytes at offsets 0/32/64/72/80, with `_Static_assert` on sizeof and on
   every `__builtin_offsetof`. The listing `evidence/step_secp_v15g.s` strides
   `#REC_SIZE` (88) and reads/writes len and esc at 72/80, no longer reusing
   offset 64 as a length. `scripts/check_layout.py` compiles the struct TU and
   compares its offsets against the listing's `.equ` constants; the package
   cannot be regenerated if they disagree.

2. Labels. The step is `A' = s*phi^-slot*(A + J)`. `v15o_label` therefore adds
   the selected row's exponent pair `(ju, jv) mod r`, then multiplies both
   exponents by the unit multiplier `s*lambda^-slot` held as a table constant
   (`R`, `lambda^2`, `lambda` for slot 0,1,2) in Montgomery form. Montgomery one
   is `R mod r`; the negation is unconditional `r - m`, valid because `m` is a
   unit; the sign comes from `v15o_align`, the single canonicalisation site.

3. Index. `idx = 3*(mix(key) & 2047) + slot` with rows generated as
   `row[k][u] = phi^-u J_k` and coefficients scaled by `lambda^-u`. The address
   contains the input-dependent slot, and that is exactly what equivariance
   needs: `step(phi^j P)` uses jumps differing from `step(P)` by `phi^j`.

4. Cost. Each row of the composition table is a distinct symbol; the ledger
   recomputes them from the shipped listing and asserts `SUM(rows) == W` as exact
   rationals, with a planted-+1 control.

5. Status. `evidence/label-invariant-v15p.json` is the gate and it FAILS: the
   shipped `fmul` and `rmont` bodies are still incorrect against big-integer
   probes (5/349 and 96/469 mismatches), so the recovery invariant does not hold
   and no correctness claim is made. The instruction counts are static counts of
   the emitted symbols and are unaffected in kind by that defect.

6. Scope. Success is the v15m fixture census (2643/6400, CP-95 ~0.398, toy curve
   q25, dp1, 16 trails, mu = 1.107), a statement about a different program; the
   transfer premise is not asserted.
