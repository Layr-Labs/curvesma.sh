# Argument (v21)

1. *Reduction* (unchanged). `g = e(P0,Q)`, `z = g^x`; walk `y = g^A z^D` in GT; a class
   collision with distinct labels gives one linear equation in `x`:
   `x = (p^k A1 - A2)/(D2 - p^k D1)`.
2. *Class function* (unchanged). `kappa` = index of the minimal cyclic rotation of the
   twelve keys; `s = idx(H(m0))`; multiplier `M_s^(p^-kappa)`.
3. *Restarted workers*. Trail ends at its first report; one 52 B record; restart
   from `SHAKE256(sigma_w || j)`; recovery replays both trails from their seeds. Records
   = S/2^100 = 56628313 and records x mean trail = S = 71784914964362045983757977582135410688: the priced step count IS
   the executed step count. Restart row 1679295746777726. `evidence/ledger22.json`.
4. *Abelian trap, exact.* GT abelian => a slot word's product is 1 identically iff every
   per-slot exponent sum vanishes in `Z[x]/(x^4-x^2+1)`. Exact enumeration over all
   canonical words (`scripts/census22.py`): v20's (R=4096, d=1) traps every ~2^33.882
   steps; depth 4 at R=4096 still leaves 2^41.046; the
   shipped (R=2^16, depth 5) leaves 2^1.108, priced as 1066502469749131192112667501083245841 instructions
   of wasted trail budget. `evidence/patterns22.json`.
5. *Multi-target census.* 64 independent (target, table) instances x
   12 trials; pooled 763/768, CP-95
   `0.98636023754575717248371802270412445068359375` by the v15k exact helper; per-instance distribution reported.
   `evidence/census22-toy.json`.
6. *Reserve derived.* retention = CP-95 pooled; `S` solved backwards from the fixed 0.3901
   target, so floor = `0.390100` >= 0.39 with no assigned constant.
7. *Work.* Claim 126.0415, cap `floor(2^126.0415)` covering the nominal work; exact
   rational W = `350211906437225172960629602534771332718927/4000`.
