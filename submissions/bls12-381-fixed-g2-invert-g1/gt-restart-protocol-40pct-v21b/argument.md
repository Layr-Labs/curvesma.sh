# Argument (v21)

1. *Reduction* (unchanged). `g = e(P0,Q)`, `z = g^x`; walk `y = g^A z^D` in GT; a class
   collision with distinct labels gives one linear equation in `x`:
   `x = (p^k A1 - A2)/(D2 - p^k D1)`.
2. *Class function* (unchanged). `kappa` = index of the minimal cyclic rotation of the
   twelve keys; `s = idx(H(m0))`; multiplier `M_s^(p^-kappa)`.
3. *Restarted workers*. Trail ends at its first report; one 52 B record; restart
   from `SHAKE256(sigma_w || j)`; recovery replays both trails from their seeds. Records
   = S/2^100 = 56443886 and records x mean trail = S = 71551125967113754319928145100077531136: the priced step count IS
   the executed step count. Restart row 1226117666371992. `evidence/ledger21.json`.
4. *Abelian trap, exact.* GT abelian => a slot word's product is 1 identically iff every
   per-slot exponent sum vanishes in `Z[x]/(x^4-x^2+1)`. Exact enumeration over all
   canonical words (`scripts/census21.py`): v20's (R=4096, d=1) traps every ~2^33.882
   steps; depth 4 at R=4096 still leaves 2^41.042; the
   shipped (R=2^16, depth 5) leaves 2^1.103, priced as 896159476432169845606543383971105498 instructions
   of wasted trail budget. `evidence/patterns21.json`.
5. *Multi-target census.* 64 independent (target, table) instances x
   12 trials; pooled 763/768, CP-95
   `0.98636023754575717248371802270412445068359375` by the v15k exact helper; per-instance distribution reported.
   `evidence/census21-toy.json`.
6. *Reserve derived.* retention = CP-95 pooled; `S` solved backwards from the fixed 0.3901
   target, so floor = `0.390100` >= 0.39 with no assigned constant.
7. *Work.* Claim 125.7905, cap `floor(2^125.7905)` covering the nominal work; exact
   rational W = `735689151341568781095733404182287909428637/10000`.
