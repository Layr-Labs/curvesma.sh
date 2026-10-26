# hashring-rho-w16ub — sixteen-slot window with upper-bound collision terms

Parent: `6423c2728dc03d5f2d647c7f22adf1dac0643a66`. Source: accepted
`hashring-rho-sel1exact` (`aba345aa`, 125.82) via `b735ab01`.
Claim: `time_log2 = 125.817` (125817/1000), success at least 0.39.

One algorithmic change: look-back window 24 to 16 class keys. Membership row
92 to 68 instructions (16 unrolled LDR/CMP/B.EQ = 48 plus the 20-instruction
hit path), charge 1,749 to 1,725 (9-instruction margin kept over the 1,716
maximum), `log2 W` 125.81838 to 125.81679, claimed at thousandths.

Review-driven correction (answers verdict on `60d50884`): every collision
term is stated as a conservative upper bound, never an exact probability.
Conditional on a window holding `d` distinct 64-bit hash words (`d <= 16`),
a fresh non-member candidate with uniform hash matches with probability
`d*2^-64 <= 16*2^-64 = 2^-60`; immediately after seeding all 16 slots hold
one key so `d = 1` and the probability is `2^-64`. Per accepted step the
collision-rejection probability is at most the union bound over at most 8
probes, `<= 8*2^-60 = 2^-57` as arithmetic on the bound; probing stops at
the first acceptance and collision events may overlap, and for any `q`,
`1-(1-q)^8 <= 8q`, so the union bound holds regardless of dependence or
stopping. The slot update keeps its priced ADD/CMP/CSEL sequence with no
claim about masking. All rounded bounds face outward; the only exact
comparisons are the integer claim checks `W^1000 <= 2^125817` and
`W^1000 > 2^125816`.
