# Directed quadrature refinement (inherited, e3940d11, used unchanged)

Mesh 384 (was 48), Dickman spacing 1/1024 (was 1/128); PLAN16/BOOT16 scale every iA by eight so nA is unchanged; TRUNC = 0; retention and descent share the mesh. Eight entries re-optimised at the fine mesh (see source-core.md section 8).

**Validity.** Each old interval splits into eight whose log-measures sum exactly to the old one; feasible subcells stay feasible, sum(t)=m crossings discarded. Left endpoints give lower rho; the upper bracket is consistent. Same model, same heuristic; eight probes matched full rows exactly in binary64. Evidence-only speedup, not a charged kernel.

**Limitations.** rho coefficients in 200-digit Decimal, evaluated in binary64 (numerical evidence, not an interval certificate). Model and extrapolation remain H1-H5.
