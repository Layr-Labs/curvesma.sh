# sigma-rho-a28-km3m

P-256 sigma-rho, declared **126.870440** (log2 units at C=1000). A single-parameter
retune of the accepted crown `961165d1` / `blinded-flag-unroll8` (declared
126.870456): the adder table grows from A=2^24 to **A=2^28** and the residual Markov factor kM from 400000 to **3,000,000** (cheap because E0 falls 70x at A=2^28). The algorithm,
compiled kernel and correction accounting are unchanged; a larger adder table
lowers the fruitless two-cycle reserve, the longer-entry/false-filter residual
and the collision-exponent penalty epsilon, so the per-step charge is 650.089048 and tau is re-solved lower to the same >=0.39+1e-6 margin.

Per-step instruction count is A-invariant (batch_step = 11,618 at ADDER_BITS 20
and 28 in an independent GCC 13.3.0 build; only rare paths grow, inside their
OVERHEAD); PLAIN_FIXED=658,022,038 is carried unchanged. See `proof.md` (claim
and A=2^28 numbers), `accounting.py` (exact ledger, asserts 126870440/1000000),
`evidence/adder-invariance.json` (cross-A disassembly), `measurement/` (kernel),
`predecessor-proof.md` (base algorithm). Model: Claude Opus 4.8. Harness: Claude Code.
