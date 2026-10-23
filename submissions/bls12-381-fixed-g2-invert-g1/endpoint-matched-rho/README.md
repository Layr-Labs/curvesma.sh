# endpoint-matched-rho

Candidate for `bls12-381-fixed-g2-invert-g1`: given nonidentity `Q` in G2 and `z`
in GT, return `P` in G1 with `e(P,Q) = z`.

Mechanism: `g = e(P0,Q)` generates GT, so `z = g^a`; the attack runs `k = 2^20`
fixed-length walks of `L = 240218821685709510680521902710375` steps in GT with
the iteration `x -> x * g^{e(x)}`, `e(x) = 1 + (c11 mod 32)`, started at
`g^{u_i} z^{v_i}` for uniform `u_i, v_i`. Collisions are found by sorting the
`k` final elements: equal endpoints mean the two walks merged, and the merge
gives `a = (u_i + c_i - u_j - c_j) * (v_j - v_i)^{-1} mod r`. The returned
`P = [a]P0` is checked with one pairing.

- total work: 2^127.602539 target-group-operation equivalents (dense Fp12
  multiplication, `C = 20000`), with `k*L = 2^127.566048` multiplications fixed
  before the instance is seen;
- success: 0.45393 under the standard random-mapping model, declared as 0.39;
- peak memory: 2^31 bytes; preprocessing 2^29.581231 equivalents, included;
- advice: none.

`proof.md` holds the algorithm, the correctness argument, the success analysis,
the complete ledger and the disclosed heuristics. `claim.json` binds the
protected `target-profile.json` and `classical-group-work-v2.json` by SHA-256.
