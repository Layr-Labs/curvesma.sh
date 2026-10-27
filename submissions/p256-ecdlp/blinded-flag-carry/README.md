# blinded-flag-carry

P-256 sigma rho, declared **126.879312** (log2 units at C=1000). Pure
cost-ledger repricing of kaankolcu's 126.879316 flag-carry crown
(4ff4847): three allowances move to their listed or measured support
with no shared fate. Blinding-product term 2^11 -> 256 per batch round
(two listed 115-instruction branch-free k_mul rows plus the rho reads).
Rare-event and look-ahead overhead 1,024 -> 256 per step (forced-path
own control at most 191 per check). Counting-method allowance 1,024 ->
768 per round (28x the measured 27-per-round gap). Same kernel, counts,
meter, CORR, premises and success model; DELTA_CORR bound 5.07e-14 and
the 1.0390e-6 margin hold verbatim. C_COMMON 631.620007 -> 631.618053,
per-step total 654.098811 -> 654.096780 (rare reserve 0.005494 ->
0.005418 inside it). The flag-carry base is mpjunior92's 00db6aff
(blinded-corrections, 126.89249; extends aab57cd4) with the six product
finishes branching on the carry flag, 5*(2^20-1)+(2^20-2) fewer
instructions per 2^20-walk round, measured exactly. Fixed path 631.04 +
metered corrections 22.47 per walk step (base 637.04).
Files: `proof.md`, `accounting.py`, `measurement/` (gen8.py, headers),
`evidence/measurement/branchy.json` (`flag_carry`).
