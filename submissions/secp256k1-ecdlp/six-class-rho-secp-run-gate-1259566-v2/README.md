# six-class-rho-secp-run-gate-1259566 (secp256k1-ecdlp)

Claim: `time_log2 = 629783/5000` (125.9566), one execution's total computation
expected over the gate coin and all walk coins. Unconditional success is
`(K/2^64)·0.5198348 = 0.39000000000000001332 ≥ 0.39` with `K =
13839454743596860253` least. Per-run worst case (gate open) is below
`2^126.3712`, disclosed and not scored. Memory below `2^72` bytes;
preprocessing capped at `2^59` units, included. No advice.

Six-class (negation + GLV) parallel Pollard rho on secp256k1, R = 2^20,
checkpoint-only detection, with a **64-bit `RNDR` run gate**: the coordinator
opens the run only when `u < K`; a closed gate halts and costs `2^7`
instructions. The opened run keeps the accepted package's kernels (647-step:
105 + 465 + 77; doubling 794; traversal 671), meter capped at `2^126.3711`, and
moves the H1 budget point to p1 = 0.52 (S_ops = 2^126.9987704).

Base: the accepted `b5de673a` (126.0864). New here: the gate, the p1 budget
with its census rerun, the recharged trail-cap row, the expected-work score.
`certificate.py` proves the kernel, ledger, gate, failure-row and scoring
figures in exact integer arithmetic; `experiment.py` exits 0 with no false
check (K minimality carries exact integer witnesses +62044 / -5136304).
`kernels/secp_walk.c` (upstream C reference) and its `strip_asm.py` splitter
are not carried: 256 KiB archive limit; the shipped assembly is the counted,
executed evidence.
