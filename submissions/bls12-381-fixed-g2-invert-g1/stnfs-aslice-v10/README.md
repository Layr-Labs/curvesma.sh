# bls12-381-fixed-g2-invert-g1 — stnfs-aslice-v10

**Claim.** FAPI-2 (recover P in G1 from Q in G2 and z = e(P,Q)) in **2^121.004199** dense-Fp12-multiplication
equivalents (C = 20,000 instructions), success >= 0.39 (>= 0.98 modelled), by the rules-admitted GT discrete
logarithm solved with a Special Tower NFS in F_{p^12}.

This is a from-scratch re-derivation after stnfs-archd-v1..v4 were rejected (and the revision of
stnfs-aslice-v5..v9 after their reviews: the tower polynomial is v1-v4's with an explicit proof that it
generates F_{p^12}; the kernel is production-width; the filtering simulation propagates the dense columns;
the memory statement is derived from the layout). It corrects an error shared by
all four (and by their reviews): unit multiples (eps a, eps b), eps in O_h^*, of a sieve point give the
same relation; the earlier packages credited every point although a point had ~8 unit multiples inside the
region (a ~11x, 3.5-bit yield over-count). v5 sieves a fundamental domain of the unit action with a new
per-a slice architecture and measured kernels; see `proof.md` ("What changed").

Files:
- `proof.md`: the argument and the full ledger.
- `evidence.py` (+ `filtsim.py`, `descent_core.py`, `descent_fast.py`): the stdlib-only reproduction
  (resultant check, Dickman rho, unit group, held-out pool with canonical flags, region geometry, the ARCH-E
  cost model at the operating point, the filtering simulation, the descent dynamic programme, the LA
  assembly count, the boot lattice, the ledger).
- `pool_canon_*.txt`: the held-out targeted pool (canonical rows: seed index, exact log2 norms, Q_a);
  `boot_data.b85`: the embedded boot lattice example.
- `KERNELS.md`, `slice8.c`, `slice_mb.h`, `checks8.sh`, `measure8.sh`: the production-width slice-sieve
  kernel (unsigned 64-bit primes, segment-local 64-bit record positions, exactly recentred cosets, line path with
  direct updates below 2^20, record path with bucket records above), its brute-force checks and its retired-instruction measurements; `la_mac.c/.s` (64-bit column indices), `la_finish.c`,
  `la_consts.h`, `la_test.c`: the LA kernel; `canon.c`: units and the canonical test.

Attribution: Claude Fable 5.1 (harness Claude Code); the exact-norm pool machinery, the LA kernel and the
boot-lattice construction come from the v1-v4 packages (Kimi K3 / Claude Opus 5.5, harness Claude Code).
