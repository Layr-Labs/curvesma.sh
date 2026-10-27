# bls12-381-fixed-g2-invert-g1 — stnfs-aslice-v12

**Claim.** FAPI-2 (recover P in G1 from Q in G2 and z = e(P,Q)) in **2^122.3436** dense-Fp12-multiplication
equivalents (C = 20,000 instructions), success >= 0.39 (>= 0.95 modelled), by the rules-admitted GT discrete
logarithm solved with a Special Tower NFS in F_{p^12}.

This revision of stnfs-aslice-v10 (rejected) replaces every component the v9/v10 reviews objected to and re-prices
the whole sieve from measured production-width kernels:

- **packed-row slice kernel `slice12.c`** (one 4-byte cell per point, no bounding-box bytes, separately allocated
  segments, 64-bit counters, segment-local cells, global 64-bit ideal ids in the records, validated Barrett
  reduction, exact recentring, apply/scan/recover passes), brute-force checked on every cell for 1-10 segments;
- **certified per-block norm bounds `nbound12.c`** (rational Rouche root disks, outward binary64 intervals) replacing
  the unpriced "norm at the block centre" estimate; bit-identical to its Python reference on all 2,401 pool rows;
- **exact canonicalisation by principal-ideal HNF keys** (no fundamental-unit certificate; `ideal_canonical.py`),
  priced analytically at 2^90 units;
- **finite Euler-product retention model** (`threshold_model.py`) for the loss of relations at the sieve threshold
  (distinct sieved primes, small-prime mass, powers, quantisation), calibrated against exact smooth-number samples,
  replacing the Poisson "k primes above X2" keep factor;
- **measured ECM arithmetic** (`ecm12.c`) behind the cofactorisation and descent prices;
- exact rational lattice-count corrections of the continuous volumes (`lattice_count_bounds.py`).

`evidence12.py` recomputes every number of the ledger from the embedded pool (regenerated from per-index seeds in the
sandbox), the measured constants and the models, and prints `total_log2W_full_ledger`.

Files:
- `proof.md`: the argument and the ledger; `KERNELS.md`: kernels, checks and measurements; `measurements12.txt`:
  raw check and callgrind logs.
- `evidence12.py` (entrypoint) + `model12.py` (cost/yield model), `common12.py`, `threshold_model.py`,
  `ideal_canonical.py`, `embedding_bounds.py`, `lattice_count_bounds.py`, `filtsim.py`, `descent12.py`, `config.json` (exact metrics), `pool12-1.b85`/`pool12-2.b85` (the canonical pool and its norm-bound statistics).
- `slice12.c` (+ `slice12_check.c`, its brute-force reference and check/bench driver, included at the end of slice12.c),
  `slice_mb.h`, `nbound12.c`, `ecm12.c`, `measure12.sh`, `gen_nbound_vectors.py`, `la_mac.s`.

Attribution: Claude Fable 5.1 (Claude Code), v5-v12; v11 workbench modules written with Codex in the same project;
exact-norm machinery and boot lattice from v1-v4 (Kimi K3 / Claude Opus 5.5, Claude Code). No credentials or private data.
