# stnfs-aslice-v14

Exploratory STNFS claim for bls12-381-fixed-g2-invert-g1: W <= 2^123 reference
units (dense Fp12 multiplications, C = 20,000), success modelled >= 0.599
(required 0.39).

What is new relative to v13 (all parameters one consistent set across code,
prose and manifest):

* **Window-sieved descent** (DESCENT.md, descent14.py): every node including the
  boot is a bucket-sieved special-q lattice node with a certified envelope
  screen, a bounded child window and residual ECM. The v13 exact-norm-per-trial
  screen (honestly priced at two-log 2^126.37 units) is replaced by an algorithm:
  congruence functionals, hit-lattice walks, thresholds and survivors, validated
  by exhaustive toy differentials (stats14.json; sampled in-sandbox).
  Two-log recurrence 2^110.42 units (measured walk anchor; the production leaves/hit
extrapolation is disclosed with both-end sensitivities); H5 x128 < 2^117.42;
metered cap 2^119; modelled success 0.599.
* **Lean record routing** (STORAGE.md, store14.py): per-line bounded buffers,
  fixed pages, two-level directories: 20 instructions per record + 8 per cell
  (apply/recovery charged separately), replacing 256 + 8.
* **Block selection** at certified bound cut 1285: 0.808 of cells, 0.962 of
  yield, measured over all pool rows (stats14.json).
* **Linear algebra**: 16 blocks, four hops, 88 instructions/entry.
* **One ledger** (model14.py, recomputed by evidence14.py): relcol 2^122.183,
  LA 2^120.869, descent cap 2^119, others 2^90/2^88/2^80/2^72/2^11:
  total 2^122.726, claimed 123 with 0.274 bits reserve.

Run the declared experiment with the shared verifier's credential-free runner:
`evidence14.py` with `{"samples": 30, "seed": 1}` (about 10 s in the pinned
sandbox). Nothing here may run on a host or with credentials.
