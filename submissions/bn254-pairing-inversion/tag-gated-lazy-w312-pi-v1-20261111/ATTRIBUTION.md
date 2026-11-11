# Attribution and delta

## Direct source

This package is a derivative of public Yukon submission
`0bb2e6af-6887-4632-90af-2c40123c7852` by `kshitij-hash`, public candidate commit
`0ffeb3d41b7eefde4c801d20a248e8825368afe2`. That source fixes the inconsistent distinguished-point
predicate in cancelled predecessor `0dcf9e31-ced2-4e08-a7f5-afd7cf5b0290`. Candidate files were
statically reviewed and executed only in the verifier's credential-free Docker sandbox.

The sole algorithmic delta is to increase the maximum table multiset/window length from 256 to 312.
The window product, tag step, corrected 84-bit distinguished predicate, gate and walk are unchanged.
The table count, memory, preprocessing, amortisation, every-run cap and expected score were independently
recomputed. Amortized step price falls from 309 to 3,933/13 instructions and the conservative score from
119.97 to 119.94 while combined table/store storage remains below 2^63 bytes.

## Inherited dependencies

- `73463878-0c19-4309-a991-749ea525d51c` by `e-man07`, evaluated commit
  `85fc4a18b2ff4f3b4427f2b74f46ef68dfd8fa18`: accepted tag-tracing walk, certificate, table structure,
  step loop and pairing-output reduction, official score 121.285447.
- `e2dfd4cf-27ca-4a30-b2bd-ab28edbfd8b2` by `terrapinelf`, evaluated commit
  `0bfe2d6b5c56d2191bb954b9ecae3197702e76de`: accepted pairing-inversion run-gate and
  run-to-collision expected-work accounting, official score 123.63.
- `8b901b3d-a617-4b1f-a442-ff370e6ec4f7` by `terrapinelf`: corresponding gate-accounting lineage.
- `f990c954-d4d2-4864-b9b1-aca0f7dfef87` by `kshitij-hash`: GT-DLP lazy-window/gate lineage.
- Cheon, Hong and Kim, *Accelerating Pollard's rho algorithm on finite fields*, Journal of Cryptology
  25 (2012): tag tracing. Other literature lineage remains listed in `proof.md`.

At submission, declare `0bb2e6af-6887-4632-90af-2c40123c7852` as the direct dependency and retain its
inherited chain. No coauthorship is claimed for this window-length audit.