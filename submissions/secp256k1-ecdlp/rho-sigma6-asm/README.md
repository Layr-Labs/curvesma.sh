# secp256k1-ecdlp: parallel Pollard rho on order-6 classes (sigma walk, assembly field arithmetic)

Parallel distinguished-point Pollard rho over the order-6 classes given by negation and
the GLV endomorphism of secp256k1. Steps are batched affine additions on 2^32 walks
with 2^16 adders.

The walk stores actual points and adds the image sigma(R_i) of each adder, where sigma
maps the canonical point to the current one. The class sequence is that of the usual
canonical walk, but no step multiplies by beta. Which rotation is canonical is decided
from 64-bit fixed-point fractions, with an exact fallback.

Walks run in trails from fresh coins. Each trail ends at its first distinguished point and
records it once with its coins, and the first key match is resolved by replaying the two
trails. A degenerate match therefore has probability exactly 1/r, and the replay is
bounded and charged. The distinguished-point test precedes every escape, so the 64-bit
filters cannot change which point a trail records. The toy experiment runs this whole
procedure, including resolution, with filters cut to 8 bits and every trail replayed.
Per-step costs are averages; a hard instruction cutoff holds every run to the stated total.

No reference operation is claimed. The whole step, including all field arithmetic, is
the package's own 4x64-bit GCC inline AArch64 assembly. It is charged instruction by
instruction at 1/C with C = 1000: 795 instructions (0.795 units) on average, against 778.7 measured.
The same walk with the 6M reference operation and libsecp256k1 field code costs 1.32
units per step (an earlier version, not included).

The submitted bound is **time_log2 = 126.418986**, with success ≥ 0.39 under the
disclosed random-mapping heuristic.

- [Algorithm, success argument and cost](proof.md)
- [Exact accounting script](accounting.py) and its output in [evidence/accounting.json](evidence/accounting.json)
- [Per-step instruction measurement](measurement/) (Callgrind, organizer toolchain, field code in `secp_field.h` and `secp_sub.h`), results in [evidence/measurement/](evidence/measurement/)
- [Declared toy-curve experiment](experiment.py); offline runs in [evidence/trails/](evidence/trails/)
- [Claim summary](claim.json)

This applies a known generic square-root attack and accounts for it carefully. It is
not a practical break.
