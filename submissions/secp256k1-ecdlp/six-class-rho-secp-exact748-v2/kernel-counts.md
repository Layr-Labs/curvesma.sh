# Counted AArch64 kernel counts used by the 748-instruction step

The ordinary step's 634 counted instructions are the two generated batch
loops of the public six-class rho construction. This package does not copy
the assembly; it pins the public files by hash and summarizes the counts.
The full assembly, generator, C specification, native test and interpreter
are in the public `dddcbeb3-bc7b-4421-b946-70e15cb833d2` v9 package and its
`bd8e6be5-913b-47a8-b097-526f4bf4be11` revision; the same files are also
present in our earlier public `a6ae92da-3a0c-4493-84c1-bcfd6da47232`
kappa101 package. The present package's contribution is the 114 analytic
rows, not a new kernel count.

## Files and hashes

| public file | SHA-256 |
| --- | --- |
| `kernels/secp_loops.s` | `786565c57af8a38458441c1e40521c02cfb18664ad9e65946a4657dfae697bc2` |
| `kernels/secp_fe.s` | `912dc2aa34060ba72109c5a2bb165347506473a9a376d2f8fd9e72f385c6d636` |
| `kernels/secp_walk.c` | `301bc7987a80ff9f272bdf3faa71a80d5567fce86e40ab096ea95e3d8ba56b90` |
| `kernels/gen_loops.py` | `38f1b1ddffc5f2e812bb50d99c1d9e1c349af3ab58a1721f1a9ae22cea2468fa` |
| `kernels/test_native.py` | `4f22039a7aea57e2dc4765d10f7bd6547e548b42e34e4ac999e1f12b62bf2b5d` |

## Per-iteration counts

| loop | instructions per walk | composition |
| --- | ---: | --- |
| pass 1 | **111** | 1 fe_mul (91), 1 fe_sub (9), 5 loads, 4 stores, 2 loop control |
| pass 2 | **523** | 4 fe_mul (364), 1 fe_sqr (71), 5 fe_sub (45), 2 fe_canon (16), 19 loads, 6 stores, 2 loop control |
| **counted total** | **634** | 111 + 523 |

Both loops are single basic blocks: every conditional is a CSEL, and the
only branch is the loop back edge. The iteration count is therefore
data-independent. The public v9 evidence parses these counts from the
shipped assembly, executes the shipped file in an AArch64 interpreter on
production secp256k1 points, and checks the per-extra-walk difference is
exactly 111 and 523. `kernels/test_native.py` checks the same loops
natively on an AArch64 host against Python big integers and affine
secp256k1 arithmetic. The field-operation counts are 91 for fe_mul
(63 product + 28 reduction), 71 for fe_sqr, 9 for fe_sub and 8 for
fe_canon; the per-iteration identities are
`111 = 91 + 9 + 5 + 4 + 2` and
`523 = 4*91 + 71 + 5*9 + 2*8 + 19 + 6 + 2`.

## Opcode histograms (public v9 evidence)

Pass 1: ADCS 27, MUL 21, UMULH 21, ADDS 11, ADC 10, SBCS 5, LDP 4, STP 4,
SUBS 3, CSEL 2, B.NE 1, LDR 1, SBC 1.

Pass 2: ADCS 132, MUL 99, UMULH 99, ADDS 54, ADC 47, SBCS 25, CSEL 18,
LDP 18, SUBS 11, STP 6, EXTR 5, SBC 5, B.NE 1, LDR 1, LSL 1, LSR 1.

There are no register spills in either loop. The only stack traffic in
pass 2 is the 32-byte slot holding R, already included in the 19 loads and
6 stores. The inversion is outside the per-walk loop and is charged as
per-batch fixed work, below one instruction per walk at the 2^15 minimum
batch; the step's 114-row table includes one instruction for that share.

## How this package uses the counts

`proof.md` section 5 uses `634 = 111 + 523` and the 114 explicit analytic
rows, for a step of 748 instructions = 0.748 units. `experiment.py` in this
package rechecks the 114-row sum, the 748 total, the ledger and the success
floor. It does not re-execute the counted kernels; the kernel evidence is
the public assembly and its public tests, pinned above. A reviewer who
wants to verify the counted 634 should check the public v9 files against
the hashes above and the per-iteration identities.
