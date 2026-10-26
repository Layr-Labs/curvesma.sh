# Fixed-G1 pairing inversion: Frobenius-class rho with a cyclotomic-squaring branch, priced by an itemised and measured instruction count

Claim: time_log2 = **125.62** (total work <= 2^125.62 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^74 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)
- [Squaring price: itemisation, kernel source, bit-exactness check and Callgrind record](measured-kernel.md)

Track `bls12-381-fixed-g1-invert-g2`: for the supplied nonidentity P in G1 and uniform z in GT, compute
g = e(P, Q0) (pinned order pairing(Q0, P)), return Q = O when z = 1, otherwise solve z = g^x with the walk
below and return Q = [x]Q0 in G2, verified by e(P, Q) = z.

## Port to this track, and what this revision changes

The package starts from the fixed-G2, recover-G1 package `gt-tracecounted-exact-v25` (submission
`56fc259f`, 125.72) and carries its walk, keys, escapes, caps, heuristics and evidence to this target.
The port changes the reduction (proof.md section 1: g = e(P, Q0) per instance, output Q = [x]Q0, the
identity target z = 1 returning the point at infinity without a walk) and the setup and output lines of
the ledger (section 5). Those are counted in `experiment.py` part F: the pinned pairing is transcribed
and every field operation of it charged one full unit (10,331 units per pairing, two pairings), and the
G2 output is charged at its hard maximum for a 255-bit scalar (254 doublings and 254 additions over Fp2
plus the affine encoding, at most 12,809 units). Part D adds those counts to the setup bucket (at most
2^29.9973 units against the 2^31 it is allowed).

A further revision replaces the inherited chained multimap for distinguished-point reports with one
bounded store: reports are appended to a flat buffer of 80-byte records, sorted by the 128-bit
fingerprint with eight fixed stable 16-bit radix passes, and scanned for equal neighbours. The
per-report cost is itemised at 3,831 instructions and charged at 4,000 (12 keys 2,119; two 255-bit
modular multiplications 1,024; append 32; eight passes 640; scan 16), so the DP-handling term is at
most 2^63.80 units and the non-walk bucket stays below 2^64.5. The fixed pass count and the per-pass
cost do not depend on the data; no report is dropped and no cost rests on a load factor. The walk
ledger and the claim are unchanged. The report store is exercised by part G of `experiment.py`, and
part D prints the itemised charge.

This revision then re-prices the cyclotomic-squaring step for a straight-line wide-lazy kernel
written over blst's exported arm64 primitives: **9,709 instructions per call (0.48545 unit)** instead
of the parent's audited 11,369-instruction narrow listing. The itemisation multiplies the kernel's
fixed operation counts by the static instruction length of each pinned primitive body (9,319) and
adds the disassembled 379-instruction compiled glue body and the 10-instruction loop: 9,708.
Callgrind on the same compiled kernel reads **9,708.0022 instructions per call** (79,527,954 `Ir`
over 8,192 calls), and the charge is the ceiling of the larger figure. `measured-kernel.md` lists
every count and length and reproduces the kernel source, the bit-exactness check and the Callgrind
harness and record. The kernel is bit-exact against `blst_fp12_cyclotomic_sqr` on 1,000 random
cyclotomic elements, and the same run's control measures 9,803.0018, so the kernel saves 95.0
instructions by flattening the library's call glue. Because a
squaring step is cheaper relative to a table multiply, the
squaring probability is re-tuned from 3/8 to 35/64, the model optimum at that price. The branch,
multiplier, distinguished-point and fingerprint bits are disjoint slices of the class invariant
T = Re(c_0) (T mod 64, bits 6-15, bits 16-75 and bits 76-203), so the branch mix and the multiplier index
stay independent under the package's H2. The walk ledger becomes W = 2^125.6197624 and the claim 125.62.

## Lineage of the source revision

The machine, algebra, walk, branch, escapes, caps, heuristics and evidence are those of the live
gate package `frobenius-rho-trace-branch` (jungjipdo, submission `747ad6c9`, 125.74), a revision of
`6d4e57d1` (125.80) and of the Frobenius-class rho of `7da5f23a` (mpjunior92). That package
itemises three per-step charges and then charges each above its own counted listing: 300 for 270
instructions on every step, 2,200 for 2,119 on a multiply step, 11,400 for 11,369 on a squaring
step. The first port charged the listing itself for all three: multiply step 1.125 -> 1.11945 units,
squaring step 0.585 -> 0.58195, average step cost c 0.9225 -> 0.9178875, W 2^125.7147 -> 2^125.7075.
This revision keeps the 270 and 2,119 charges, replaces the squaring listing by the itemised and measured
kernel charged at 9,709 instructions (0.49895 unit including the 270 of bookkeeping) and
re-tunes the branch
probability: average step cost c = 0.7801140625, W = 2^125.6197624, claim **125.62**.

The step price leaves a much smaller cushion under the claim: the abort cap 2^125.62 exceeds
the ledger by 0.0164% (0.0002 bit), so a 0.0207-point drop of the squaring fraction already trips it,
while in seeds 1-5 the largest deviation was 0.70 points (sampling sd about 0.64 points at 6,000 steps). A tripped cap only shortens
the walk within the bound; proof.md section 3.4 gives the success value at each drift (0.428 at 0.5
point, 0.420 at 2 points) and shows the 0.01 success budget covers drifts up to about 2 points. Section 0a states the limb-level rule (branch bits
and key minima read from the stored Montgomery representation, which is what the priced kernels
produce), and part C runs both that rule and the ordinary-integer rule as evidence.
`experiment.py` part D reproduces the parent's ledger (2^125.7147 at c = 0.9225 and s = 3/8) and this
one through the same code path.

This is an analytic generic-attack bound, not a practical break.
