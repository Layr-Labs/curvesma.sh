# bn254-gt-dlp: gated segmented distinguished-point Pollard rho

**What this is.** A generic (square-root) attack on the BN254 target-group discrete
logarithm, priced end to end: one dense Fp12 multiplication per walk step, a
bounded distinguished-point table, a fresh-coin run gate, and every one-time cost
charged to this target. The submitted score is `log2 E[W] = 126.888233004` reference
units, with success probability at least 0.39 over uniform challenges and fresh
coins.

**What this is not.** It is not the best known attack on this target and it is not
close to the current research frontier of this track. Published submissions on the
same benchmark reach `2^88.4` with special tower number field sieve (STNFS)
arguments. Nothing here reproduces, checks or improves those claims, and this
submission does not assert optimality of anything. It is a fully specified,
fully charged generic baseline whose arithmetic is executable and whose
collision-time assumption is measured at toy size.

**Scope of the claim.** Target: the protected profile `target-profile.json`
(Ethereum alt_bn128, GT ⊂ Fp12 of prime order r, fixed generator
g = e(P0, Q0)). Input: only z = g^x with x uniform in [0, r). Output: x. The
identity target (z = 1) is included and is answered as x = 0.

---

## 1. Algorithm

All arithmetic is in the profile's tower. The profile defines
Fp12 = Fp[w]/(w^12 - 18w^6 + 82) with `u = w^6 - 9`, so with v = w^2:

    Fp2  = Fp[u]/(u^2 + 1)
    Fp6  = Fp2[v]/(v^3 - xi),        xi = 9 + u
    Fp12 = Fp6[w]/(w^2 - v)

and the w-basis coefficients c_0..c_11 map to the tower by

    a_m = (c_m + 9 c_{m+6}) + c_{m+6} u,      m = 0..5,

where m = j + 2k selects the Fp6 half j and the v^k slot. `evidence.py`
converts both ways and cross-checks the tower product against an independent
schoolbook multiplication in the w-basis polynomial ring on deterministic test
elements.

### 1.1 Walk

Fix a 32-entry table. For j = 0..31 let

    (c_j, d_j) = (H("bn254gt-c" || j) mod r, H("bn254gt-d" || j) mod r)
    M_j       = g^{c_j} · z^{d_j}

with H the integer reduction of the stated byte string (a public, deterministic
pseudorandom function; no secret advice). A state is (Y, a, b) with the invariant
Y = g^a z^b. One step is

    j  <- top five bits of the first Fp coefficient of Y
    Y  <- Y · M_j            (one dense Fp12 multiplication)
    a  <- a + c_j mod r
    b  <- b + d_j mod r

The partition is a function of Y only, so the walk is deterministic given its
start; it is the standard 32-adding walk (Teske, Math. Comp. 70 (2001) 809-825).

### 1.2 Distinguished points and segments

Choose W = 12,700 table slots and T = ceil(1.6 · floor(sqrt(r))) total steps.
Then the mean segment length is L = T/W = 2^113.843877. Define the distinguished
point test on a step's Y by

    DP(Y)  iff  the 128-bit integer t(Y) formed by the low two words of the
                imaginary part of the highest Fp2 coefficient of Y satisfies
                t(Y) < D,   D = 18256   (probability 2^-113.843917 per step).

A *segment* starts at a fresh uniformly random state
Y_0 = g^{a_0} z^{b_0} (a_0, b_0 uniform mod r) and iterates until the first
distinguished point. The pair (128-bit fingerprint of that DP, a, b) is stored in
the table; if the fingerprint already occurs, the attack stops and recovers x.
At most W segments are run; the loop also stops at the step budget T. The table
holds 80 bytes per entry (16-byte fingerprint plus two 32-byte scalar labels), so
peak storage is 12,700·80 + 32·448 + 4,096 = **1,034,432 bytes < 2^20**.

Recovery on a fingerprint match. Let (a, b) be the current labels and (a', b')
the stored labels of a segment ending at the same Y:

    g^a z^b = g^{a'} z^{b'}  =>  a + b x = a' + b' x (mod r)
    x = (a - a') / (b' - b)  mod r,      defined when b' != b.

The returned x is then checked by a full exponentiation g^x == z, comparing all
twelve coefficients; a failed check is a failed trial and never an accepted
answer. If b' = b the trial fails (probability exactly 1/r, see §3). The identity
input is answered as x = 0 immediately.

`evidence.py` implements exactly this walk (same partition, same table, same
recovery) on the real BN254 GT and on prime-order toy groups.

---

## 2. Cost

The unit is the protected one: R/20,000 for R scalar AArch64 instruction
equivalents. Costs are per target; nothing is amortised over other targets.

| item | charge |
|---|---|
| one walk step: one dense Fp12 multiplication | 22,996 |
| partition (5-bit extract, table index) | 24 |
| two 4-limb label updates mod r (load, ADDS/ADC, SUBS/SBC, select) | 128 |
| distinguished-point test (2-word compare against D) | 8 |
| loop counter, branch, record access margin | 40 |
| **per step** | **23,196** |
| table build: 64 GT exponentiations + 32 records | 2^30 |
| DP handling: ≤ 12,700 fingerprints, serialisation, table ops | 2^32 |
| recovery, inversion mod r, g^x == z verification, output | 2^22 |
| **one-time total charged** | **2^36** |
| **segments** | 12,700 |
| **steps T** | 236,714,811,010,862,895,392,714,499,903,555,895,296 |

    R_run        = T · 23,196 + 2^36  = 2^141.98 instructions
    W_run        = R_run / 20,000     = 2^127.690294 units
    E[W]         = p_open · W_run + (1 - p_open) · 1 = 2^126.888233004
    p_open       = 0.39 / 0.68 = 0.5735294...
    overall success = p_open · 0.68 = 0.39

The gate is a fresh 64-bit coin: with probability p_open the run is opened and
executes the T-step budget; otherwise it stops at once and costs one unit (the
closed branch does no group arithmetic). Expected work is therefore
p_open·W_run + (1-p_open)·1 and overall success is exactly 0.39.

`machine.md` itemises the instruction counts and records the kernel measurement.
The 22,996 figure is the benchmark's own reference measurement for a **dense Fp12
multiplication** (`tracks/docs/COST_MODEL.md`: blst `blst_fp12_mul`, 22,996
AArch64 instructions, the value C = 20,000 rounds from). It is applied to a
254-bit BN254 field although it was measured on a 381-bit BLS12-381 field, so it
is the conservative side of the analogy: a BN254 multiplication has four 64-bit
limbs per Fp element rather than six. We did **not** implement that assembly
kernel. The portable C kernel shipped here (`kernel.c`) executes the same tower
multiplication and measures **52,313** instructions per multiplication
(callgrind `Ir`, x86-64, gcc 14.2.0 `-O3 -fno-tree-vectorize -fno-builtin`);
charging that measured figure instead would move the claim to
`log2 E[W] = 128.067`, i.e. 1.18 bit worse. The submitted number is the reference
charge, not the C figure.

Memory is reported, not scored. Preprocessing and advice: none beyond the public
profile and the 32 table entries, whose construction is charged above. Public
profile constants, the walk design and the algorithm-selection effort are the only
excluded items.

---

## 3. Correctness

**Invariant.** Y_0 = g^{a_0} z^{b_0} by construction, and each step multiplies Y
by M_j = g^{c_j} z^{d_j} while adding (c_j, d_j) to the labels, so
Y_i = g^{a_i} z^{b_i} for every i by induction on exponent addition. `evidence.py`
runs the real walk on the real GT for 200 steps and re-derives g^{a + b x}
independently every 50 steps; the check is exact (0 bad steps).

**Recovery.** If two states coincide, g^{a} z^{b} = g^{a'} z^{b'}, and with z = g^x
this is a - a' = x(b' - b) mod r. r is prime, so a nonzero b' - b has a unique
inverse and x is the unique discrete logarithm. Every returned x is verified by
g^x == z before output.

**Degeneracy.** Conditional on the true x and on the group elements visited, each
d_j that has not been conditioned remains uniform, so a collision between two
segments whose label difference b' - b is a nonzero combination of independent
uniform d_j is zero with probability exactly 1/r. With r ≈ 2^253.6 this is
negligible and is charged as a possible failed trial.

**Distinct starts.** Starts are drawn uniformly and independently; a repeated
start would be a wasted segment, whose probability is at most W^2/(2r) ≈ 2^-213.

---

## 4. Success probability, and what is not proved

Let T = u·floor(sqrt(r)) with u = 1.6. Under the standard random-map model for a
partitioned r-adding walk, the number of colliding segment pairs is
approximately Poisson with mean T^2/(2r), so

    P(success at T) ≈ 1 - exp(-u^2/2) = 1 - exp(-1.28) = 0.721963.

We charge **0.68**, not 0.721963, and the gate converts 0.68 into an overall
success of exactly 0.39. The discount is deliberate: 0.68 is below the model value
and below the larger of the two measured toy rates.

**Measured support (scope-limited).** `evidence.py part=toys` runs the identical
walk, the identical one-DP-per-segment bookkeeping and the identical recovery on
prime-order subgroups of F_q^* (n ≈ 2^27 and n ≈ 2^33), with a fixed step budget
of 2·sqrt(n) and 256 table slots, recording the step count at first success:

| n | trials | successes | fraction successful by u = 1.2 / 1.4 / 1.6 / 1.8 / 2.0 |
|---|---|---|---|
| ≈ 2^27 | 40 | 32 | 0.500 / 0.575 / 0.675 / 0.750 / 0.800 |
| ≈ 2^33 | 40 | 35 | 0.475 / 0.600 / 0.725 / 0.825 / 0.875 |
| model 1-exp(-u²/2) | — | — | 0.513 / 0.625 / 0.722 / 0.802 / 0.865 |

The measured rates bracket the model value at u = 1.6 (0.675 and 0.725 against
0.722); the smaller group is pessimistic, as expected from finite-size effects.
These runs are evidence at their stated size only. They are not a production
discrete-logarithm execution and they do not prove the collision-time behaviour
of this walk at 2^253.6.

**Disclosed heuristics.**

- H1 (walk statistics). The 32-adding walk with this partition and these
  increments behaves like a random map on the 2^253.6-element group, so a
  collision appears after Θ(sqrt(r)) steps with the birthday constant above.
  SHA-256-derived increments do not prove this.
- H2 (independence across segments). Segment starts and increments are treated as
  independent draws over uniform challenges and fresh coins; conditioning on the
  target x is handled as in §3.
- H3 (gate). The coin is fresh per target and independent of the walk.
- H4 (arithmetic price). One dense Fp12 multiplication costs 22,996 instruction
  equivalents. This is the benchmark's reference measurement for the operation,
  applied conservatively to a smaller field; our own C kernel is 2.27× slower and
  the alternative score is stated in §2.

H1 is the only assumption that the score depends on for success (the work bound T
holds regardless). If H1 failed badly, the run would fail and the claim's success
bound would be wrong, not the work bound.

---

## 5. Comparison, limitations, provenance

- **Not a break and not a record.** 2^126.9 is a generic square-root attack.
  Published STNFS submissions on this same track claim about 2^88.4; the promoted
  baseline on the benchmark board is 2^135.216266. This entry improves on that
  promoted baseline and is dominated by the STNFS line, which we did not
  reproduce.
- **No tail claim.** The 0.39 is an average-case statement over uniform targets
  and fresh coins, as the rules require; individual runs can exceed the charged
  expectation only through the one-unit closed branch.
- **No amortisation.** The table, the 32 multipliers, all segments, all failures
  and the verification are charged once for this one target.
- **No disguised oracle.** The partition uses five bits of a computed element;
  hashing, serialisation and table access are charged.

**Provenance.** The rho method is Pollard (1978); the r-adding walk and its
statistics are Teske (2001); distinguished points and parallel collision search
are van Oorschot and Wiener (1999); the batched cycle-detection alternative is
Brent (1980). The tower, the field constants and the generator are the protected
profile's. No prior submission's code, data or measured constants were used: the
walk, the arithmetic, the kernel and the ledger here are written for this entry,
and the only externally sourced number is the benchmark's published
reference-operation measurement cited in §2 and `machine.md`.

## 6. Reproduction

```
echo '{}'                        | python3 -I -B evidence.py   # checks + ledger + short toys
echo '{"part":"toys","trials":40,"bits":[26,32]}' | python3 -I -B evidence.py
echo '{"part":"rebuild"}'        | python3 -I -B evidence.py   # recompute the declarations
echo '{"part":"kernel"}'         | python3 -I -B evidence.py   # replay the recorded counts
gcc -O3 -fno-tree-vectorize -fno-builtin -o kernel kernel.c && ./kernel g
valgrind --tool=callgrind ./kernel 4096 ; valgrind --tool=callgrind ./kernel 8192
```

The last two commands reproduce `kernel-counts.json` (the difference of the two
`Ir` summaries divided by 4096); `./kernel g` checks g^r = 1 in Fp12 for the
profile's generator using the shipped kernel.
