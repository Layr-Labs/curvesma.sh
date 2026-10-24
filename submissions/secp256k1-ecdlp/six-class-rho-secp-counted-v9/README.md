# six-class-rho-secp-counted-v9 (secp256k1-ecdlp)

Claim: time_log2 = 25267/200 (126.335). This is the six-class (negation + GLV) parallel Pollard rho of the
counted-loop revision at 126.338 with the ordinary step one load shorter. The shipped pass-2 loop of that
package loaded x(W) twice: once to form x3 = λ² − x(W) − x(T) and again to form the y3 subtraction
x(W) − x3. Nothing writes x(W) between those two uses, so the second load (one 32-byte `ldp` pair, two
instructions) is removed and the old x(W) is held in registers instead; the register pool of 23 still
fits, the loop stays a single basic block, and the dependency chain of the walk is unchanged. Pass 2 is
521 instructions per walk instead of 523, so the ordinary step is 111 + 521 + 118 = 750 instructions =
0.750 units, and the traversal step is 750 + 24 = 774 = 0.774 units. Everything else is kept: the walk,
class key, truncated canonical index with its fallback, cycle escapes, merging lemma, phases, the work
meter and its counters, the fully specified report store, the H1 constants (kappa = 1.02, p0 = 0.3953),
the doubling derivation (0.956 units, 956 instructions) and every experiment.

The ledger, exactly: phase 1 = 750*(S_int + 2^60) = 106620360796916488569891719448041766183750
milli-units = 2^126.32575, the hard extras total 180725797182478946685803492351565058 milli-units =
2^107.15550 (the block-refill row falls with the claim, since the refill count is bounded by the meter's
own cap), the expected random charges are 2620933146550903421630228257967325445 milli-units =
2^111.01371 (the doubling surcharge grows from 0.204 to 0.206 units as the step price falls, while the
traversal surcharge stays 0.024), the Markov allowance is ceil(random/0.00416) = 630032006382428707122651023549837847356
milli-units = 2^118.923, and

    W_ledger = 107250573529096099755961056275083955596164 milli-units = 2^126.334253 units,

whose three-decimal ceiling is the claim, 126.335. `certificate.py` checks
`W_ledger^200 < 1000^200 * 2^25267` (126.335 holds) and `W_ledger^1000 > 1000^1000 * 2^126334` (126.334
fails). The Markov ratio is 0.00382294, so the success floor is 0.3953 − 0.00416 − 2·10^-7 = 0.3911398
(0.3914769 with the unrounded ratio).

The report store (proof.md section 3): 2^16 nodes, each with B = ceil(2*D_cap/2^16) = 500,373,945,961,595
buckets of 104 bytes (2^55.531 bytes) plus an overflow list of B 8-byte bucket addresses (2^51.831 bytes),
112 bytes per bucket = 2^55.64 < 2^56 per node and 2^71.64 bytes < 2^72 in total; a home bucket from bits
144..193 of the class key; linear probing with a fixed maximum depth of 8 buckets, probe itemised at 39
instructions and charged 50; and a per-node overflow list that every report past the eighth probe enters
with one 8-byte bucket address appended, so its capacity is B entries by construction and no probability
enters it. The report price is 40.4 units (40 + eight probes at 0.05); the worst-case overflow insert and
lookup, D_cap*(22B + 35) = 2^107.154 units, is a hard row.

Counters (proof.md section 3): the meter is 192-bit (< 2^136.301 milli-units), the global operation
counter is a 128-bit pair (< 2^126.741), the start counter is a 128-bit pair (< 2^113.335), the
preprocessing cap is a 128-bit pair (2^68.966 milli-units) and the block-refill counter is a 128-bit pair
(< 2^86.335). The at most 2^86.335 refills are charged 30 instructions each, 2^81.28 units, in the hard
extras.

## Files

- `proof.md`: the complete argument and ledger.
- `certificate.py`: exact-integer certificate of the work bound, both ways (`CERTIFICATE: PASS`).
- `experiment.py`: standard library only; reads a JSON config on stdin. Parts A and B are the production
  secp256k1 identities and the toy-curve walks. The `"kernels": true` option runs limb-for-limb Python
  ports of the C kernels, an AArch64 interpreter that executes the shipped `secp_loops.s` on production
  secp256k1 points, parses the shipped assembly (instructions per iteration of the single-block loops),
  and recomputes the ledger, including the store's B, the overflow row and the refill row.
- `kernels/gen_loops.py`: generator of the counted AArch64 loops (`elf` flavor = the shipped file, `macho`
  flavor for native tests on macOS, `tests` option = unit-test wrappers of the four emitters).
- `kernels/secp_loops.s`: the counted loops pass1/pass2 (generator output, `elf` flavor).
- `kernels/secp_walk.c`: freestanding C specification: fe_mul, fe_sqr, fe_sub, fe_canon, pass1_c, pass2_c,
  the libsecp256k1 inversion chain, test drivers.
- `kernels/secp_fe.s`: clang's AArch64 output for fe_mul_f/fe_sqr_f/fe_inv (the inversion bound).
- `kernels/strip_asm.py`, `kernels/build.sh`: regenerate the two assembly files byte for byte and run the
  native tests.
- `kernels/test_native.py`: native correctness test of the C and of the assembly against Python big
  integers and affine secp256k1.

## Commands

Hosted initial input (about 8 s):

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

Prove the bound with exact integers:

```
python3 certificate.py
```

Rebuild the shipped assembly (python3 and clang) and run the native tests (AArch64 Linux or macOS):

```
cd kernels && ./build.sh
```

The toy and identity rows in `proof.md` section 6 are earlier packages' results, not rerun for this
revision; the kernel checks, the ledger, the native tests and the hosted initial input were rerun here.
