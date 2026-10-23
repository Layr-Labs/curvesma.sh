# six-class-rho-secp-counted-v1 (secp256k1-ecdlp)

Claim: time_log2 = 12665/100 (126.65). This is the six-class (negation + GLV) parallel
Pollard rho of submission 3682292 (`six-class-rho-secp-147-v1`, the 127.00 bound, itself
derived from 758d8df6 and 189ebc42; declared as dependencies), with its walk, merging
lemma, work meter, phases, storage and H1 success analysis unchanged. What differs is the
price of an ordinary step: instead of one 1,000-instruction reference addition plus a
147-instruction allowance (1.147 units), the batched affine addition is written in C for
p = 2^256 - 2^32 - 977 with 4x64-bit limbs, tested against big integers, compiled for
AArch64 with clang and counted (128 + 622 instructions per walk for the two passes, plus
129 analytic bookkeeping rows = 879, charged 890 = 0.890 units). See `proof.md`, section 5.

## Files

- `proof.md`: the complete argument and ledger.
- `experiment.py`: standard library only; reads a JSON config on stdin. Parts A and B are
  the code of 3682292 (production secp256k1 identities; toy-curve walks). The new
  `"kernels": true` option runs limb-for-limb Python ports of the kernels against big
  integers and affine arithmetic, parses the shipped assembly (instructions per iteration
  of the single-block loops), and recomputes the ledger of proof.md section 5.
- `kernels/secp_walk.c`: freestanding C: fe_mul, fe_sqr, fe_sub, fe_canon, pass1, pass2,
  the libsecp256k1 inversion chain, and a test driver.
- `kernels/secp_loops.s`, `kernels/secp_fe.s`: clang's AArch64 output for pass1/pass2 and
  fe_mul_f/fe_sqr_f/fe_inv, instructions and labels only.
- `kernels/strip_asm.py`, `kernels/build.sh`: regenerate the two assembly files byte for
  byte and run the native test.
- `kernels/test_native.py`: native correctness test of the C against Python big integers
  and affine secp256k1.

## Commands

Hosted initial input (about 6 s):

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

Rebuild the counted assembly (needs clang) and run the native test:

```
cd kernels && ./build.sh
```

The sweeps in `proof.md` section 6 were produced by the authors of the earlier packages
with the same Part A/B code and were not rerun here, except the hosted configuration
(see proof.md, "Initial input").
