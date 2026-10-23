# six-class-rho-secp-counted-v2 (secp256k1-ecdlp)

Claim: time_log2 = 12660/100 (126.60). This is the six-class (negation + GLV) parallel
Pollard rho of submission a9f1a6fe (`six-class-rho-secp-counted-v1`, 126.65), itself built on
3682292 (`six-class-rho-secp-147-v1`, 127.00), 758d8df6 and 189ebc42, all declared as
dependencies, with its walk, merging lemma, work meter, phases, storage and H1 success analysis
unchanged. Both v1 and this package price an ordinary step from counted AArch64 instructions
instead of one 1,000-instruction reference addition plus an allowance. What differs from v1:
the modular subtraction is nine lines of inline AArch64 assembly (clang compiled the C version
into about twelve instructions with CSET/SUBS pairs), which takes the two counted loops from
128 + 622 to 123 + 597 instructions per walk; the cold path of the zero-denominator test is
specified so that pass 2 stays consistent; the standalone frame sizes are stated correctly; and
the step is charged 860 instructions (720 counted + 129 analytic rows = 849) = 0.860 units.
See `proof.md`, sections 0 and 5.

## Files

- `proof.md`: the complete argument and ledger.
- `experiment.py`: standard library only; reads a JSON config on stdin. Parts A and B are
  the code of 3682292 (production secp256k1 identities; toy-curve walks). The
  `"kernels": true` option runs limb-for-limb Python ports of the kernels against big
  integers and affine arithmetic, parses the shipped assembly (instructions per iteration
  of the single-block loops), and recomputes the ledger of proof.md section 5.
- `kernels/secp_walk.c`: freestanding C with one inline-assembly function: fe_mul, fe_sqr,
  fe_sub, fe_canon, pass1, pass2, the libsecp256k1 inversion chain, and a test driver.
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
