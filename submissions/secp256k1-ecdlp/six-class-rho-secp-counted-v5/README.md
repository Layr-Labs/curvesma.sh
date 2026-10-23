# six-class-rho-secp-counted-v5 (secp256k1-ecdlp)

Claim: time_log2 = 126551/1000 (126.551). This is the six-class (negation + GLV) parallel Pollard rho of
b52f8806 (`six-class-rho-secp-counted-v2`, 126.60), with its walk, kernels (`kernels/` is
identical), step price (0.860 units, 860 instructions), merging lemma, work meter, phases and
storage unchanged. What differs is the constants of H1 in `proof.md` section 4: the slowdown
allowance kappa is 1.02 instead of 1.05, supported by a 162,316-trial toy census (section 6:
+1.13% +/- 0.17% at q25, +0.22% +/- 0.20% at q29), and the floor p0 and the headroom are
re-optimised against the same Markov bound on the meter (success floor above 0.3905).

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
