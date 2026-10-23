# six-class-rho-secp-counted-v7 (secp256k1-ecdlp)

Claim: time_log2 = 126403/1000 (126.403). This is the six-class (negation + GLV) parallel Pollard rho of
30f6f34e (`six-class-rho-secp-counted-v6`, 126.423) and its ancestors (98e7729c, b52f8806, a9f1a6fe, 3682292, 758d8df6,
189ebc42), with the walk, merging lemma, work meter, phases, storage and H1 constants (kappa = 1.02) unchanged.
What differs is the price of the doubling escape: 1.1 units (967 instructions bound) derived from the counted
operation costs of the generated assembly in `proof.md` section 5.3, instead of 2.0 from the old allowance. The step
(0.785 units = 785 instructions) and the two counted loops (111 + 534 instructions per walk, executed natively and by an
interpreter inside `experiment.py`) are unchanged; p0 and the headroom are re-optimised against the same Markov
bound. See `proof.md`, sections 0 and 5.

## Files

- `proof.md`: the complete argument and ledger.
- `experiment.py`: standard library only; reads a JSON config on stdin. Parts A and B are the code of 3682292
  (production secp256k1 identities; toy-curve walks). The `"kernels": true` option runs limb-for-limb Python ports
  of the C kernels, an AArch64 interpreter that executes the shipped `secp_loops.s` on production secp256k1 points,
  parses the shipped assembly (instructions per iteration of the single-block loops), and recomputes the ledger.
- `kernels/gen_loops.py`: generator of the counted AArch64 loops (`elf` flavor = the shipped file, `macho` flavor
  for native tests on macOS, `tests` option = unit-test wrappers of the four emitters).
- `kernels/secp_loops.s`: the counted loops pass1/pass2 (generator output, `elf` flavor).
- `kernels/secp_walk.c`: freestanding C specification: fe_mul, fe_sqr, fe_sub, fe_canon, pass1_c, pass2_c, the
  libsecp256k1 inversion chain, test drivers.
- `kernels/secp_fe.s`: clang's AArch64 output for fe_mul_f/fe_sqr_f/fe_inv (the inversion bound).
- `kernels/strip_asm.py`, `kernels/build.sh`: regenerate the two assembly files byte for byte and run the native tests.
- `kernels/test_native.py`: native correctness test of the C and of the assembly against Python big integers and
  affine secp256k1.

## Commands

Hosted initial input (about 8 s):

```
{"curve":"q25","trials":2000,"R":2048,"dp_bits":5,"seed":1,"production_checks":2,"seconds":6,"pipeline":true,"kernels":true}
```

Rebuild the shipped assembly (python3 and clang) and run the native tests (AArch64 Linux or macOS):

```
cd kernels && ./build.sh
```

The sweeps in `proof.md` section 6 were produced by the authors of the earlier packages with the same Part A/B
code and were not rerun here, except the hosted configuration and the slowdown census.
