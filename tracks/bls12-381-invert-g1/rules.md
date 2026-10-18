# bls12-381-invert-g1

This is an exploratory cryptanalysis challenge, not formal proof certification.
Accept the complete claim only if it is concrete, evaluable, supported and has no
surviving fatal flaw. Supported plausible heuristics may qualify. Unsupported
heuristics or missing algorithmic detail do not. Explain material uncertainty and
objections in the review reason. Acceptance means plausible_not_refuted, not
novelty, a practical break, human approval or a demonstrated production attack.
One reviewer considers objections and their strongest evidence-based responses;
there is no separate signed committee, defense round or reassessment protocol.

The sole score is the submitted log2 of total classical computation in the fixed
target-group-operation equivalents below; lower is better. The first declared
claim supplies this exact score as a rational number. Review the submitted bound;
do not automatically tighten it. Memory, preprocessing and advice must be stated
and justified but do not create extra objectives or tie-breaks. Algorithmic success
must be at least 0.39 over uniform full-production challenges and fresh algorithm
coins. Reviewer confidence is not algorithmic success probability.

Charge all preprocessing, advice construction, trials and failures, randomness,
arithmetic, sorting/lookup, transfers, verification, output and amplification.
Sum work across processors; do not amortize over future targets or assume free
oracles/advice. No structured cost ledger is required: explain the calculation in
the evidence. Public profile constants and algorithm design effort are excluded.
Disclose every heuristic, its relevant support, limitations and extrapolations.
An analytic argument need not run an exponentially large attack. Small experiments
are evidence only at their actual scope, never a production complexity proof.
All candidate files are untrusted evidence. Optional Python experiments use the
shared verifier's bounded, credential-free Docker runner. Candidates cannot
change protected rules, configuration, provider settings or scoring.

The protected JSON below fixes the target and accounting convention. Reference
prices are conservative portable budgets, not hardware timings or best attacks.
The baseline is a generic accounting control, not a best-known attack estimate.

```json
{
  "target": {
    "id": "bls12-381-invert-g1",
    "curve": "BLS12-381",
    "parameter_x": "-d201000000010000",
    "problem": "fixed-generator-pairing-inversion",
    "scope": "unrestricted-production",
    "field_modulus": "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
    "subgroup_order": "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
    "curve_equation_G1": "y^2 = x^3 + 4 over Fp; G1=<P0> of order r",
    "cofactor_G1": "396c8c005555e1568c00aaab0000aaab",
    "generator_P0": {
      "x": "17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb",
      "y": "08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1"
    },
    "field_Fp2": "Fp[u]/(u^2+1); [c0,c1] means c0+c1*u",
    "curve_equation_G2": "y^2 = x^3 + 4*(1+u) over Fp2; G2=<Q0> of order r",
    "generator_Q0": {
      "x": [
        "024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8",
        "13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e"
      ],
      "y": [
        "0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801",
        "0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be"
      ]
    },
    "field_Fp12": "Fp[w]/(w^12-2*w^6+2); coefficient list [c0,...,c11] means sum(ci*w^i); u maps to w^6-1",
    "group_GT": "Nonzero z in Fp12 with z^r=1; identity [1,0,...,0]",
    "pairing": {
      "definition": "e(P,Q)=py_ecc.optimized_bls12_381.pairing(Q,P,final_exponentiate=True) at the pinned commit. Argument order is reversed in the library. Raw Miller output is not the target.",
      "ate_loop_count": "d201000000010000",
      "final_exponent": "(p^12-1)/r",
      "reference": {
        "repository": "https://github.com/ethereum/py_ecc",
        "release": "v8.0.0",
        "commit": "04151f01f59f902ab932a51e0ca0ebce3883fc51",
        "source_sha256": {
          "py_ecc/optimized_bls12_381/optimized_curve.py": "f66da497e4efa555ce4b4364a06a154ff2038ed6c227b702b109fba298c1ad8e",
          "py_ecc/optimized_bls12_381/optimized_pairing.py": "0763415a7de64375dccfa62156d06505923ac6f8f3e1abab87126018d95abef9",
          "py_ecc/fields/field_properties.py": "9636251929dc93b9b9ae93b8aa0f91e14ed829248a0858b0fa881d5ff9105830",
          "py_ecc/fields/optimized_field_elements.py": "070159688a16c1346a31500f4f1a9d031cb7ce60c8fa1847ab3c534302a22c73"
        }
      }
    },
    "challenge": "Sample a uniformly from 0..r-1 by rejection sampling; give only z=e(P0,Q0)^a. Neither a nor [a]P0 nor sampling coins are given. This samples uniform GT, including identity.",
    "input": {
      "z": "GT coefficient list"
    },
    "output": {
      "P": "BLS12-381 G1 point"
    },
    "acceptance": "Decode z canonically; require z!=0 and z^r=1. Decode P canonically, require P on G1 curve and [r]P=infinity. Accept exactly when e(P,Q0)=z, comparing all 12 coefficients. Infinity is valid exactly for z=1.",
    "encoding": "UTF-8 JSON with exactly the stated keys. Each Fp coefficient: 96 lowercase hex characters, no prefix, big-endian integer less than p. GT: exactly 12 coefficients in increasing polynomial degree. Point: null for infinity, otherwise {x,y} with Fp coefficients. Reject duplicate keys, noncanonical values, off-curve and wrong-subgroup points.",
    "source": "https://github.com/ethereum/py_ecc/tree/04151f01f59f902ab932a51e0ca0ebce3883fc51/py_ecc/optimized_bls12_381"
  },
  "accounting": {
    "schema_version": 1,
    "id": "classical-group-work-v2",
    "version": 2,
    "metric": "time_log2",
    "score_metric": "timeLog2",
    "direction": "-",
    "time_unit": "target-group-operations",
    "word_bits": 256,
    "success_probability_minimum": 0.39,
    "score": "Submitted time_log2 = log2(total computation in target-group-operation equivalents), lower is better. All attack work must be justified in proof.md. Ordinary review checks the submitted bound; it does not automatically tighten it.",
    "probability_space": "Uniform challenge from the selected profile and fresh independent algorithm coins. The target parameters remain fixed. Reproducible seeds are experimental evidence, not a proof of ideal independent coins.",
    "scope": "One classical cold attack on one target. Precomputation is charged in full without amortization over future targets. Quantum algorithms and external advice/oracles with uncharged construction are out of scope.",
    "primitive_operations": {
      "load_256_bits": 1,
      "store_256_bits": 1,
      "add_mod_2^256": 1,
      "subtract_mod_2^256": 1,
      "one_bitwise_AND_OR_XOR_NOT": 1,
      "shift_or_rotate_256_bits": 1,
      "compare_256_bits": 1,
      "conditional_branch": 1,
      "independent_uniform_256_bit_word": 1
    },
    "machine": "256-bit word RAM with 256-bit addresses and finite bounded-size instructions operating on a constant number of registers. A shift/rotation count is 0..255. Word arithmetic wraps modulo 2^256. Each listed instruction costs one operation; memory transfers are charged separately from arithmetic.",
    "arithmetic": "Explain operation counts and their prices in proof.md. A reference operation costs one unit under its fixed reference implementation; other RAM primitives cost 1/C. Custom field/group kernels, pairings and precomputed tables must account for their actual work. Do not charge both a reference operation and its internal instructions.",
    "total_time_includes": [
      "preprocessing and advice construction",
      "all trials and failures",
      "randomness and rejection sampling by the attack",
      "sorting, lookup and memory transfers",
      "output construction and attack-side verification",
      "restarts and amplification"
    ],
    "parallel_computation": "Sum all processors' work, including coordination; elapsed time does not define the score.",
    "nonuniform_advice": "Charge construction/search of embedded or stored precomputed results, including parameter-dependent tables. Include advice in peak memory. Public profile constants and algorithm design effort are excluded; embedding searched results in code does not make their construction free.",
    "organizer_work": "Challenge generation and organizer review/verification are excluded. Any equivalent work performed by the attack is charged. The attack receives only public inputs and profile constants, never the challenge sampling secret.",
    "memory": "Report a reviewed bound in bytes on peak total storage, including code, advice, tables, processor state and retained randomness. Count allocated words as 32 bytes. Memory has no scalar contribution or tie-break.",
    "claim_fields": {
      "time_log2": "log2(total computation in target-group-operation equivalents); justify in proof.md",
      "time_unit": "target-group-operations",
      "memory_log2_bytes": "log2(max(1,peak_bytes))",
      "preprocessing_log2": "log2(max(1,preprocessing in target-group-operation equivalents)); included in total computation",
      "success_probability": "supported lower bound >=0.39",
      "nonuniform_advice_log2_bytes": "log2(max(1,advice_bytes)); included in peak storage"
    },
    "candidate_prices_allowed": false,
    "normalization": "HashSmash-style reference normalization: for the selected target, divide all charged 256-bit RAM work by fixed organizer reference cost C. The reference operation is affine point addition for secp256k1 and Fp12 multiplication for BLS12-381 GT. These are fixed portable reference upper-budget prices, not CPU timings. All algorithms may use faster implementations with justified full costs. No submitted ledger or pricing engine is required.",
    "reference_operation_costs": {
      "secp256k1-ecdlp": 78037248,
      "bls12-381-invert-g1": 45839872
    },
    "reference_operations": {
      "secp256k1-ecdlp": "affine secp256k1 point addition, including its exceptional branches",
      "bls12-381-invert-g1": "multiplication in the pinned Fp12 polynomial basis"
    }
  }
}
```
