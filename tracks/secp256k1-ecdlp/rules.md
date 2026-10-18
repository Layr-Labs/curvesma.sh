# secp256k1-ecdlp

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
    "id": "secp256k1-ecdlp",
    "curve": "secp256k1",
    "problem": "ecdlp",
    "scope": "unrestricted-production",
    "field_modulus": "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f",
    "subgroup_order": "fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",
    "curve_equation": "y^2 = x^3 + 7 over Fp",
    "cofactor": 1,
    "generator_P": {
      "x": "79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",
      "y": "483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8"
    },
    "challenge": "Sample x uniformly from 0..r-1 by rejection sampling; give only Q=[x]P. The scalar and sampling coins stay private.",
    "input": {
      "Q": "secp256k1 point"
    },
    "output": {
      "x": "32-byte scalar"
    },
    "acceptance": "Decode Q canonically, require Q in <P>; decode 0<=x<r; accept exactly when [x]P=Q.",
    "encoding": "UTF-8 JSON with exactly the stated keys. Scalar: 64 lowercase hex characters, no prefix. Point: null for infinity, otherwise {x,y}, each 64 lowercase hex characters representing a big-endian integer less than p. Reject duplicate keys, noncanonical values and off-curve points.",
    "source": "https://www.secg.org/sec2-v2.pdf#page=13"
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
