# plain-rho-v6: secp256k1-ecdlp candidate

v5's negation walk is replaced by the plain additive walk on the full
group with mixed-radix constants c_j = jQ + 2^j mod r. Fruitless
cycles are excluded by an exact theorem; all cycle-detection
machinery is deleted.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **16058573/125000 = 128.468584** |
| success_probability (lower bound) | 0.39 (modeled 0.4066) |
| memory_log2_bytes | 105.1 (over-report) |
| preprocessing_log2 | 22.0 (included) |
| nonuniform_advice_log2_bytes | 0 |

## Why the plain walk

The v5 review found a real error: the census assumed per-symbol even
multiplicity, but 2^j + 2^j = 2^(j+1) gives length-3 carry families
that the enumeration never saw. v6 removes the negation map, so the
composed map of any L steps is C -> C + [sum c_{w_u}]P: with c_j =
jQ + 2^j, a return forces 32B = rho*A with B >= 2A - impossible for
any length up to 2^220.00 (proof.md). No census, no look-back,
no ratchet, no tail bounds.

The j*Q term is not decoration: pure 2^j constants do not mix at
production scale (the scalar difference walk stays in a 2^143 band
and cross-walker collisions essentially never happen). With the j*Q
term the constant differences span 0.9688 of the modulus,
the standard parallel-rho birthday model applies, and the toy
measures 0.867 of the ideal.

## Reproduce

```
echo '{"trials": 200, "dbits": 0, "walkers": 8}' | \
  python3 experiment.py
```

Stdlib only; deterministic; single JSON object on stdout. The
declared run also records the dbits=6 sparse-DP census.
