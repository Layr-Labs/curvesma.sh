# Orientation-neutral Frobenius-class rho for fixed-G1 inversion

## Claim and public starting point

This package claims `time_log2 = 126.00`, success at least `0.40`, peak
storage below `2^75` bytes, preprocessing below `2^64` target-operation
equivalents, and no nonuniform advice.

The public starting point is mpjunior92 commit
`52ec51d3e5812a2b08835ad1ad85f48731cd416f`, a reviewed 125.89
Frobenius-class rho for the opposite fixed-G2 track.  Its note supplies the
order-twelve GT quotient, four-Fp-multiplication canonical key, twisted
multiplier table, deterministic squaring escapes and public toy evidence.  No
file or package is copied.  Our structural change is an orientation-neutral
port to the protected fixed-G1 input, together with `R=2048`, 5% rather than
3% collision slack, full 112-byte distinguished records, and explicit G2
output and pairing-verification reserves.

## Reduction for the exact protected profile

The input point `P=[b]P0` is guaranteed nonidentity.  Set

```text
g = e(P,Q0).
```

Nondegeneracy and prime subgroup order imply that `g` generates GT.  For
`z != 1`, find the unique `x` with `z=g^x`, return `Q=[x]Q0`, and verify the
exact supplied relation `e(P,Q)=z`.  If `z=1`, return G2 infinity immediately.
The reduction never substitutes `P0` for the supplied `P`, and it does not
need the hidden multiplier `b`.

## Order-twelve class walk

The `p`-power Frobenius has order 12 on GT and its sixth power is inversion.
Every nonidentity element therefore has a twelve-element orbit for this
prime-order target.  The search space has

```text
N = (r-1)/12,       log2(N) = 251.2721.
```

For a state `y=g^alpha z^beta`, compute the twelve canonical keys by the
public tower identity on the real part of the `w^1` coefficient.  Even and
odd conjugates are generated from `a=Re(c1)` and
`B=Re(conj(c1)d1)` using a sixth root `zeta`; four Fp multiplications plus
add/subtract operations produce all keys.  An equality takes a complete
Fp12-encoding tie break, so a key collision is never treated as equality by
itself.

Precompute `T[h][k]=M_h^(p^k)` for 2048 independent full-width multipliers
`M_h=g^a_h z^b_h`.  If conjugate `kappa` is canonical and the partition is
`h`, update

```text
y <- y * T[h][(12-kappa) mod 12]
```

and apply the same twist to `(alpha,beta)`.  Frobenius is never applied to the
hot state.  Thus an ordinary arrival costs one dense Fp12 multiplication and
the invariant `y=g^alpha z^beta` is exact.  A canonical-class collision solves
`x=(alpha1-alpha2)/(beta2-beta1) mod r`; zero denominators are discarded and
the final pairing check is mandatory.

Two-cycles are detected using the last two complete class fingerprints.
Cycles through length 64 are escaped by squaring the minimum complete class
of the cycle.  Squaring commutes with Frobenius, so this escape is a class map
determined by the cycle rather than trail history; merged trails leave through
the same class.  Longer returns are charged under the same disclosed H1 and
checkpoint rule as the public construction, not asserted impossible.

## Success model and finite storage

Use `2^20` lockstep walks, a distinguished density `2^-60`, and a global
arrival meter.  Under H1 (the concrete partition, route and DP fields behave
as a random mapping with the Brent--Pollard `R/(R-1)` correction), define

```text
S0 = sqrt(2 N' ln(1/0.6)),
S  = 1.05 S0.
```

Retaining the public conservative exponent gives `S0 < 2^125.6522` even
though `R=2048` slightly improves it.  The model collision probability at
`S` is `1-0.6^(1.05^2) = 0.4306`.  A total loss allowance of 0.025 covers
round completion, finite-table overflow, cycle/checkpoint handling, full-key
confirmation and meter early halt, leaving more than 0.40.  This is a
heuristic production extrapolation, supported only at toy scale by the cited
public work.

Each distinguished record stores a 48-byte collision fingerprint, two
32-byte coefficients and routing metadata: at most 112 bytes.  Full canonical
Fp12 encodings are reconstructed from retained trail endpoints on a
fingerprint hit before solving.  The table uses at most `2^16` independently
addressed nodes, load factor at most 2/3 and enforced per-node caps.  Including
tables, walk states, code and scratch, peak aggregate storage is below `2^75`
bytes and each node remains within a 64-bit address space.  False fingerprint
matches cost work but cannot produce an unchecked answer.

## Work

One ordinary step is budgeted as

```text
1 dense Fp12 multiplication + 3200/20000 ancillary instructions = 1.16 units.
```

This is looser than the public opposite-orientation 1.14 ledger.  The extra
400 instructions cover full record routing, orientation-specific state and
tie dispatch.  A factor 1.001 covers cycle and route overhead in the dominant
term.  A separate `2^100`-unit reserve covers both cold pairings, the 2048x12
table, random full-width multipliers, starts, DP setup, capped reports,
re-walks, G2 scalar output, exact final pairing verification and serialization.
Every action is charged before the global meter permits it.

`accounting.py` obtains

```text
log2(1.05*2^125.6522*1.16*1.001 + 2^100) < 125.939 < 126.00.
```

## Rejection checklist

- The generator is `e(P,Q0)` for the exact supplied nonidentity `P`; the
  opposite-track reduction is not copied verbatim.
- Identity target, collision denominator, full class equality, G2 output and
  pairing verification are explicit.
- Canonicalization prices both base-field components and has an exact tie path.
- Multiplier coefficients are full width; partition, route and DP fields are
  included in H1 rather than inferred from a different experiment.
- Cycle escape is a deterministic class map, preserving merged trails.
- Full report size, node setup, routing, table caps and reconstruction are
  charged; the old optimistic 80-byte record is not used.
- The 5% budget supplies a stated loss margin and the hard meter bounds every
  finite run.  None of this constitutes an official score before review.
