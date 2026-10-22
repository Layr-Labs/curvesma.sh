# bls12-381-g1-ecdlp: Cheon p−1 attack at n=2^20

Status: ready for exploratory review. Analytic classical attack using the
auxiliary powers and the number-theoretic precondition stated in the track
rules. Not a practical break or full-size execution.

## Target and precondition

Let r be the BLS12-381 G1 subgroup order from the protected profile. The track
rules and `target-profile.json` state that d = 2^20 divides r−1, that A_d is
supplied at the scoring parameter n = 2^20, and that Cheon-type methods are
admissible. Explicitly (r−1) is divisible by 2^32, hence by d = 2^20. Scoring
asks for log2(W(2^20)).

Inputs at n = 2^20: public P0, T = A_1 = [s]P0, and A_i = [s^i mod r]P0 for
i = 1..n. In particular A_d = [s^d]P0 is available. Recover s ∈ [1,r).

## Algorithm (Cheon p−1 BSGS, addition form)

Follow Cheon’s p−1 algorithm (JoC 2010 / Eurocrypt 2006) with the group-addition
implementation below. Let ξ be a primitive root modulo r (found offline from the
known factorisation of r−1; public profile constant work). Write

```text
s = ξ^{((r-1)/d)·k1 + k2}   with 0 ≤ k1 < d, 0 ≤ k2 < (r-1)/d.
```

Then s^d ≡ ξ^{d·k2} (mod r). Set ζ = ξ^d mod r.

**Phase 1 — recover k2.** Let m1 = ⌈√((r−1)/d)⌉. Build a hash table of baby
steps in G1:

```text
B ← [ζ]P0                 # one scalar mul
R ← A_d                   # = [s^d]P0 = [ζ^{k2}]P0
for v = 0..m1-1:
    store (R, v) keyed by R;  R ← R + (−B)   # i.e. walk A_d · ζ^{−v}
G ← [ζ^{m1}]P0            # one scalar mul
S ← O (identity)
for u = 0..m1-1:
    if S hits the table at v: k2 ← m1·u + v; break
    S ← S + G
```

This is the standard BSGS match for ζ^{k2} = ζ^{m1·u + v}. Cost: 2 scalar
muls to build B and G, then at most 2·m1 G1 additions and 2·m1 hash
insert/lookup operations.

**Phase 2 — recover k1.** After k2 is known, set

```text
T' ← [ξ^{−k2}] · T = [ξ^{((r-1)/d)·k1}] P0.
```

Let δ = ξ^{(r−1)/d} mod r and m2 = ⌈√d⌉. Run the same BSGS for k1 against
base [δ]P0 and target T', costing 2 scalar muls and at most 2·m2 G1 additions.

**Output.** s ← ξ^{((r−1)/d)·k1 + k2} mod r. Verify [s]P0 = T.

Success probability: 1 on every valid challenge (deterministic BSGS covers the
full ranges of k1, k2). Hence ≥ 0.39. Heuristics: none for correctness; the only
heuristic is that a random search finds a primitive root ξ quickly, which is
standard and is charged in setup.

## Cost W(n) at n = 2^20

Reference operation: batched incomplete affine G1 addition, C = 2,000.
Each G1 addition below is one reference unit. Scalar multiplication by a
≤255-bit scalar is charged as at most 512 reference additions (double-and-add
with a generous margin). Hash-table work is counted in instruction equivalents
and divided by C.

Exact integer parameters:

```text
d   = 2^20 = 1048576
m1  = ceil(sqrt((r-1)/d)) = 223621884914696835096185700517675008
m2  = ceil(sqrt(d)) = 1024
```

| Component | Reference G1 additions (upper bound) |
| --- | ---: |
| Setup (≤32 scalar muls for primitive-root trials and scratch) | 32 × 512 = 16384 |
| Phase-1 bases B, G | 2 × 512 = 1024 |
| Phase-1 baby + giant walks | 2 × m1 |
| Phase-2 bases | 2 × 512 = 1024 |
| Phase-2 baby + giant walks | 2 × m2 |
| Verification [s]P0 and [r]T checks | 2 × 512 = 1024 |
| Hash insert/lookup instruction budget (250 Ir × (2 m1 + 2 m2)) / C | (500·(m1+m2))/2000 = (m1+m2)/4 |

```text
H_total = 16384 + 1024 + 2·m1 + 1024 + 2·m2 + 1024 + (m1+m2)/4
        = (9/4)·m1 + (9/4)·m2 + 19456
```

```text
log2(H_total) = 118.754195   (rounded UP to six decimals from the exact rational)
```

So W(2^20) ≤ 2^118.754195 and the submitted score is **118.754195**.
Preprocessing (setup + bases) ≤ 2^15; reported `preprocessing_log2 = 15.0`.
Peak memory: one hash table of m1 G1 keys with 96-byte records plus scratch
⇒ `memory_log2_bytes = 125.428545` (rounded up). Nonuniform advice: 0.

Applicability: the bound holds for every n ≥ d = 2^20 at which A_d is present
and d | (r−1); the scoring point n = 2^20 satisfies both. Using more than A_1
and A_d is unnecessary for this claim.

## Comparison

A generic automorphism-aware Pollard rho on BLS12-381 G1 without powers sits
near 2^127 group-equivalents. Cheon’s reduction replaces that with a BSGS of
size √((r−1)/d) ≈ 2^117.43 plus a √d ≈ 2^10 epilogue, which is exactly the
improvement the track’s auxiliary powers are designed to enable.

## Experiments

Optional toy experiment `experiment.py` runs Cheon’s p−1 BSGS on a small prime
order group (integers mod r under addition as a stand-in) with a chosen d | (r−1),
checking full recovery. Toy scope only; does not establish production G1 cost.

## References

- J. Cheon, Discrete logarithm problems with auxiliary inputs, JoC 2010.
- J. Cheon, Security analysis of the Strong Diffie–Hellman problem, Eurocrypt 2006.
- Track rules: “d=2^20 divides r−1 … Cheon’s attack using these inputs is admissible”.
