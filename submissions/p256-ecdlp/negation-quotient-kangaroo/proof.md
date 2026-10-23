# p256-ecdlp: negation-quotient kangaroo with a charged instruction ledger

Status: ready for exploratory review; no qualification, emitted score or practical
break is claimed. The normative inputs/outputs are those of the protected
`target-profile.json`; `claim.json` binds that profile and
`classical-group-work-v2` by SHA-256. No restricted challenge distribution,
uncharged advice, random-oracle premise or quantum step is used.

## 1. Target

The protected profile fixes NIST P-256 (SEC 2 v2.0 Section 2.4.2):

- `p = ffffffff00000001000000000000000000000000ffffffffffffffffffffffff`
- `a = ffffffff00000001000000000000000000000000fffffffffffffffffffffffc` (so `a = -3`)
- `b = 5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b`
- subgroup order `r = ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551`
  (prime), cofactor 1, generator `P` as published.

The input is only `Q = [x]P` with `x` uniform in `[0, r)`. The identity target is
included: `x = 0` exactly when `Q` is the point at infinity, and the algorithm
returns `0` for that canonical input without a walk. The output is `x`; the
attack-side check is `[x]P = Q` after canonical decoding and the `[r]Q = infinity`
subgroup check.

## 2. Mechanism: kangaroo walk on the negation quotient

Let `R -> -R` be the order-two automorphism of the group. For every integer `k`,
`x([k]P) = x([-k]P)`, so the pair `{k, -k}` is a single class. Define the
canonical representative

```text
can(R) = R            if y(R) <= (p-1)/2
can(R) = -R           otherwise          (negate the y-coordinate mod p)
```

and let `w = ceil(r/2)` be the number of classes, `sqrt(w) < 2^127.5`.

**Jump function.** Fix `JN = 2^8` public jump sizes `s_0..s_{JN-1}` with mean
`m = floor(sqrt(w)/2)`, and the public points `[s_j]P`. The index of a state is
`idx(R) = can(R).x mod JN` (the low byte of the canonical abscissa). One step is

```text
R  <-  can( R + [ s_idx(R) ] P ).
```

The index depends only on the canonical class, so both animals run the same
deterministic walk on the quotient, and a step is one mixed Jacobian addition
followed by a conditional y-negation.

**Tame animal.** Start at `A0 = [a0]P` with `a0 = floor(r/2)`. Perform
`T = 2^128` steps. Before each step, insert the 32-byte canonical abscissa into
an open-addressing hash table with load factor at most 1/2, keyed by that
abscissa; the entry retains the tame scalar position `k_t`.

**Wild animal.** Start at `Q`. At each step, look the canonical abscissa up in
the table. If it is present with tame position `k_t` and the wild has accumulated
offset `w_off`, then

```text
x + w_off  ==  +- k_t   (mod r)
```

so the two candidates are `x = (k_t - w_off) mod r` and `x = (-k_t - w_off) mod r`.
Both are verified by computing `[x]P` and comparing with `Q`; the matching one is
output. The wild is allowed `2T` steps. If no verified candidate appears inside
that budget the attack reports failure; the budget is charged in full either way.

**Correctness.** `Q = [x]P`, and the wild's position after `w_off` total jump
distance is `[x + w_off]P`, whose class is tracked exactly by the walk (each step
is an addition of the public point `[s_j]P` composed with the canonical
projection). If the wild's class equals a stored tame class, then
`x + w_off = +- k_t (mod r)`, so one of the two candidates equals `x`, and the
independent check `[x]P = Q` accepts it. No other candidate is accepted. The
identity case `x = 0` is handled by the canonical input check before the walk.

## 3. Success probability

The unknown distance `d = x - a0` lies in `(-r/2, r/2)` and is uniform, because
`x` is uniform in `[0, r)`. The quotient identifies `d` with `-d`, so `d` is
uniform over `w = ceil(r/2)` classes. The tame trail of `T >= sqrt(w)` stored
positions spans the interval with mean jump `m = floor(sqrt(w)/2)`, and the wild
walks forward through it.

The standard analysis of Pollard's kangaroo (Pollard 1978; van Oorschot and
Wiener, *Parallel collision search with cryptanalytic applications*, J. Cryptology
12(1), 1999, Section 3; Galbraith, *Mathematics of Public Key Cryptography*,
Section 14.2.5) gives an expected total number of steps of at most `2*sqrt(w)` for
a fully stored trail and mean jump `sqrt(w)/2`. Here `sqrt(w) < 2^127.5`, so

```text
E[steps]  <=  2 * sqrt(w)  <  2^128.5
```

The charged budget is `B = 3T` steps (`T` tame plus `2T` wild), so
`B = 3*2^128 = 1020847100762815390390123822295304634368 < 2^129.585`. Markov's
inequality over the uniform challenge distribution gives

```text
Pr[failure]  <=  E[steps] / B  <  2^128.5 / (3*2^128)  =  2^-1.085  <  0.472
Pr[success]  >=  0.528  >  0.39.
```

The claim states the weaker bound `0.5` to keep margin. The expected-step estimate
is the single heuristic premise of this package and is disclosed in Section 6;
the collision-to-candidate resolution, the verification and the budget are exact.
No restarts are required by the stated bound, and none are charged.

## 4. Instruction ledger

All work is counted in AArch64 instruction equivalents and divided by the track's
fixed constant `C = 1000` from `classical-group-work-v2`. This package does **not**
borrow the secp256k1 six-multiplication addition as a P-256 price: every P-256
field and point operation below is charged through its own instruction bound, and
the borrowed constant is used only as the normalization divisor the track fixes.

### 4.1 Primitive bounds

| Kernel | Instructions | Derivation |
| --- | ---: | --- |
| 4x64-bit Montgomery multiplication, CIOS | `M = 256` | 16 `MUL`/`UMULH` pairs = 32; 16 carry-accumulate steps <= 64; 4 Montgomery rounds of 4 products, accumulation and `q` <= 96; final conditional subtraction <= 16; 12 operand loads/stores and addressing <= 48 |
| modular addition/subtraction | `A = 32` | 4-limb add with carry, conditional subtract of `p`, 8 loads/stores |
| canonical y-negation | `8A = 256` | 4-limb comparison against `(p-1)/2`, conditional 4-limb negation, selection |
| mixed Jacobian addition, `a = -3`, `Z2 = 1` | `8M + 3S <= 12M = 3072` | standard `8M+3S`; charged at the general `12M+4S` rate for margin, `S = M` |
| per-step overhead | `256` | `x mod 256` (2), jump-point load (8), one hash probe (16), 32-byte store (8), loop control and addressing (8), slack |
| **one step** | `J = 3584` | `3072 + 256 + 256` |

`M = 256` is conservative against the cost model's own pinned measurement of 179
instructions for the 5x52 secp256k1 multiplication: a 4x64 CIOS multiply executes
16 64x64->128 products plus four reduction rounds and a final subtraction, which
is the same order with fewer operand limbs.

### 4.2 Charged components

`T = 2^128 = 340282366920938463463374607431768211456` and
`B = 3T = 1020847100762815390390123822295304634368` steps.

| Component | Instructions |
| --- | ---: |
| jump table: 256 scalar multiplications, each <= 384 steps of `J` | 352321536 |
| tame trail: `T` steps of `J` | 1219572003044643453052734593035457269858304 |
| tame inserts: `T` table insertions at 64 instructions | 21778071482940061661655974875633165533184 |
| wild search: `2T` steps of `J` | 2439144006089286906105469186070914539716608 |
| verification: two candidates, each one scalar multiplication `<= 2^21` | 4194304 |
| decoding, curve/subgroup checks, output | 1048576 |
| **preprocessing subtotal** | **1241350074527583514714390567911090787713024** |
| **online subtotal** | **2439144006089286906105469186070914544959488** |
| **total** | **3680494080616870420819859753982005332672512** |

Dividing by `C = 1000` gives `W = 3.680494080616870420819859753982005332672512e39`
and

```text
time_log2          = 131.435096   (rounded upward at six decimals)
preprocessing_log2 = 129.867106   (rounded upward at six decimals)
```

The preprocessing subtotal (jump table, whole tame trail and its insertions) is
included in the total with no amortisation or reuse credit.

### 4.3 Memory and advice

The trail holds `2^128` abscissas of 32 bytes, `2^133` bytes. The open-addressing
table at load factor at most 1/2 occupies `2^134` bytes. Code, public constants,
jump points, walk state and scratch are bounded by `2^20` bytes. Peak storage is
`2^134 + 2^20` bytes, `memory_log2_bytes = 134.000001`. The stored trail is
nonuniform advice of `2^133` bytes, `nonuniform_advice_log2_bytes = 133.000000`,
constructed by the charged preprocessing and included in peak storage. This
memory is astronomical and mandatory for this construction; it is not a practical
attack.

### 4.4 Accounting completeness

Every item in the cost model's `total_time_includes` is charged exactly once:

- preprocessing and advice construction: the jump table and the full tame trail;
- all trials and failures: the complete `2T` wild budget is charged regardless of
  when the collision is found, and the failure branch costs no more;
- randomness and rejection sampling: the construction is deterministic and the
  challenge is sampled by the organizer, so no attack-side sampling is charged;
- sorting, lookup and memory transfers: one charged probe per step and one
  charged insertion per tame point;
- output construction and attack-side verification: two scalar multiplications
  and the `[r]Q`/`[x]P=Q` checks;
- restarts and amplification: none are needed for the stated success bound.

Work is not split across processors, so no cross-processor summation applies. No
work is amortised over future targets, and no oracle or advice is left uncharged.
Public profile constants and algorithm design effort are excluded as the rules
allow.

## 5. Evidence

`experiment.py` runs the same quotient walk, the same collision-to-candidate
resolution and the same verification on the additive cyclic group `Z_n` for every
scalar `0 <= x < n`, with `n <= 5000` and a deterministic seed. It exhausts all
`n` scalars and asserts recovery. Its scope is the toy group only: it exercises
the quotient bookkeeping, the two-candidate resolution and the identity element,
and it does not establish production cost, memory or the success bound.

The production figures are analytic: no full-size execution is claimed, and the
exponential table is not built. No organizer fixture is needed beyond the
protected profile constants, which are used exactly as published.

## 6. Heuristic disclosure and limitations

Heuristic premise (the only one): the expected number of kangaroo steps for a
fully stored trail with mean jump `sqrt(w)/2` is at most `2*sqrt(w)`. Support:
the published analyses cited in Section 3, which is the standard result for this
walk; the trail is fully stored, so no distinguished-point loss factor applies.
Limitations and extrapolations:

- the estimate is an expectation, not a worst case; Section 3 converts it into a
  success bound with Markov's inequality rather than asserting success with
  probability 1;
- the instruction bounds are analytic upper bounds, not measurements, and are
  conservative relative to the pinned reference measurement;
- the `r`-adding walk's mixing is the standard random-walk heuristic;
- peak memory `2^134` bytes is not physically realisable, so this is an
  accounting bound, not a practical break;
- the submitted bound is the total work of this concrete construction; it is not
  a best-known-attack estimate and does not claim novelty of the kangaroo method
  or a weakened hardness assumption.

## 7. Submitted bound

```text
time_log2                       = 131.435096
time_unit                       = target-group-operations
success_probability             >= 0.5
preprocessing_log2              = 129.867106   (included in total)
memory_log2_bytes               = 134.000001
nonuniform_advice_log2_bytes    = 133.000000   (included in peak storage)
```

References: J. M. Pollard, *Monte Carlo methods for index computation (mod p)*,
Math. Comp. 32 (1978); P. C. van Oorschot and M. J. Wiener, *Parallel collision
search with cryptanalytic applications*, J. Cryptology 12(1) (1999); S. Galbraith,
*Mathematics of Public Key Cryptography*, Section 14.2.5; E. Brier and M. Joye,
*Weierstrass elliptic curves and side-channel attacks*, PKC 2002; SEC 2 v2.0
Section 2.4.2. The construction and accounting above are self-contained; review
does not need to fetch the citations to reconstruct them.
