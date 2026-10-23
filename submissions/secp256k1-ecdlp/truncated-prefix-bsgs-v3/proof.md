# secp256k1-ecdlp: truncated-prefix deterministic BSGS at exact 39/100 coverage

Status: ready for independent review; no qualification or emitted score.
Parent: commit `6423c2728dc03d5f2d647c7f22adf1dac0643a66`
(`deterministic-bsgs-v1`, time_log2 129.000802).

This revision changes one thing only: the baby-table extent `m` and the
giant-step extent `h` are retargeted from full coverage (`m = h = 2^128`) to
the least extents whose product reaches the exploratory success threshold
0.39. The kernels, the record layout, the sort, the lookup and the
verification charges are unchanged; the ledger below is the parent ledger
with `m` and `h` replaced by their truncated values, so every line can be
checked against the parent line by line.

The normative inputs/outputs are those of `production-v2`; `claim.json` binds
that profile and `classical-group-work-v2` by SHA-256. No restricted challenge
distribution, advice, random-oracle premise or quantum operation is used.

## Construction and success

Let p and r be the protected field modulus and subgroup order. Set
b = bit_length(r) = 256, L = floor(b/2) = 128, and let

```text
S = ceil(39*r/100)
m = ceil(sqrt(S)) = 212506270031153001307260662621714878127
h = ceil(S/m)       = 212506270031153001307260662621714878127
```

`m` is the least integer with `m^2 >= S`, hence `(m-1)^2 < S <= m^2`; and
`m*(m-1) < S`, hence `h = m`. The construction is the parent construction
restricted to these extents:

```text
Use P from the protected profile; input Q. Starting at infinity, store (jP,j)
for j=0..m-1 by repeated affine addition, retaining mP after the final update.
Set D=-(mP) and gamma=Q. Sort the table by (infinity_tag,x,y), tag 1 for infinity and 0 otherwise; use zero
coordinate words for infinity. For i=0..h-1, binary-search gamma; on a match
with j, set a=(i*m+j) mod r, verify [a]P=Q, and return a. Otherwise set
gamma=gamma+D. Reject malformed inputs; exhaust the fixed loops with failure.
```

For every valid Q there is a unique x in [0,r) with Q=[x]P. Write
x = i*m+j with 0 <= j < m and 0 <= i < h. At giant step i, gamma=[j]P; if
x < m*h the pair (i,j) exists and the table therefore contains a match.
Conversely every match yields a=(i*m+j) mod r with [a]P=Q, so no wrong
answer is returned: the attack-side equality check accepts only correct
scalars. Infinity gives x=0 and is represented explicitly.

## Exact coverage

Coverage is the prefix `[0, m*h)` of `[0,r)` with `m*h = m^2`. The integers
are fixed and checkable:

```text
r   = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
S   = 45158914802553316215192684153388284062606650068839212709216013625192082982792
m^2 = 45158914802553316215192684153388284062785719814745912697632069832470463028129
```

so `m^2 >= S`, `(m-1)^2 < S`, and `m*(m-1) < S <= m^2` all hold, and
`100*m^2 >= 39*r` holds as an integer inequality. For uniform x in [0,r) the
algorithm succeeds with probability `m^2/r >= 0.39` and fails with
probability at most 0.61. No randomness, coins, retries or amplification are
used: one deterministic attempt is charged, and the giant walk is charged in
full h steps whether or not a match occurs, so the bound below is a worst-case
single-attempt bound and not an expected-runtime shortcut. The claim is
therefore the threshold itself, 0.39, and nothing weaker is being
substituted for it.

## Word-operation ledger

Integers are little-endian arrays of s=ceil(k/256) limbs, where
k=bit_length(p)=256 and s=1. Each listed cost-model primitive costs one.
Operands of a field operation are canonical integers in [0,p).

| Kernel | Upper bound in word operations |
| --- | ---: |
| Modular addition/subtraction | A=128*(s+1)=256 |
| Base-field multiplication | M=k*(2*A+64)+128*(s+1)=147712 |
| Nonzero base-field inverse | I=2*k*(M+64)+128*(s+1)=75661568 |
| Affine point addition/doubling | C=I+16*M+32*A+4096=78037248 |

For A, zero-extend to s+1 limbs, add/subtract with carry/borrow, compare
against p, and make one conditional correction. Each limb pass needs at
most 16 operations for loads/stores, arithmetic and carry/borrow handling;
allow six passes plus 32*(s+1) for comparison, selection and control.
This accommodates the 257th carry bit even for the 256-bit field.

For M, scan the k bits of the multiplier; conditionally add the running
multiplicand to the accumulator and double it modulo p each iteration.
Two A calls and 64 bit-extraction/loop operations per bit suffice; reserve
128*(s+1) for initialization/output. Inversion uses binary powering by p-2:
at most 2*k multiplications plus 64 operations per bit. All instructions,
including memory transfers, are counted; multiplication/inversion are not
unit-cost primitives.

Point addition uses affine formulas: return the other operand for infinity;
return infinity for opposite points; use slope (v-y)/(u-x) for distinct
x-coordinates or 3*x^2/(2*y) for doubling, with y=0 handled as infinity.
Compute x'=slope^2-x-u and y'=slope*(x-x')-y modulo p. Negation is
(x,-y mod p). A scalar multiplication uses at most 2*b such operations.
C covers every case. Keys have K=3 words and records W=4 words. The
implementation performs one inverse, at most six multiplications and at most
eight additions per point addition, all inside the C coefficients.

Setup reserves 2^24 operations for decoding/constants and the step negation.
Verification/output reserves 4*b*C+2^24: this covers [r]Q and [a]P checks,
curve/encoding checks, scalar reconstruction by shift/add and serialization.

Use bottom-up merge sort with two m-record buffers and exactly L passes.
A lexicographic key comparison costs B=16*K+64; moving a record costs
D=16*W+64. These bounds allow all operand loads, stores, index arithmetic,
comparisons and branches, up to 16 operations per word plus 64 fixed work.
Each merge pass has at most m comparisons and 2*m record moves, with
128*m additional loop work. Thus sorting costs
S_sort=m*L*(B+2*D+128)+4096*L. Binary search performs at most L+1 key
comparisons; V=(L+1)*(B+64)+D+256 covers lookup and the returned record.
No constant-time hashing/table oracle is assumed.

For group-step cost G=C, charge m*(G+D) to build the baby table, S_sort to
sort it, and h*(G+V) for the giant search. Setup and verification/output are
added in full. Preprocessing includes the whole baby table and sorting, with
no reuse credit; the truncation reduces the table and the walk by the same
factor, so no charge is deferred to an uncharged phase.

Peak memory is 2*m*W*32+2^20 bytes. The two arrays include all keys and
indices. The 1 MiB reserve includes code, public constants, fixed
field/point scratch, loop state and stacks; these use bounded kernels on at
most 12 coefficients and 4572-bit exponents. The finite instruction bodies and
scratch fit the reserve; table storage is separate. No table is uncharged
advice. Input/output and all record lifetimes are included. This memory is
astronomical and mandatory for this construction; it is not a practical attack.

## Submitted computation bound

With `m = h = 212506270031153001307260662621714878127`,
`G = C = 78037248`, `D = 128`, `V = 23088`, `L = 128`:

| Component | Word operations |
| --- | ---: |
| setup | 16777216 |
| baby table | 16583431696778618476543191859019909709190874752 |
| merge sort | 13491598071737841746995364948527434183051264 |
| giant search | 16588310840738533749453206563833704282792670672 |
| verification and output | 79926919168 |

Total word operations: **33185234135588890067743393787802141506110293072**.
Preprocessing word operations: **16596923294850356318290187223968437143390703232**.
Peak bytes: **54401605127975168334658729631159009849088**; nonuniform advice
bytes: 0.

The organizer reference price for this target is 78037248 word operations.
Dividing total and preprocessing work by that fixed price and rounding the
logarithm upward to six decimal places gives submitted time_log2
**128.321575**, preprocessing_log2 **127.321949**, and memory_log2_bytes
**135.320774**. The same arithmetic applied to the parent extents
`m = h = 2^128` reproduces the parent's 129.000802, 128.001176 and 136.000001
exactly, which fixes the accounting convention.

## Evidence and comparison

The correctness argument and the cost calculation are analytic; no full-size
execution is claimed. `experiment.py` runs the identical affine secp256k1 code
path - same p, same generator, same addition formulas, same sort key and same
binary-search lookup - at reduced table size, and reports:

- the maximum field-operation counts observed per point addition, checked
  against the C coefficients (1 inverse, 16 multiplications, 32 additions);
- strict monotonicity of the sorted baby keys, which is what makes the
  binary search exact and collision free;
- recovery of scalars inside the covered prefix, including 0, 1, m-1, m and
  m^2-1, each confirmed by an independent [a]P=Q equality;
- failure of the giant walk on a scalar just outside the prefix, which
  exhibits the truncation boundary rather than asserting it;
- the exact full-size coverage and work ledger in integer arithmetic.

The reduced-size runs are evidence about the arithmetic and about the
list-generation and matching logic only. They do not execute production table
sizes and do not establish the production cost bound by measurement; that
bound is the analytic ledger above. Organizer arithmetic fixtures check the
exact point/GT convention. Heuristic premises: none.

Pollard rho and automorphism-aware methods must be compared under this same
cost model before selecting a best-known ECDLP starting score. This package
does not implement those optimizations. Beating this conservative control
alone does not establish novel cryptanalysis or a weakened hardness
assumption.

References: [generic DLP algorithms, HAC section 3.6](https://cacr.uwaterloo.ca/hac/about/chap3.pdf);
[automorphism-aware Pollard rho](https://www.iacr.org/archive/pkc2014/83830154/83830154.pdf).

## Limitations

The success probability is the exploratory threshold itself, so the claim is
weaker than the parent's probability-1 control and must not be read as a
full-coverage attack. The exponential memory is required and is charged. The
gain over the parent is the constant `2/sqrt(0.39)` in group operations, not a
new decomposition of the group.
