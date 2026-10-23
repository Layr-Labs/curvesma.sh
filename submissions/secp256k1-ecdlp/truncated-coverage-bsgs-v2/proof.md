# secp256k1-ecdlp: truncated-coverage deterministic BSGS with a rebalanced table

Status: ready for independent review; no qualification or emitted score.
This is a conservative accounting control, not a best-known attack claim.
The normative inputs/outputs are those of `production-v2`; `claim.json`
binds that profile and `classical-group-work-v2` by SHA-256. No restricted
challenge distribution, advice, random-oracle premise or quantum operation is used.

## Construction and success

Let p and r be the protected field modulus and subgroup order, b=bit_length(r)=256,
and L=128. The baby table holds m records for the fixed integer

```text
m = 215157619019306985972612376563775820358   (2^127 < m < 2^128, L=ceil(log2 m)=128)
```

and the giant walk performs H steps for

```text
H = ceil(2*r/(5*m)) = 215269326301525377570715387988135775797.
```

All table indices/addresses and i*m+j fit 256 bits: m<2^128, H<2^128 and
H*m<r<2^256. Both integers are fixed constants of the construction, not
inputs, coins or advice.

```text
Use P from the protected profile; input Q. Starting at infinity, store (jP,j)
for j=0..m-1 by repeated affine addition, retaining mP after the final update.
Set D=-(mP) and gamma=Q. Sort the table by (infinity_tag,x,y), tag 1 for infinity and 0 otherwise; use zero
coordinate words for infinity. For i=0..H-1, binary-search gamma; on a match
with j, set a=(i*m+j) mod r, verify [a]P=Q, and return a. Otherwise set
gamma=gamma+D. Reject malformed inputs; exhaust the fixed loop and return no
output when no match occurs.
```

For every valid Q there is a unique x in [0,r) with Q=[x]P. Write x=i*m+j with
0<=j<m and 0<=i<ceil(r/m). At giant step i, gamma=[x-i*m]P=[j]P, so the table
contains a match whenever 0<=i<H, i.e. whenever x<H*m. Any match yields
[a]P=Q after reduction modulo r; the independent attack-side equality check
accepts it. Infinity gives x=0 and is represented explicitly.

The covered set is exactly the interval [0,H*m) of scalars. Because
H=ceil(2*r/(5*m)), we have 5*H*m-2*r = 642482851172190102134522932175858387956 > 0,
hence H*m >= 2*r/5 and the success probability over uniform x in [0,r) is
H*m/r = 0.4000000000000000000000000000000000000011... >= 2/5 = 0.4, which is
above the 0.39 acceptance floor. Also H*m<r, so the walk is genuinely
truncated: the last 3*r/5 of the scalar interval is not covered and is not
claimed. For inputs outside the covered interval the fixed loop runs all H
steps and fails; that worst-case work is included in the bound below, so
failures are charged, not amortized. There is no randomness, retry,
restart or amplification, and the bound is a finite worst-case work budget
rather than an expected-runtime exponent standing in for total work.

The table size is the analytic minimizer of the charged bound under the
coverage constraint. With L fixed by m in (2^127,2^128], the per-record cost is
a=C+D+496*L and the per-giant-step cost is b=C+V, so the bound
a*m + b*H is minimized at m=floor(sqrt(2*r*b/(5*a))) = floor(sqrt(2*r*(C+V)/(5*(C+D+496*L))))
= 215157619019306985972612376563775820358, which is the value used. The
truncation moves the balance point from m=2^128 (where a full walk costs
2^129.000802) down to m~0.632*2^128 with a 2/5-length walk.

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
C covers every case. Keys have K=3 words and records W=4 words.

Setup reserves 2^24 operations for decoding/constants and the step negation.
Verification/output reserves 4*b*C+2^24: this covers [r]Q and [a]P checks,
curve/encoding checks, scalar reconstruction by shift/add and serialization.

Use bottom-up merge sort with two m-record buffers and exactly L=128 passes;
L=ceil(log2 m), so 128 passes sort any m<2^128 records. A lexicographic key
comparison costs B=16*K+64; moving a record costs D=16*W+64. These bounds
allow all operand loads, stores, index arithmetic, comparisons and branches,
up to 16 operations per word plus 64 fixed work. Each merge pass has at most
m comparisons and 2*m record moves, with 128*m additional loop work. Thus
sorting costs S=m*L*(B+2*D+128)+4096*L. Binary search over m records performs
at most L+1=129 key comparisons; V=(L+1)*(B+64)+D+256 covers lookup and the
returned record. No constant-time hashing/table oracle is assumed.

For group-step cost G=C, charge m*(G+D) to build the baby table, S to sort it,
and H*(G+V) for the truncated giant search; the loop always executes H steps,
so the same bound covers success and failure. Setup and verification/output
are added in full. The accounting below gives these exact integer components
and their sum T. The reported logarithms are rounded up to six decimal places;
exact integers remain the audit authority. Preprocessing includes the whole
baby table and sorting, with no reuse credit.

Peak memory is 2*m*W*32+2^20 bytes. The two arrays include all keys and
indices. The 1 MiB reserve includes code, public constants, fixed field/point
scratch, loop state and stacks; these use bounded kernels on at most 12
coefficients and 4572-bit exponents. The finite instruction bodies and
scratch fit the reserve; table storage is separate. No table is uncharged
advice. Input/output and all record lifetimes are included. This memory is
astronomical and mandatory for this construction; it is not a practical attack.

## Evidence and comparison

The correctness argument and cost calculation are analytic; no full-size execution is
claimed. Organizer arithmetic fixtures check the exact point/GT convention;
organizer exhaustive finite-group checks exercise table coverage, last-block
handling and modular output recovery. The accompanying toy experiment checks
the truncated coverage set exactly on a small additive group: every covered
scalar is recovered, every uncovered scalar exhausts the fixed loop, and the
coverage fraction is at least 2/5. It is evidence only at that scope and does
not certify production cost, scientific approval or an executed production
break. Heuristic premises: none.

The submitted bound is below the full-coverage generic control at
time_log2 129.000802 because it stops the giant walk at the 2/5 coverage floor
and rebalances the table against the walk; it does not change the cost model,
the primitive prices, the success floor or the accounting of preprocessing,
failures, sorting, lookup, verification and output. Pollard rho and
automorphism-aware methods must be compared under this same cost model before
selecting a best-known ECDLP starting score. This package does not implement
those optimizations. Beating this deliberately conservative control alone does
not establish novel cryptanalysis or a weakened hardness assumption.

References: [generic DLP algorithms, HAC section 3.6](https://cacr.uwaterloo.ca/hac/about/chap3.pdf);
[automorphism-aware Pollard rho](https://www.iacr.org/archive/pkc2014/83830154/83830154.pdf).
The needed construction and accounting are given above; review need not
fetch citations to reconstruct them.

## Submitted computation bound

The following integer word-operation bounds follow from the implementation formulas above. They are ordinary proof evidence, not a required submission format.

| Component | Word operations |
| --- | ---: |
| setup | 16777216 |
| baby table | 16790336014674410523771477732160961672985700608 |
| merge sort | 13659926916297761925429214563280999283412992 |
| giant search | 16803995941590708285696906946724242672334487792 |
| verification and output | 79926919168 |

Total word operations: **33607991883181416571393813893448485424547297776**. Preprocessing word operations: **16803995941590708285696906946724242672269113600**. Peak bytes: **55080350468942588408988768400326611060224**; nonuniform advice bytes: 0.

The organizer reference price for this target is 78037248 word operations. Dividing total and preprocessing work by that fixed price gives submitted time_log2 **128.339838** and preprocessing_log2 **127.339838**, rounded upward. Memory_log2_bytes remains **135.338662**. All table access, sorting, failure paths, verification and output costs are included; reference-operation normalization changes units, not the work or algorithm.
