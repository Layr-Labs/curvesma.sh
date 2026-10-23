# secp256k1-ecdlp: deterministic negation-map baby-step/giant-step

Status: ready for independent review; no qualification or emitted score.
This revision keeps every kernel, ledger constant and accounting convention of
the deterministic generic baseline and changes one thing: the baby table is
keyed by the x-coordinate class {R, -R}, so each stored entry covers two
residues and the giant stride doubles. The selected folder is self-contained.
The normative inputs/outputs are those of `production-v2`; `claim.json` binds
that profile and `classical-group-work-v2` by SHA-256. No restricted challenge
distribution, advice, random-oracle premise, heuristic or quantum operation is
used.

## Construction

Let p and r be the protected field modulus and subgroup order, b=bit_length(r)=256.
Set m=isqrt(floor(r/2)) =
240615969168004511545033772477625056926 (about 2^127.5),
stride n=2m-1=481231938336009023090067544955250113851,
L=ceil(log2 m)=128 merge passes and
h=ceil((r+m-1)/n)=240615969168004511545033772477625056929.
All indices, i*n and the reconstructed scalars fit 257-bit signed words of the
ledger below (s+1 limbs) and are reduced modulo r before verification.

```text
Use P from the protected profile; input Q. Starting at infinity, compute jP for
j=0..m-1 by repeated affine addition and store the record (key(jP), j), where
key(R)=(1,0,0) for infinity and (0,x(R),0) otherwise (y is not stored in the key).
Retain U=(m-1)P. Compute E=U+U+P=(2m-1)P, Dn=-E, gamma=Q+U.
Sort the table by key. For i=0..h-1: binary-search key(gamma); on a match with j,
set c=i*n-(m-1); for a in {(c+j) mod r, (c-j) mod r} verify [a]P=Q and return
the first a that verifies. Otherwise set gamma=gamma+Dn. Reject malformed
inputs; exhaust the fixed loops with failure.
```

## Success

For every valid Q there is a unique x in [0,r) with Q=[x]P. Write
x+(m-1)=i*n+u with 0<=u<=n-1=2m-2, so t=u-(m-1) lies in [-(m-1),m-1]. Because
x+(m-1) < r+m-1 <= h*n, the index satisfies 0<=i<h. At giant step i,
gamma=Q+[(m-1)-i*n]P=[t]P. Put j=|t| in [0,m-1]; then gamma is jP or -jP, and
both have key(jP) (negation is (x,-y), and infinity is its own negative). The
complete table therefore contains a match no later than step i.

Soundness of any match: on the curve y^2=x^3+7, two finite points share an
x-coordinate exactly when they are equal or opposite, and the infinity tag is
disjoint from finite keys. So a match at any step i' with entry j' implies
gamma=[+/-j']P, i.e. Q=[c'+/-j']P with c'=i'*n-(m-1). One of the two candidates
equals x mod r and passes the independent [a]P=Q check, which never accepts a
wrong scalar because P has prime order r. Duplicate keys cannot occur for
1<=j<=m-1<r/2 since jP=+/-j'P would force j=j' (mod r); j=0 is the unique
infinity entry. Infinity input gives x=0 through the j=0 entry.

The algorithm succeeds for every valid input: probability 1, hence >=0.39 over
the uniform challenge distribution. No randomness, retries or amplification.
A finite worst-case work budget is used, not an expected runtime.

## Word-operation ledger

Unchanged from the baseline: s=1, k=256 and

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
C covers every case, including the doubling U+U. Keys have K=3 words and
records W=4 words (the third key word is zero, which only over-charges
comparisons). B=16*K+64=112 per key comparison and D=16*W+64=128 per record
move. These bounds allow all operand loads, stores, index arithmetic,
comparisons and branches, up to 16 operations per word plus 64 fixed work.

Setup reserves 2^24 operations for decoding/constants and the negation of E,
plus 3*C for U+U, (2U)+P and Q+U: setup = 2^24+3C.
Verification/output reserves 6*b*C+2^24: up to two [a]P checks and the [r]Q
subgroup check, each at most 2*b point operations, plus curve/encoding checks,
signed reconstruction c+/-j by shift/add, modular reduction and serialization.

Bottom-up merge sort with two m-record buffers and L=ceil(log2 m)=128 passes;
each pass has at most m comparisons, 2m record moves and 128m loop work, so
S=m*L*(B+2*D+128)+4096*L. Binary search over m sorted records performs at most
L+1 comparisons; V=(L+1)*(B+64)+D+256 covers lookup and the returned record.
No constant-time hashing/table oracle is assumed.

With G=C, charge m*(G+D) for the baby table, S for sorting, h*(G+V) for the
full giant search (every step pays one lookup and one addition, including the
last), plus setup and verification in full. Preprocessing is the baby table and
sort with no reuse credit.

Peak memory is 2*m*W*32+2^20 bytes: the two sort buffers of m records plus the
same 1 MiB reserve for code, constants, scratch, loop state and stacks. The
memory is astronomical and mandatory; this is not a practical attack.

## Why the work drops

The baseline pays m+h point operations with m=2^128 and h=ceil(r/2^128), about
2*sqrt(r). Here each table entry answers for both t and -t, so a giant step
advances by n=2m-1 instead of m, and the balanced choice m ~ sqrt(r/2) gives
m+h ~ 2*sqrt(r/2) = sqrt(2r) point operations: a factor sqrt(2) less work and
half the table. Lookups, sorting and moves remain a ~1e-3 relative overhead
against C, as in the baseline. This is the standard negation-map speed-up of
generic algorithms; the secp256k1 order-3 endomorphism is not used here.

## Evidence and comparison

The correctness argument and cost calculation are analytic; no full-size
execution is claimed. `experiment.py` exhaustively runs the same signed-window
indexing on the additive group Z/nZ with the class {v,-v} standing in for the
x-coordinate key; it checks coverage, the last block and both-sign recovery
only at toy scale and does not certify production cost. Heuristic premises: none.

References: [generic DLP algorithms, HAC section 3.6](https://cacr.uwaterloo.ca/hac/about/chap3.pdf);
[negation map and automorphisms in generic ECDLP, Bos-Costello-Miele (PKC 2014)](https://www.iacr.org/archive/pkc2014/83830154/83830154.pdf).
The construction and accounting are complete above.

## Submitted computation bound

| Component | Word operations |
| --- | ---: |
| setup | 250888960 |
| baby table | 18777038857567975237136141435534878154355666176 |
| merge sort | 15276226650538270428971104147059459614642176 |
| giant search | 18782563400220072620721215410950964425896868144 |
| verification and output | 119881990144 |

Total word operations: **37574878484438586128286327950632902160000055600**.
Preprocessing word operations: **18792315084218513507565112539681937613970308352**.
Peak bytes: **61597688107009154955528645754272015621632**; nonuniform advice bytes: 0.

The organizer reference price for this target is 78037248 word operations.
Dividing total and preprocessing work by that price gives time_log2
**128.500802** (exact value 128.5008012...) and preprocessing_log2
**127.501176**, rounded upward. Memory_log2_bytes is **135.500001**. All table
access, sorting, failure paths, verification and output costs are included.
