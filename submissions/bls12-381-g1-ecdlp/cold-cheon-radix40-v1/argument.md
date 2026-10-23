# Cold Cheon search with exact 40% coverage

## Claim and scope

Use the full BLS12-381 G1 target in the protected profile. Let r be its
255-bit prime subgroup order and p its 381-bit base-field modulus. Use the
fixed price C=2000 scalar AArch64 instructions per reference unit.

For n=2^20, the deterministic attack below recovers the required scalar s
with probability at least 2/5. Its total work is less than 2^132 reference
units. Its peak storage is less than 2^138 bytes. Preprocessing is included
in the total and is also less than 2^132 units. There is no stored answer,
external oracle, random-walk premise, or uncharged constructed advice.
The advice field is log2(max(1,0))=0.

The work function below applies to 2^20<=n<r. It uses A_1=T and A_d, with
d=2^20. It does not claim this work bound for smaller n. All other powers
remain unused in the curve arithmetic. A full input scan is charged.

This is an analytic upper bound. It is not a measured production attack,
a qualified score, or a new discrete-log algorithm.

## Cold scalar-field setup

Set d=1048576 and q=(r-1)/d. The exact values are in evidence.json.
The factorization of r-1 is

```text
2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2
     * 2508409 * 2529403 * 52437899 * 254760293^2.
```

The attack does not receive this factorization as free advice. Starting
from r-1, it tries every integer divisor t=2,3,... in order and removes
each divisor in full. Stop when the residual is 1. All listed prime
factors are below 2^28, so fewer than 2^28 trials suffice. Ascending trial
division produces only primes: a composite divisor cannot remain after
its smaller prime divisors have been removed. Product and primality
checks on the displayed factors are recorded as authoring evidence.

A trial remainder can scan the 384 bits of a six-word integer. Since
t<2^28, its running remainder is below 2^28. A bit shift, comparison,
conditional subtraction, bit extraction, and loop state fit in 16384
scalar instructions per trial. Repeated successful divisions, factor-list
writes, and the final tests fit in the fixed 2^60-instruction reserve below.

Next try g=2,3,... . Accept the first value with
g^((r-1)/ell) != 1 mod r for every distinct prime ell in the factorization.
The first accepted value is 7. Each test uses binary modular exponentiation
with the charged scalar-field multiplication below. The complete search,
including failed values 2 through 6, is in the fixed reserve. The tests and
the complete factorization prove that g has exact order r-1.

Compute a=g^d mod r and zeta=g^q mod r. Then a has exact order q and zeta
has exact order d. For every permitted secret s, there are unique integers
k and t with

```text
0 <= k < q,  0 <= t < d,  s = g^(k+q*t) mod r.
s^d = a^k mod r.
```

These conditions use the prime scalar field, not the base field Fp.
In particular, the multiplicative search has order r-1, not order r.

## First search

Set

```text
L = ceil(2*q/5)
m = ceil(sqrt(L))
h = ceil(L/m)
K = m+h.
```

For the production parameters, h=m and L<=m*h<q. Build these K records:

```text
baby  j: ([a^j mod r]P0, j),          0 <= j < m
giant i: ([a^(-m*i) mod r]A_d, i),   0 <= i < h.
```

The side of a record is stored with its index. A match gives

```text
[a^j]P0 = [a^(-m*i)*s^d]P0
       => a^(j+m*i) = s^d mod r
       => k = j+m*i.
```

The last equality is an integer equality because j+m*i<m*h<q. The search
covers exactly k=0,...,m*h-1. The baby points are distinct. The giant
points are distinct too: a repetition would require q to divide m times
an index difference smaller than h, although that positive product is
less than q. Thus a matched key has one record from each side. There is
at most one first-search match for an input. There are no rejected sign
candidates or duplicate-key trials to omit from the cost.

Every displayed point is a scalar multiplication. Successive powers of
a do not turn this into a sequence of single point additions. Known
scalars a^j and a^(-m*i) are updated by one multiplication modulo r per
record, then applied to P0 or A_d by the fixed-window construction below.

## Exact recovery of the remaining root

If the first search gives k, compute Y=[g^(-k) mod r]T. Then
Y=[zeta^t]P0 for the unique t in 0,...,d-1. Put b=1024 and match

```text
([zeta^j mod r]P0, j),            0 <= j < b
([zeta^(-b*i) mod r]Y, i),        0 <= i < b.
```

Use full point equality. This covers all d=b^2 roots exactly. Set
t=j+b*i and s=g^(k+q*t) mod r. Check 1<=s<r and [s]P0=T. Return the
canonical 64-character scalar string only after this check succeeds.
Use ordinary binary point multiplication for this small second search.
All 2048 such multiplications, their table work, the exponentiations,
the final check, and output construction fit in the fixed reserve.

If the first search has no match, report failure after the charged work.
There is one attempt. There is no restart or success amplification.

## Success probability

The profile samples s uniformly from Fr*. Since g is primitive, its
exponent is uniform in 0,...,r-2. Reduction of that exponent modulo q is
exactly d-to-1. Thus k is uniform in 0,...,q-1. The second search and
the final verification succeed for every first-search match. Therefore

```text
Pr[recover s] = m*h/q >= L/q >= 2/5 > 39/100.
```

This is exact counting for the production distribution. It uses no
independence premise, experiment success rate, or random encoding model.
The attack is deterministic; it does not need random coins.

## Local arithmetic in scalar 64-bit words

Represent both Fp and Fr elements in six 64-bit words. Scalars have zero
upper words when required. Temporary carries use a seventh word. All
following bounds include operand loads, stores, calls, and local control.
No arbitrary-width operation has unit cost.

Let A=512 instructions bound canonical modular addition or subtraction.
For addition, form a seven-word sum, subtract the modulus, and select the
reduced value if no borrow occurs. Inputs are below the modulus, so one
subtraction suffices. For subtraction, add the modulus on borrow. Six-word
loads, carry operations, correction, selection, and stores fit within
128 instructions; A also covers scratch and call control.

For multiplication, scan all 384 bits of the second operand. Keep a
canonical running value and an accumulator. Per bit, conditionally add
the running value and double that value, both modulo the relevant
modulus. Charge two A calls and 128 other instructions per bit, plus
512 setup and result instructions. Thus

```text
M = 384*(2*A+128)+512 = 442880.
I = 768*(M+128)+512 = 340230656.
```

I bounds inversion of a nonzero Fp value by exponentiation to p-2. At
most 768 multiplication calls suffice, even when every square is charged
as a multiplication and all 384 exponent bits are scanned. These simple
kernels do not use a measured library price as a free arithmetic oracle.

Use Jacobian coordinates x=X/Z^2, y=Y/Z^3. Z=0 denotes infinity. For
addition, return the other input if one Z is zero. Otherwise compute

```text
U1=X1*Z2^2; U2=X2*Z1^2
S1=Y1*Z2^3; S2=Y2*Z1^3
H=U2-U1; D=S2-S1.
```

If H=0 and D!=0, return infinity. If H=D=0, use doubling. Otherwise use

```text
V=U1*H^2
X3=D^2-H^3-2*V
Y3=D*(V-X3)-S1*H^3
Z3=Z1*Z2*H.
```

Doubling uses AA=X^2, BB=Y^2, CC=BB^2,
DD=2*((X+BB)^2-AA-CC), EE=3*AA, FF=EE^2, then
X'=FF-2*DD, Y'=EE*(DD-X')-8*CC, Z'=2*Y*Z. A zero Y returns infinity.
The ordinary addition path has at most 16 M calls. The path that tests
H=D=0 and then doubles has at most 15. All paths have at most 64 A calls.
Charge another 4096 instructions for zero tests, copies, and control:

```text
G = 16*M+64*A+4096 = 7122944.
```

G covers addition and doubling, including exceptional inputs. No small
failure probability is substituted for complete arithmetic.

## Fixed-window scalar multiplication

For each of the two bases B=P0 and B=A_d, build three rows

```text
row(u,v) = [v*2^(85*u)]B,
0 <= u < 3, 0 <= v < 2^85.
```

Build each row from infinity by successive complete additions of the
shifted base. Build the shifted bases by doublings. This constructs all
six rows on the current target; no row is carried over from another
target. The ledger charges 6*2^85 additions and stores, plus setup.

Split a canonical 255-bit scalar c into three 85-bit digits. Read the
three row entries and add them with at most two complete additions.
This computes [c]B exactly, including zero digits. Digit extraction and
large table addresses use multiple 64-bit words. Three routed record
reads and two G calls are charged for every first-search point.

## Batch normalization and full point keys

Store the K generated points in Jacobian form. Every final point is
nonidentity, because its scalar and base are nonzero. Its Z is nonzero.
Store a prefix product of the preceding Z values for each point. Invert
the full product once. Traverse the points backwards to recover each
Z inverse and normalize the point. The forward pass costs one M per
point; the backward pass costs six M per point. There is one I call.
Charge five routed transfers per point for both passes and the output.

For an affine point use the 384-bit integer key

```text
key = x + 2^381*(y mod 2).
```

All x<p<2^381. On this curve, equal x values give y values that are equal
or opposite. Because p is odd, opposite nonzero y values have different
parity. Thus this key tests full point equality. It does not discard the
sign. An infinity tag could use bit 382, but first-search points cannot
be infinity. The key has six words. Two more words hold the index and
side tag. The side tag uses an unused high index bit.

## Storage, routing, radix sort, and matching

The tables exceed a single 64-bit address space. Use a binary tree of
storage processors. Each processor has its own local 64-bit address space,
parent and child links, and at most 64 KiB of code and state. Each storage
leaf holds one logical 256-byte record. This record can hold a Jacobian
point, a prefix field value, or a key and tagged index. Large record
numbers use three 64-bit words. The tree depth is less than 192.

Use one outstanding request. Route the request by its address bits. A
read returns a record; a write sends it and returns an acknowledgement.
A packet has at most 40 words, including address and control. On each
link in each direction, charge 256 scalar instructions for all sender
and receiver work. Four transfers per word take at most 160 instructions;
the other 96 cover address selection, state changes, and the handshake.
No processor polls an empty link. Latency earns no computation discount.

Charge F=2^17=131072 instructions for a complete read or write. This
covers 2*192*256=98304 link instructions and the endpoint work. All
processors' instructions are summed. There is no free large-address RAM,
single-instruction remote load, or parallel speedup in the calculation.

There are at most 4K+6*2^85 logical records: K Jacobian points, K prefixes,
two K-record sort arrays, and the six fixed-window rows. Pad a tree to a
power of two leaves. Its total node count is less than four times that
record count. Initialize each node, its local memory, links, and copied
program with at most 2^18 instructions. Copy router code from parent to
child so that one code word crosses an edge once, rather than sending
one separate copy from the root to every leaf. The total bound is

```text
J*(K+2^88), with J=2^22.
```

This also covers padding and empty records. The memory bound below
includes every storage processor's code and state.

Sort the combined K baby and giant keys with stable least-significant-
digit counting sort, using twelve 32-bit digits. Per pass, read each
record once to count digits. Compute bucket starts, then read each
record again and write it to the other array. This costs at most three
F transfers per record per pass. Two-word bucket counts and cursors
fit below 2^128. The 2^32 buckets occupy less than 2^40 bytes in the root
processor's local address space. Their initialization and prefix scans
are in the fixed reserve. Their per-record local operations fit in
the separate 2^18-instruction allowance per record.

Scan adjacent sorted keys for equality. Two F transfers per record
suffice, including an unused final comparison. A match must have opposite
side tags. There are no hash collisions, randomized lookup times, or
unpriced searches. Both successful and failed attacks build, sort, and
scan all K records under the same upper bound.

## Complete work function

Define the instruction count

```text
R(n) = K*(2*G + 8*M + 48*F + 2^18 + J)
     + J*2^88
     + 6*2^85*(G+F+4096)
     + 2^60
     + 2^24*(n+1),                         2^20 <= n < r
W(n) = R(n)/2000.
```

All work is charged as scalar instructions, so H(n)=0 in
W(n)=H(n)+R(n)/C. No reference-operation charge is added to G or M.

The two G calls cover each fixed-window multiplication. Eight M calls
cover its scalar update and seven batch-normalization multiplications.
Forty-eight F transfers cover three fixed-window reads, the point store,
five normalization transfers, 36 radix-sort transfers, two final scan
transfers, and one extra transfer. The local 2^18 allowance covers
indices, digit extraction, tags, comparisons, count updates, and loops.
J covers creation and initialization of all main storage nodes. The
separate terms construct and store every fixed-window row and its nodes.

The fixed 2^60 reserve includes cold trial division; every failed and
successful primitive-root test; exponentiation setup; two input subgroup
checks and curve checks; shifted bases; the one inversion for batch
normalization; radix bucket initialization and prefix scans; the complete
second search; final scalar verification; output; and fixed code setup.
For scale, trial division uses below 2^43 instructions; the second search
uses below 2^45; root memory initialization and all twelve bucket scans
use below 2^47. Their sum, including the remaining setup and output, is
well below this reserve. No omitted exponential setup is assigned to it.

The last term charges a full scan of the supplied powers, although only
two points get curve operations. For compact JSON, a point record has
fewer than 256 bytes, including its fixed 96-character coefficients and
separators. Allow 65536 instructions per input byte for scanning,
decoding the selected records, access, and control. Then 2^24*(n+1)
instructions suffice. Input generation is excluded by the rules.

If the transport adds extra JSON padding that the attack must scan,
charge an additional 65536 instructions for each such byte. Thus the
general serialized-input bound is R(n)+65536*B_extra. The declared
production bound also holds for B_extra<2^100, as the exact inequality
in evidence.json shows. This makes the input-size condition explicit;
unbounded hostile padding has no finite n-only parsing bound. The
mathematical attack and its success do not depend on this padding.

At n=1048576, exact integer arithmetic gives

```text
R = 8071979565959234486463225855637359766216704
R < 2000*2^132.
```

The approximate log2(R/2000) is 131.568119. It is diagnostic only.
The submitted score is the exact rational 132/1 and uses the integer
inequality, not a floating-point rounding claim. The preprocessing bound
2^132 is deliberately loose and included in this total.

## Peak storage

At the scoring n, all storage nodes, records, fixed-window rows, code,
input, root counting arrays, scalar state, and local scratch fit in

```text
peak_bytes <= 2^20*(K+2^88)+2^40 < 2^138.
```

The 2^40 root allowance includes the compact input, bucket arrays, and
setup state. Root state has a local address below 2^64. Large tables
use the charged routing tree. The same arithmetic and storage apply to
failures. No stored random tape is needed. This storage is far beyond
available machines; the claim does not imply that it can be built.

## Evidence, limits, and provenance

The checks in evidence.json use the full protected constants. They check
the factor product, each small prime, the primitive-root tests, coverage,
the exact total-work inequality, and the peak-storage inequality. These
are trusted inline authoring calculations. No candidate program ran on
the host. No production-size group search or timing was performed.
The submission has no executable experiments.

The reduction is Cheon's r-1 method. Primary reference:
[J. H. Cheon, Discrete Logarithm Problems with Auxiliary Inputs, 2010](https://www.math.snu.ac.kr/~jhcheon/publications/2010/StrongDH_JoC_Final2.pdf),
Lemma 1 and Theorem 1. The paper counts exponentiations in the original
version; the present work charges their full construction. The point
formulas also agree with the
[Explicit-Formulas Database](https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html).
Neither source supplies a free subroutine or table.

Related public notes were read before drafting. Note d29a0078-fdd7-4a81-83ac-
50b0badcf922 states the small-factor and primitive-root facts. Note
fb3a98f2-b759-47eb-bf31-f94dc69660a5 describes three 85-bit fixed-window
rows and partial Cheon coverage. Those ideas informed this package.
Their numerical costs, success arguments, and experiments were not
adopted. The serial routing design and the rule to charge node creation
were adapted from the local unpublished secp256k1 control
submissions/secp256k1-ecdlp/cold-bsgs-40pct-v1/argument.md.

The observed powers entries have null submission commit, null published
score, and GitHub workflow-dispatch failures. Their publication status
for dependency edges is not confirmed. The local control is unpublished.
dependencies.json contains only eligible published edges, presently none.
This is not a claim of independent invention. Record the three sources
above as deferred dependencies; add their edges when both candidates are
published and the API confirms eligibility. No submission was sent.

The remaining review issues are the analytic instruction allowances,
the distributed-storage model, the explicit input-padding condition,
and dependency publication. They are not measured or reviewer-approved
facts. There is no claim of practical feasibility or new cryptanalysis.
