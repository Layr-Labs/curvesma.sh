# Indexed baby-step Cheon for the supplied \(2^{20}\) power

The scalar field order is the prime

\[
r = 52435875175126190479447740508185965837690552500527637822603658699938581184513.
\]

\(r\) is odd and has bit length 255, so \(r \ge 2^{254}+1\) and \(r-1 \ge 2^{254}\). Let \(n = 2^{20}\) and \(m = (r-1)/n\). Then \(2^{20}\) divides \(r-1\) and

\[
m = 50006747412801924209068050869165387952509453297164571592906626415194112,
\]

so \(m \ge 2^{234}\). The public input at this \(n\) contains \(T = [s]P_0\) and \(A_n = [s^n]P_0\) in the order-\(r\) subgroup \(G_1 = \langle P_0 \rangle\) of \(y^2 = x^3 + 4\) over \(\mathbb{F}_p\). The algorithm reads those two points and does not need the intermediate powers. Generation of the powers is not charged. Scanning and canonical decoding of the whole supplied sequence is charged below.

## Factorisation and a primitive root

\[
\begin{align*}
r-1 = 2^{32} &\cdot 3 \cdot 11 \cdot 19 \cdot 10177 \cdot 125527 \cdot 859267 \cdot 906349^{2}\\
&\cdot 2508409 \cdot 2529403 \cdot 52437899 \cdot 254760293^{2}.
\end{align*}
\]

Every prime on this list is at most \(254760293\). Each is proved prime by trial division up to its square root. The product of the multiset equals \(r-1\). Both checks are finite and are repeated by `experiment.py`. For every prime \(q\) on the list, \(7^{(r-1)/q} \not\equiv 1 \pmod r\), so \(7\) is a primitive root modulo \(r\). No heuristic is used.

Let \(\omega = 7^{2^{20}} \bmod r\) and \(\kappa = 7^{m} \bmod r\). Then \(\omega\) has order \(m\) and \(\kappa\) has order \(2^{20}\). The unique subgroup of \(m\)-th roots of unity in \(\mathbb{F}_r^\times\) is \(\langle \omega \rangle\), and the unique subgroup of \(2^{20}\)-th roots of unity is \(\langle \kappa \rangle\).

## Reduction

For every \(s \in \{1,\ldots,r-1\}\), \(s^{r-1} \equiv 1\), so \((s^{2^{20}})^m \equiv 1\). Hence there is a unique \(\lambda \in \{0,\ldots,m-1\}\) with

\[
s^{2^{20}} \equiv \omega^\lambda \pmod r.
\]

The supplied point is \(A_n = [s^{2^{20}}]P_0 = [\omega^\lambda]P_0\). The kernel of the map \(x \mapsto x^{2^{20}}\) on \(\mathbb{F}_r^\times\) is exactly \(\langle \kappa \rangle\), so the full preimage is the coset

\[
s = 7^\lambda \cdot \kappa^j \bmod r, \qquad j \in \{0,\ldots,2^{20}-1\}.
\]

One such \(j\) satisfies \([s]P_0 = T\), and it is unique: \([a]P_0 = [b]P_0\) if and only if \(a \equiv b \pmod r\), and \(\kappa\) has order \(2^{20}\).

## Baby index that survives sorting

Let

\[
N = 223621884914696853578958571250384896.
\]

This integer satisfies \((N-1)^2 < m \le N^2\), so \(N = \lceil \sqrt{m} \rceil\). Direct comparison of decimals gives

\[
2^{117} = 166153499473114484112975882535043072 < N
\]

and

\[
N < 332306998946228968225951765070086144 = 2^{118},
\]

because all three integers have 36 digits and the leading triples satisfy \(166 < 223 < 332\). Thus \(\lceil \log_2 N \rceil = 118\). Also

\[
N < 228461061775532415655341838485684224 = 11 \cdot 2^{114},
\]

because both have 36 digits and \(223 < 228\).

Write a candidate exponent as \(j + N \ell\) with \(0 \le j < N\) and \(0 \le \ell < N\). Since \(m \le N^2\), every residue modulo \(m\) has at least one such representative, and the canonical \(\lambda\) has \(j = \lambda \bmod N\) and \(\ell = \lfloor \lambda / N \rfloor \le N-1\).

The exact difference \(N^2 - m = 218380746987008646073201729736704\) has 33 digits, so it is strictly smaller than \(N\). Therefore the giant exponents \(0, N, \ldots, (N-1)N\) are strictly increasing and all lie in \(\{0,\ldots,m-1\}\). They are distinct modulo \(m\). Baby exponents in a window of length \(N\) are likewise distinct modulo \(m\).

Define the baby points and the giant points by

\[
R_j = [\omega^{-j}] A_n, \qquad G_\ell = [\omega^{N\ell}] P_0.
\]

Then \(R_j = [\omega^{\lambda - j}] P_0\). A match \(R_j = G_\ell\) means \(\lambda \equiv j + N\ell \pmod m\). Distinct baby indices give distinct points and distinct giant indices give distinct points, so no record is duplicated, and the canonical pair always matches. The search uses the first match it finds. Since \(j + N\ell \le N^2 - 1 < m + N < 2m\), any match satisfies \(j + N\ell \in \{\lambda,\ \lambda + m\}\). Either value is a valid logarithm of \(A_n\) in the base \(7\): \((7^{j+N\ell})^{2^{20}} = \omega^{j+N\ell} = \omega^{\lambda} = s^{2^{20}}\), so \(s \cdot 7^{-(j+N\ell)}\) lies in the kernel \(\langle\kappa\rangle\). Every match therefore returns the same coset \(7^{\lambda}\langle\kappa\rangle\), and the second stage below searches that whole coset.

The baby table is an array of \(N\) records of 104 bytes:

| offset | contents |
|---:|---|
| 0..47 | canonical affine \(x \in \{0,\ldots,p-1\}\), 48-byte big-endian |
| 48..95 | canonical affine \(y\), same encoding |
| 96..103 | baby index \(j\), little-endian integer |

The sort key is the 96-byte concatenation \((x,y)\). Every swap exchanges the whole 104-byte record, so the integer at offset 96 moves with its point. After the sort, the position of a record is the rank of its key and is not equal to \(j\) in general. Binary search compares keys only. On a hit it returns the integer stored at offset 96 of the matched record. That integer is the baby index \(j\) used in \(\lambda = j + N\ell\). Dropping it and using the rank would return a different integer and would not be a discrete logarithm.

The same layout, with the coset index in place of \(j\), is used for the order-\(2^{20}\) search below.

## Fixed-base evaluation

Every baby scalar \(\omega^{-j} \bmod r\) and every giant scalar \(\omega^{N\ell} \bmod r\) lies in \(\{1,\ldots,r-1\}\) and therefore has a 255-bit representation

\[
k = d_0 + d_1 2^{85} + d_2 2^{170}, \qquad 0 \le d_i < 2^{85}.
\]

For a fixed base \(B \in \{A_n, P_0\}\) the implementation precomputes, for each column \(i \in \{0,1,2\}\) and each digit \(u \in \{1,\ldots,2^{85}-1\}\),

\[
T_i[u] = [u \cdot 2^{85 i}] B
\]

in canonical affine coordinates. Column bases advance by 85 Jacobian doublings. Within a column the points are successive mixed additions of the normalised column base, starting from that base. Digit 0 is not stored. Building at most \(2^{85}\) additions per column is an upper bound on the \(2^{85}-2\) loop iterations.

A scalar multiplication reads at most one affine point from each of the three columns and combines them with the scheduled mixed addition below. The identity is a flag, not a curve point. The first nonzero digit installs that affine point as a Jacobian point with \(Z = 1\). Each later digit performs one scheduled addition. At most two additions occur; the cost bound charges three.

Charging a scalar multiplication as one group operation is not used. The order-6 automorphism is not used to shrink the search.

## Field operations and the scheduled addition

\(C = 2000\) instruction equivalents is one reference group operation. The measured `blst_fp_mul` body is 383 instructions. A squaring is implemented by that multiplication, so it costs at most 383 instructions as well. The envelope below charges **512** instructions for every base-field multiplication or squaring, which is strictly above 383. A base-field addition or subtraction of two reduced 6-limb values is six loads, six loads, an `adds`/`adcs` chain, a conditional subtraction of \(p\), and six stores. That program is under 80 instructions; the envelope charges **128**. A 4-limb Montgomery product modulo \(r\) is under 400 instructions; the envelope charges **2048** for every \(\mathbb{F}_r\) multiplication. No reference addition \(H\) is billed on top of these instructions. Total work is \(W = R/2000\).

Jacobian doubling, \(a = 0\), the straight-line program `dbl-2009-l`:

```
A = X^2
B = Y^2
C = B^2
t1 = X+B
t2 = t1^2
t3 = t2-A
t4 = t3-C
D = t4+t4
e1 = A+A
E = e1+A
F = E^2
d2 = D+D
X3 = F-d2
t5 = D-X3
t6 = E*t5
c2 = C+C
c4 = c2+c2
c8 = c4+c4
Y3 = t6-c8
yz = Y*Z
Z3 = yz+yz
```

This is \(2M+5S+14A\). At 383 and 128 instructions plus 256 for copies and control, the derived cost is \(7\cdot 383 + 14\cdot 128 + 256 = 4729\). The envelope charges \(8\cdot 512 + 20\cdot 128 + 1024 = 7680 > 4729\). Every point doubled in this algorithm lies in the odd-order subgroup \(\langle P_0 \rangle\) and is not infinity, so \(Y \ne 0\) and the doubling formula applies.

Mixed addition of a Jacobian point and an affine point is not run blindly. First compare coordinates:

```
Z2 = Z1^2
Z3 = Z2*Z1
X2Z = X2*Z2
Y2Z = Y2*Z3
```

which is \(3M+1S\), then two integer comparisons. If the \(x\)-coordinates differ, run `madd-2007-bl`:

```
Z1Z1 = Z1^2
U2 = X2*Z1Z1
Z1Z1Z1 = Z1*Z1Z1
S2 = Y2*Z1Z1Z1
H = U2-X1
HH = H^2
h2 = HH+HH
I = h2+h2
J = H*I
s2y = S2-Y1
r = s2y+s2y
V = X1*I
r2 = r^2
r2j = r2-J
v2 = V+V
X3 = r2j-v2
vx = V-X3
rv = r*vx
yj = Y1*J
yj2 = yj+yj
Y3 = rv-yj2
zh = Z1+H
zh2 = zh^2
zt = zh2-Z1Z1
Z3 = zt-HH
```

This is \(7M+4S+14A\). The test plus this addition is \(10M+5S+14A\) plus two comparisons. At 383 and 128, plus 128 for the comparisons and 512 for control, loads of the 96-byte affine point, and address arithmetic, the derived cost is \(15\cdot 383 + 14\cdot 128 + 640 = 8177\). If the test says the points are equal, the implementation doubles instead. If it says they are inverses, the result is infinity and the identity flag is set. Both exceptional branches are cheaper than the distinct-point branch. The envelope charges the distinct-point branch on every scheduled addition, at

\[
16\cdot 512 + 24\cdot 128 + 2048 = 13312 > 8177
\]

instructions. The spare covers more than ten extra multiplications if a reduction or a copy was omitted from the listing. Call this a scheduled addition.

Normalisation of a non-infinity Jacobian point uses one batched inversion of the \(Z\) coordinates (Montgomery prefix products: fewer than three multiplications per point, plus one Fermat inversion per batch) followed by

```
z2 = zinv^2
z3 = z2*zinv
x = X*z2
y = Y*z3
```

and one multiplication of each coordinate by the public Montgomery constant to land in \(\{0,\ldots,p-1\}\). That is at most \(8M+1S\) plus copies. The envelope charges \(12\cdot 512 + 1024 = 7168\) instructions per normalised point, which dominates that program. Fermat inversion of one element, \(a^{p-2}\), is at most 380 squarings and 228 multiplications. The envelope charges \(611\cdot 512 + 4096\) per inversion. Eight inversions cover the six precomputation columns and two spares. Four further inversions cover the baby batch, the giant batch, and two spares. Zeros do not occur: every scalar in the two searches is in \(\{1,\ldots,r-1\}\) and the base points have order \(r\).

A key comparison loads at most twelve 64-bit limbs and branches. The envelope charges 256 instructions. A record move copies 104 bytes. The envelope charges 256 instructions. Bottom-up mergesort has \(\lceil \log_2 N \rceil = 118\) passes; each pass compares at most \(N-1\) pairs and moves each record at most once. Binary search on \(N\) keys uses at most 118 comparisons; the bound charges 120.

## Second stage

Compute the scalar \(7^\lambda \bmod r\) by square-and-multiply and the point \(U = [7^\lambda]P_0\). Let \(B = 2^{10}\). Search \(j = j_0 + B j_1\) with \(0 \le j_0, j_1 < B\) such that

\[
[\kappa^{j_0}] U = [\kappa^{-B j_1}] T.
\]

The baby records store \(j_0\) at offset 96 exactly as above. There are \(B^2 = 2^{20}\) labels and exactly one solution. Each successive baby or giant step is one variable-base scalar multiplication by the fixed scalar \(\kappa\) or \(\kappa^{-B}\): at most 255 doublings and 255 scheduled additions. The bound charges 256 of each, and charges 4096 such scalar multiplications. The real number is \(2\cdot 1024\) chain steps plus \(U\), the final check \([s]P_0 = T\), two subgroup checks \([r]T = \infty\) and \([r]A_n = \infty\), and a margin, altogether well under 4096.

Sixty-four \(\mathbb{F}_r\) exponentiations of at most 512 multiplications each cover \(\omega\), \(\omega^N\), \(\omega^{-1}\), \(\kappa\), \(\kappa^{-1}\), \(\kappa^{-B}\), \(7^\lambda\), the twelve primitive-root powerings, and a margin. The bound charges those 64 exponentiations at 2048 instructions per multiplication.

## Instruction total

Let \(Q = 2N\) be the number of baby scalars plus giant scalars. With \(G = \lceil m/N \rceil \le N\), one has \(Q \le 2N\). The \(N\)-linear envelopes are

| item | instructions |
|---|---:|
| 3 scheduled additions on each of the \(Q\) scalars | \(79872\, N\) |
| normalisation of those \(Q\) points | \(14336\, N\) |
| one \(\mathbb{F}_r\) update per scalar | \(4096\, N\) |
| three digit extractions per scalar, 512 each | \(3072\, N\) |
| writing the baby index | \(256\, N\) |
| 118 sort comparisons, 256 each | \(30208\, N\) |
| 118 sort moves, 256 each | \(30208\, N\) |
| 120 binary-search comparisons on \(\le N\) giants | \(30720\, N\) |

The coefficients sum as \(79872+30208 = 110080\), \(+30208 = 140288\), \(+30720 = 171008\), \(+14336 = 185344\), \(+4096 = 189440\), \(+3072 = 192512\), \(+256 = 192768\). The linear part is therefore at most \(192768\, N\). Since \(N < 11\cdot 2^{114}\),

\[
192768 \cdot 11 \cdot 2^{114} = 2120448 \cdot 2^{114}.
\]

And \(2120448 = 2070\cdot 1024 + 768\), so

\[
2120448 \cdot 2^{114} = 2070 \cdot 2^{124} + 768 \cdot 2^{114} = 2070 \cdot 2^{124} + 3 \cdot 2^{122},
\]

where \(768 = 3\cdot 2^8\).

Every remaining charged instruction fits in the following caps: precomputation additions \(6\cdot 13312\cdot 2^{85} < 2^{102}\); precomputation normalisation \(< 2^{101}\); 512 doublings \(< 2^{22}\); eight Fermat inversions \(< 2^{22}\); four further inversions \(< 2^{22}\); 4096 variable-base scalar multiplications \(< 2^{35}\); 64 field exponentiations \(< 2^{26}\); scanning \(2^{20}\) input points at 4096 instructions each, exactly \(2^{32}\); and a \(2^{48}\) absolute pad. The sum of these caps is \(< 2^{102}+2^{101}+2^{49} < 2^{103} < 2^{104}\).

Hence the total instruction count \(R\) satisfies

\[
R < 2070\cdot 2^{124} + 3\cdot 2^{122} + 2^{104}.
\]

Now \(3\cdot 2^{122} = (3/4)2^{124}\) and \(2^{104} < (1/4)2^{124}\), so \(R < 2071 \cdot 2^{124}\). Therefore

\[
W = R/2000 < \frac{2071}{2000} \cdot 2^{124}.
\]

And \(2071\cdot 100 = 207100 < 214000 = 107\cdot 2000\), so \(2071/2000 < 107/100\). Finally \(107^{10} = 196715135728956532249 < 200000000000000000000 = 2\cdot 100^{10}\), hence \((107/100)^{10} < 2\) and \(107/100 < 2^{1/10}\). Combining these,

\[
W(2^{20}) < 2^{124 + 1/10} = 2^{1241/10} = 2^{124.1}.
\]

The same split shows that precomputation, the field exponentiations and the input scan are together \(< 2^{103}\) instructions, hence \(< 2^{93}\) group operations. That preprocessing is inside \(W\), not an extra objective. Nonuniform advice is empty: \(\omega\), \(\kappa\) and the window tables are computed from public parameters and from \(A_n, P_0\), and their construction is inside the envelopes above.

## Memory

Simultaneous storage of the six window columns (\(6 \cdot 2^{85}\) affine points, 96 bytes each), a Jacobian image of every baby and every giant (144 bytes each), the baby records and a mergesort buffer (104 bytes each), and \(2^{20}\) bytes of code is at most

\[
99\cdot 2^{119} + 143\cdot 2^{118} + 2^{95} < 342 \cdot 2^{118}.
\]

The first term bounds \(Q < 11\cdot 2^{115}\) Jacobian points at 144 bytes, because \(2N < 22\cdot 2^{114}\). The second bounds two record arrays. Since \(342 < 384\),

\[
342 \cdot 2^{118} < 384 \cdot 2^{118} = 48 \cdot 2^{121} = \tfrac{3}{2} \cdot 2^{126}.
\]

And \((3/2)^5 = 243/32 < 8 = 2^3\), so \(3/2 < 2^{3/5}\) and the peak is \(< 2^{126 + 3/5} = 2^{126.6}\) bytes. Giants may be streamed; counting them all only enlarges the bound.

## Success

Every step is deterministic. There are no algorithm coins. For every \(s \in \{1,\ldots,r-1\}\) a matching baby/giant pair exists, the record at offset 96 of the matched record contains its baby index, and whichever match is used, \(j + N\ell \equiv \lambda \pmod m\) returns the same coset \(7^{\lambda}\langle\kappa\rangle\). The coset search then returns the unique \(j\) with \(s = 7^{\lambda}\kappa^{j}\), and \([s]P_0 = T\) by construction. The final scalar multiplication checks that equality before output. Success probability is 1, which is at least 0.39. The bound is worst case, not expected: every giant is looked up, and every scheduled addition is charged at the distinct-point envelope.

The bound applies only at \(n = 2^{20}\). That is the scoring point, \(2^{20}\) divides \(r-1\), and \(A_n\) is the supplied power the algorithm uses.

`experiment.py` repeats the factorisation, primality, primitive-root, cover, coefficient and envelope inequalities, and it recovers every secret in the multiplicative group modulo 233 by reading the stored baby index rather than the sorted rank. That run is a finite check of this argument. It is not an execution of the production curve and it is not a measurement of the 383-instruction multiplication.
