# Width-5 wNAF scalar recurrence with blocked Montgomery normalization

Parent definition: BLS12-381 G1, order

\[
r = \mathtt{0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001},
\]

base point \(P_0\), input \(A_i = [s^i]P_0\) for \(i = 1,\ldots,n\) with scoring value \(n = 2^{20}\). The attack reads only \(A_1 = T = [s]P_0\) and \(A_n = [s^n]P_0\). It returns \(s \in \{1,\ldots,r-1\}\).

Reference price: one batched incomplete affine G1 addition costs \(1\), and \(C = 2000\) scalar AArch64 instructions make one such unit. Field multiplications are not reference additions. Each one is charged at the published `blst_fp_mul` count of \(383\) instructions (COST_MODEL.md). Squarings are charged as multiplications. No reference addition is also unpacked into its own field multiplications.

## Parameters

The BLS parameter \(X = \mathtt{0xd201000000010000}\) satisfies \(r - 1 = X^2(X-1)(X+1)\). Exact cofactors:

- \(X = 2^{16} \cdot 906349 \cdot 254760293\)
- \(X - 1 = 19 \cdot 125527 \cdot 2508409 \cdot 2529403\)
- \(X + 1 = 3 \cdot 11 \cdot 10177 \cdot 859267 \cdot 52437899\)

Hence

\[
r - 1 = 2^{32} \cdot 3 \cdot 11 \cdot 19 \cdot 10177 \cdot 125527 \cdot 859267 \cdot 906349^2 \cdot 2508409 \cdot 2529403 \cdot 52437899 \cdot 254760293^2.
\]

Every odd factor is smaller than \(2^{32}\). Deterministic Miller–Rabin with bases \(2, 7, 61\) (decisive for all odd integers below \(4\,759\,123\,141\)) certifies each odd factor prime. The power of two is the trailing zero bits of \(r-1\). In particular \(v_2(r-1) = 32\), so \(d = 2^{20}\) divides \(r-1\).

Set \(m = (r-1)/d\). Then \(2^{234} < m < 2^{235}\). The integer

\[
M = 223621884914696853578958571250384896
\]

satisfies \((M-1)^2 < m \le M^2\), so \(M = \lceil\sqrt{m}\rceil\), and \(M(M-1) < m \le M^2\), so \(\lceil m/M\rceil = M\). Also \(2^{117} < M < 2^{118}\). Thus every baby index below \(M\) fits in \(118\) bits and does not fit in \(64\) bits.

Let \(\zeta = 7\). For every prime \(q\) dividing \(r-1\), \(7^{(r-1)/q} \not\equiv 1 \pmod r\), so \(7\) is a primitive root modulo \(r\). Certification is twelve exponentiations in \(\mathbb{Z}/r\mathbb{Z}\).

## Decomposition

Write the unknown as \(\alpha = s \bmod r\) with \(1 \le \alpha \le r-1\). Because \(7\) generates \((\mathbb{Z}/r\mathbb{Z})^\times\), there is a unique exponent \(e \in \{0,\ldots,r-2\}\) with \(\alpha \equiv 7^e \pmod r\). Divide by \(m\):

\[
e = k_0 + m k_1, \qquad 0 \le k_0 < m,\quad 0 \le k_1 < d.
\]

This is ordinary division. It does not need \(\gcd(d, m) = 1\). Here \(\gcd(d, m) = 2^{12}\), and the mixed radix still covers every residue exactly once. Then

\[
\alpha^d \equiv (7^d)^{k_0} \pmod r, \qquad \alpha \cdot 7^{-k_0} \equiv (7^m)^{k_1} \pmod r,
\]

because \(7^{dm} = 7^{r-1} \equiv 1\) and \(7^d\) has order \(m\) while \(7^m\) has order \(d\).

## Two searches

**Stage \(k_0\).** Write \(k_0 = u + M v\) with \(0 \le u < M\) and \(0 \le v < M\). The bounds on \(M\) put every \(k_0 \in \{0,\ldots,m-1\}\) in exactly one such pair, and \(u + Mv < m\). In the curve group the match is

\[
[\alpha^d \cdot 7^{-d u}] P_0 = [7^{d M v}] P_0.
\]

Left side, \(u = 0\): the input point \(A_d\). The passage \(u \to u+1\) multiplies the scalar by the fixed unit \(\lambda = 7^{-d} \bmod r\). Right side, \(v = 0\): \(P_0\). The passage \(v \to v+1\) multiplies the scalar by the fixed unit \(\rho = 7^{d M} \bmod r\).

Both units lie in \(\mathbb{Z}/r\mathbb{Z}\) and are strictly less than \(r < 2^{255}\), so each is four little-endian \(64\)-bit limbs (the top limb is at most \(63\) bits). Multiplying a point by one of them is a \(4\)-limb scalar recurrence described below.

Baby points are pairwise distinct: \(7^{d(u'-u)} \equiv 1 \pmod r\) with \(|u'-u| < M < m\) forces \(u' = u\), because the order of \(7^d\) is \(m\). Giant points collide only if \(m\) divides \(M(v'-v)\). Let \(g = \gcd(M, m) = 218380746987008646073201729736704\). Then \(m/g = 228988810152649578064853576960394133503 \ge M\), and \(\gcd(M/g, m/g) = 1\), so \(m/g\) divides \(v'-v\). Distinct stage-1 giant indices differ by less than \(M \le m/g\), so the difference is \(0\). The matching pair \((u, v)\) is unique. Set \(k_0 = u + Mv\).

**Stage \(k_1\).** One \(4\)-limb scalar multiplication sends \(T\) to \(T_2 = [7^{-k_0}] T = [\alpha \cdot 7^{-k_0}] P_0\). Let \(M_2 = 2^{10}\) and write \(k_1 = u_2 + M_2 v_2\) with \(0 \le u_2, v_2 < M_2\). The same recurrence, now multiplying by \(7^{-m}\) and by \(7^{m M_2}\), finds the unique pair. Set \(k_1 = u_2 + M_2 v_2\) and

\[
s = 7^{k_0 + m k_1} \bmod r.
\]

One final \(4\)-limb scalar multiplication checks \([s] P_0 = T\). The check passes for every \(s \in \{1,\ldots,r-1\}\). There is no random coin and no rejected trial. Success probability is \(1\).

The same schedule works for every \(n\) with \(1 \le n < r\) and \(n \mid (r-1)\), reading \(A_1\) and \(A_n\), with \(d\) replaced by \(n\) and \(M\) by \(\lceil\sqrt{(r-1)/n}\rceil\). The scored point is \(n = 2^{20}\).

## Records

A baby record is \(112\) bytes:

- bytes \(0..47\): affine \(x \in \mathbb{F}_p\), six little-endian \(64\)-bit limbs;
- bytes \(48..95\): affine \(y\), same layout;
- bytes \(96..111\): the baby index, two little-endian \(64\)-bit limbs.

\(p\) is \(381\) bits, so a field element fits in six limbs and the top limb is shorter than \(64\) bits. The index integer is \(j_0 + j_1 2^{64}\). For every index produced above, \(0 \le j < M < 2^{118} < 2^{128}\), so two limbs are necessary and sufficient. An \(8\)-byte slot cannot hold \(j\).

Mergesort and the later scans move the \(112\)-byte record as one block, so the index stays attached to its point. The sort key is the \(96\)-byte coordinate prefix under unsigned byte order. Search uses that same order. Distinct points have distinct prefixes, so the key determines the stored index.

Nothing in either list is the point at infinity: every scalar that appears is a power of \(7\) times a non-zero field element, hence non-zero modulo \(r\), so every \(Z\) coordinate is non-zero.

## Curve arithmetic inside one scalar step

Points stay in Jacobian coordinates \((X:Y:Z)\) on \(y^2 = x^3 + 4\) during the addition chain. A doubling uses at most \(12\) base-field multiplications (an explicit \(a = 0\) double is \(5\) squares and \(2\) multiplications: \(X^2\), \(Y^2\), \(Y^4\), \((X+Y^2)^2\), \((3X^2)^2\), times \(Y Z\) and the slope–difference product). A Jacobian addition uses at most \(16\) base-field multiplications (an explicit formula is \(4\) squares and \(7\) multiplications). Field additions inside one curve operation are at most \(32\), each charged at \(64\) instructions (six-limb add, carry, and a conditional subtract of \(p\)). Negating a point for a negative digit is one field negation, charged inside the same \(32\)-addition envelope.

Each recurrence step multiplies the running point by a known scalar: the four stage scalars \(\lambda, \rho, 7^{-m}, 7^{mM_2}\) are fixed integers computed once inside preprocessing, and the two closing steps use scalars below \(2^{255}\) computed at run time. Every step therefore multiplies by a scalar of bit length at most \(255\), recoded as a width-\(5\) signed window (wNAF) digit string. A width-\(5\) wNAF of a \(255\)-bit scalar has at most \(44\) non-zero digits; the ledger charges \(64\). The recoding itself is \(255\) carry-propagation steps on four limbs, charged inside the per-step loop control.

Per scalar step the chain is:

- odd-multiple table of the running base: one doubling and \(15\) additions produce \(2Q, 3Q, \ldots, 31Q\), at most \(12 + 15\cdot 16 = 252\) field multiplications;
- main loop: at most \(255\) doublings and \(64\) additions, at most \(255\cdot 12 + 64\cdot 16 = 4084\) field multiplications;
- affine normalization, charged below.

So one scalar step costs at most \(252 + 4084 = 4336\) field multiplications before normalization, plus at most \((255 + 80)\cdot 32 + 64 = 10784\) field additions.

**Blocked Montgomery normalization.** A point must be affine before it is stored or compared. Instead of one Fermat inversion per point, each list is produced in blocks of \(2^{15}\) points. For each block the Montgomery trick forms the running products of the \(Z\) coordinates (\(1\) multiplication per point), inverts the last product once (\(608\) multiplications, charged as one Fermat inversion), walks the products back (\(2\) multiplications per point), and maps each point to affine coordinates (\(8\) further multiplications per point, matching the parent ledger). That is at most \(12\) field multiplications per normalized point, plus one inversion per block. The number of blocks is below \(2^{104}\), so all block inversions together cost below \(608 \cdot 2^{104} \cdot 383 < 2^{123.3}\) instructions. Because every \(Z\) is non-zero, the batch never divides by zero.

Charging \(12\) multiplications for normalization, one scalar step is at most

\[
4336 + 12 = 4348
\]

base-field multiplications, i.e. at most \(4348 \cdot 383 = 1\,665\,284\) multiplication instructions, plus at most \(10784 \cdot 64\) field-addition instructions.

Scalar-field elements (the units \(\lambda, \rho\), the root certification, and the final exponent \(k_0 + m k_1\)) are four-limb integers modulo \(r\). One schoolbook product is \(16\) limb products. The ledger charges \(1024\) instructions for a multiply-and-reduce, which covers `MUL`/`UMULH`, the carry chain, and a conditional subtract of \(r\). An exponentiation of a \(255\)-bit exponent is at most \(508\) such multiplications. The whole scalar-field setup, including the final reconstruction of \(s\), is charged as \(40\) exponentiations.

## Work at \(n = 2^{20}\)

Number of \(4\)-limb scalar steps:

\[
H_{\mathrm{steps}} = 2M + 2 \cdot 2^{10} + 2 = 2M + 2050.
\]

The two stages contribute \(M + M\) and \(2^{10} + 2^{10}\). The extra two steps are \(T \mapsto [7^{-k_0}]T\) and the final check \([s]P_0\).

Instruction total, with \(C = 2000\):

- curve chain: \(H_{\mathrm{steps}} \cdot 1\,665\,284\)
- field additions inside those operations: \(H_{\mathrm{steps}} \cdot 10784 \cdot 64\)
- block inversions for normalization: \(2^{124}\)
- mergesort of \(M\) records, \(118 = \lceil\log_2 M\rceil\) passes, \(24\) instructions to move one \(112\)-byte record (address update plus seven \(16\)-byte load/store pairs): \(118 \cdot M \cdot 24\)
- mergesort comparisons, \(64\) instructions each (full \(96\)-byte key): \(118 \cdot M \cdot 64\)
- binary search of each of the \(M\) giants through \(118\) affine keys: \(M \cdot 118 \cdot 64\)
- loop control on the two-limb index and the wNAF recoding, \(2^{12}\) instructions per scalar step: \(2^{12} \cdot H_{\mathrm{steps}}\)
- scalar-field exponentiations: \(40 \cdot 508 \cdot 1024\)
- a \(2^{24}\)-instruction lump for reading \(A_1\) and \(A_n\), reducing the final scalar, and writing the \(32\)-byte output

Let \(R\) be that sum. Every term is an integer upper bound on work the algorithm performs, including the failed comparisons inside sorting and search and the wNAF recoding of every step scalar. There is no second trial.

Numerically, the per-step instruction bound is

\[
B = 1\,665\,284 + 10784 \cdot 64 + 2^{12} = 2\,359\,556 < 2^{22},
\]

and \(H_{\mathrm{steps}} = 2M + 2050 < 2^{119}\), so the curve chain and loop control contribute

\[
H_{\mathrm{steps}} \cdot B = (2M + 2050) \cdot 2\,359\,556 < 1.056 \cdot 10^{42}.
\]

The remaining terms are below \(7.6 \cdot 10^{39}\) in total: the sort moves \(118 \cdot M \cdot 24\), its comparisons \(118 \cdot M \cdot 64\), the giant lookups \(118 \cdot M \cdot 64\), the block inversions \(2^{124}\), the scalar field \(40 \cdot 508 \cdot 1024\), and the lump \(2^{24}\). Hence

\[
R < 1.056 \cdot 10^{42} + 7.6 \cdot 10^{39} < 1.064 \cdot 10^{42} < 2000 \cdot 2^{129} = 1.3611 \cdot 10^{42},
\]

which is the required integer comparison \(R < 2000 \cdot 2^{129}\). Therefore

\[
\log_2 W(2^{20}) \le 129.
\]

Preprocessing (root certification, the two baby tables, the stage scalars and their wNAF recodings, and the scalar-field units) is part of \(R\), so it is at most \(2^{129}\) reference units and is not added again. Nonuniform advice is empty: the tables are computed inside \(R\), not embedded. \(\log_2(\max(1, 0)) = 0\).

Peak storage is the affine baby table, one mergesort buffer of the same size, one \(2^{15}\)-point inversion block of \(2^{15} \cdot 144\) bytes, and a few kilobytes of Jacobian scratch and code:

\[
2 \cdot 112 \cdot M + 2^{15} \cdot 144 + 2^{22} < 2^{127}
\]

bytes. Processors are one chain of dependent scalar steps; work is not divided away. Generation of the powers \(A_i\) is organizer work and is not charged. The attack's reads of \(A_1\) and \(A_n\) are inside the \(2^{24}\) lump.

No heuristic is used. Coverage is the mixed-radix division of a \(255\)-bit exponent, not a random subset of a coset and not an orbit under a proper subgroup of the \(d\)-th roots.
