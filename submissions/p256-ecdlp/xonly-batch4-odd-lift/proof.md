# Odd-lift BSGS with x-only conversion

NIST P-256, short Weierstrass \(y^2 = x^3 + ax + b\) over \(F_p\), \(a = -3\), cofactor 1, subgroup order \(r\) prime. The public generator is the SEC 2 point \(P\). Input \(Q = [x]P\) with \(x\) uniform in \(\{0,\ldots,r-1\}\), including the identity. Output the unique such \(x\). No auxiliary input.

Reference unit: \(C = 1000\) scalar AArch64 instruction equivalents borrowed from the secp256k1 normalization. P-256 work is charged as justified instruction counts divided by \(C\). A point addition is not one unit. No reference-operation count \(H\) is added on top of those instructions.

## Scalar partition

\(r\) is odd. If \(Q\) is infinity, return \(0\). Otherwise let \(z = x\) when \(x\) is odd and \(z = x + r\) when \(x\) is even. Then \(z\) is odd, \(0 < z < 2r < 2^{257}\), and \([z]P = Q\).

Fix \(m = 128\). Every odd integer in \((0, 2^{2m+1})\) is uniquely
\[
z = (2i+1)\,2^{m} + \sigma(2j+1),
\]
with \(i \in [0, 2^{m})\), \(j \in [0, 2^{m-1})\), \(\sigma \in \{+1,-1\}\).

The intervals \((2i \cdot 2^{m},\ (2i+2)\,2^{m})\) partition the positive integers below \(2^{2m+1}\). Inside the interval for a fixed \(i\), the minus branch is the odd integers in \((2i\cdot 2^{m},\ (2i+1)2^{m})\) and the plus branch is the odd integers in \(((2i+1)2^{m},\ (2i+2)2^{m})\). The shared endpoint \((2i+2)2^{m}\) is even. The map is therefore a bijection onto the odds below \(2^{257}\). Since \(z < 2r < 2^{257}\), exactly one triple \((i,j,\sigma)\) represents \(z\).

Babies are \(B_j = [2j+1]P\) for \(j < 2^{127}\). Their scalars lie in \((0, 2^{128})\). They are distinct modulo \(r\) because the scalars differ by less than \(r\). No baby is the negation of another: that would require two odd scalars below \(2^{128}\) to sum to \(r > 2^{255}\). Distinct points that are not negations have distinct affine \(x\)-coordinates, and negation preserves \(x\), so the \(2^{127}\) baby \(x\)-coordinates are distinct.

Giants are \(G_i = [(2i+1)2^{128}]P\). Then \(Q - G_i = \sigma B_j\) for the true triple, so the affine \(x\)-coordinate of \(Q - G_i\) equals the affine \(x\)-coordinate of \(B_j\). A hit determines \(j\). The two candidate scalars \((2i+1)2^{128} \pm (2j+1)\), reduced modulo \(r\), are tested by two scalar multiplications. Exactly one reproduces \(Q\): their difference is twice a baby, and \(r\) is an odd prime, so no baby has order 2.

## Schedules

Field elements are four 64-bit limbs in Montgomery form, representative in \([0,p)\). Public constants \(R^2 \bmod p\), \(p\), and \(n_0' = -p^{-1} \bmod 2^{64}\) are profile constants, not advice. The generator coordinates are entered once by multiplication by \(R^2\).

Jacobian coordinates stay in Montgomery form. Addition and doubling use field addition, subtraction, multiplication and squaring only. Multiplication by 2, 4 or 8 is repeated field addition. The \(a = -3\) doubling identity \(3X^2 + aZ^4 = 3(X^2 - Z^4)\) holds on Montgomery representatives because squaring and addition commute with the factor \(R\).

Baby chain. Compute the affine point \(2P\) by one doubling. Set \(A_0 = P\) with \(Z = 1\). For \(j = 0,\ldots,2^{127}-1\), store \(A_j\) and set \(A_{j+1} = A_j + 2P\) by mixed Jacobian-affine addition. The last sum is discarded. That is \(2^{127}\) mixed additions. Convert the stored points to canonical affine \(x\) in one Montgomery batch, then one Montgomery multiplication by 1 per \(x\) to leave Montgomery form.

Giant chain, in blocks of \(L = 2^{12}\). Build \(S = [2^{128}]P\) by 128 doublings from \(P\), then the stride \(2S\) by one doubling, and negate both \(y\)-coordinates (one field subtraction each). Start from \(Q + (-S)\). Each subsequent point adds the fixed affine stride \(-2S\). There are \(2^{128}\) stored giants and \(2^{128}\) stride additions, plus the initial mixed addition,  \(2^{128}+1\) mixed additions in all. Every block of \(L\) Jacobian points is converted x-only, canonicalized, looked up, and discarded. Blocks do not change the operation count; they bound live memory and replace one inversion of \(2^{128}\) field elements by \(2^{116}\) inversions.

Lookup. Baby records are pairs \((x, j)\): four limbs plus two limbs, 48 bytes. They are sorted by canonical \(x\) with the top-down mergesort that splits at \(\lfloor n/2\rfloor\). Each giant \(x\) is looked up by the half-open binary search below. A hit is verified by the two candidate scalars. Each verification scalar multiplication is at most 256 Jacobian doublings and 255 mixed additions. Two verifications are 512 doublings and 510 mixed additions. Together with the 130 setup doublings this is at most 1024 doublings and is absorbed in the additive term of the instruction total, not in the \(2^{128}\) coefficient.

## Field-operation shapes

Mixed Jacobian-affine addition, finite non-doubling branch. Four squarings \(Z_1^2\), \(H^2\), \(r^2\), \((Z_1+H)^2\). Seven multiplications \(X_2 Z_1^2\), \(Z_1^3\), \(Y_2 Z_1^3\), \(H\cdot 4H^2\), \(X_1\cdot 4H^2\), \(r(V-X_3)\), \(Y_1 J\). At most 16 field additions or subtractions implement \(H\), the two doublings that build \(4H^2\), \(S_2-Y_1\), the doubling to \(r\), \(2V\), the two subtractions that build \(X_3\), \(V-X_3\), the doubling of \(Y_1 J\), \(Y_3\), \(Z_1+H\), and the two subtractions that build \(Z_3\). If the branch discovers equal summands it also performs one doubling; the charge for every mixed addition is therefore one full mixed addition plus one full doubling. The infinity branch returns earlier and is cheaper.

Doubling, \(a = -3\). Seven squarings and two multiplications: \(Z^2\), \(Y^2\), \(X^2\), \(Z^4\), \(X^2 Y^2\), \(\alpha^2\), \(Y^4\), and \(\alpha(S-X_3)\). At most 24 field additions build \(3(X^2-Z^4)\), \(4X Y^2\), \(2S\), \(X_3\), \(S-X_3\), \(8Y^4\), \(Y_3\), and \(Z_3 = (Y+Z)^2 - Y^2 - Z^2\).

x-only batch of \(k\) finite points. Montgomery inversion of the \(k\) values of \(Z\): \(k-1\) prefix multiplications, one multiplication by the last \(Z\), one inversion, \(k-1\) unwind multiplications, then \(k\) squarings and \(k\) multiplications for \(x = X Z^{-2}\). That is \(3k-1\) multiplications, \(k\) squarings and one inversion. Canonicalization is \(k\) further multiplications and is counted separately. Babies use \(k = 2^{127}\). Giants use \(2^{116}\) blocks of \(k = 2^{12}\), which is \(3\cdot 2^{128} - 2^{116}\) multiplications and \(2^{128}\) squarings.

Inversion of one field element is Fermat exponentiation by \(p-2\) in Montgomery arithmetic. 
\[
p-2 = \texttt{0xffffffff00000001 00000000 00000000 00000000 ffffffff ffffffff fffffffd}.
\]
The eight 32-bit groups have Hamming weights 32, 1, 0, 0, 0, 32, 32, 31, total 128, and the bit length is 256. Binary exponentiation uses 255 squarings and 127 multiplications, hence 382 field multiplications at the squaring price below. One baby inversion and \(2^{116}\) giant-block inversions.

## Instruction prices

A field multiplication or squaring is CIOS on four 64-bit limbs. One limb step is the straight-line sequence MUL, UMULH, ADDS, ADCS, ADC, ADDS, ADCS, ADC, ADDS, ADC, MOV, MOV, MOV: 13 instructions. It adds a 64-bit limb and a carry of at most 65 bits into a 128-bit product. The sum is strictly less than \(2^{129}\), so the captured extra bit is 0 or 1. Each of the four outer iterations runs four multiply-accumulate limb steps and four modulus limb steps. Charging 32 instructions for every limb step, rather than 13, costs \(4 \cdot 4 \cdot 32 + 4 \cdot 4 \cdot 32 = 1024\) instructions. The word \(m = t_0 n_0' \bmod 2^{64}\), the final carry writes, a fold of a top limb in \(\{0,1\}\) by the fixed residue of \(2^{256}\), and two conditional subtractions of \(p\) add at most 96 instructions. The functional check rejects a top limb above 1 or a value still at least \(p\) after those two subtractions. Operand traffic that the arithmetic schedule does not already include is at most eight loads, four stores and a handful of address adds. The charge per field multiplication, and per squaring, is **1280** instructions.

A field addition or subtraction is four ADDS/ADC or SUBS/SBC, a conditional subtraction of \(p\), and the loads and stores of three field elements. The charge is **32** instructions.

A canonical \(x\) comparison, including the address of record \(mid\) (a shift-add by 48), four loads and four compare/branch instructions, is charged **48** instructions.

A 48-byte record move is six loads, six stores and pointer arithmetic, charged **32** instructions.

Loop overhead, block indexing and the final verification doublings are not given their own \(2^{128}\) coefficient beyond the explicit 128-instruction allowance below. The additive \(2^{25}\) instructions in the total covers 1024 doublings: a doubling is at most \(9\cdot 1280 + 24\cdot 32 = 12288\) instructions, and \(1024 \cdot 12288 < 2^{24}\).

## Operation counts and instruction total

Let \(M = 1280\) and \(A = 32\). One charged mixed addition costs
\[
(11M + 16A) + (9M + 24A) = 14592 + 12288 = 26880
\]
instructions. Mixed additions charged: \(2^{127} + 2^{128} + 511 = \tfrac32 2^{128} + 511\).

Instruction total \(T\), with coefficients of \(2^{128}\):

| Term | Coefficient of \(2^{128}\) |
| --- | ---: |
| Mixed additions, \(\tfrac32 \cdot 26880\) | 40320 |
| Batch multiplications and squarings, \(6M\) (upper bound on \(6\cdot 2^{128} - 2^{116} - 1\)) | 7680 |
| Canonicalization, \(\tfrac32 M\) | 1920 |
| Fermat inversions, \(382 M \cdot 2^{116} / 2^{128} < 120\) | 120 |
| Binary search, \(128 \cdot 48\) probes per giant | 6144 |
| Mergesort comparisons, \(126 \cdot 48 / 2\) | 3024 |
| Mergesort moves, \(127 \cdot 32 / 2\) | 2032 |
| Initial record stores, \(32 / 2\) | 16 |
| Giant block spill, baby \(X,Z\) store, loop allowance | 224 |
| **Sum** | **61480** |

The mergesort comparison count is \(n\lceil\log_2 n\rceil - 2^{\lceil\log_2 n\rceil} + 1\) at \(n = 2^{127}\), which equals \(126\cdot 2^{127} + 1\). Dividing by \(2^{128}\) leaves the coefficient \(126\cdot 48/2\) plus one comparison in the additive term. Binary search on \(2^{127}\) keys performs at most \(\lfloor\log_2(2^{127})\rfloor + 1 = 128\) comparisons; the half-open loop shrinks a power-of-two interval and the closed form was checked against the recurrence for every width the experiment runs.

Thus
\[
T < 61480 \cdot 2^{128} + 2^{25}.
\]
Since \(64000 - 61480 = 2520\) and \(2^{25} < 2520 \cdot 2^{128}\),
\[
T < 64000 \cdot 2^{128} = 1000 \cdot 2^{134}.
\]
Work \(W = T/1000 < 2^{134}\). The submitted bound is \(2^{134}\) target-group-operation equivalents. No second copy of this work is charged as reference additions.

## Preprocessing, memory, advice, success

Preprocessing is the baby chain, the baby batch and canonicalization, the one baby inversion, and the mergesort. Its \(2^{128}\) coefficient is
\[
26880/2 + 2M + M/2 + 3024 + 2032 + 16 + 32 + 64 + 1 = 21809 < 32000.
\]
The same additive \(2^{25}\) argument yields preprocessing work \(< 2^{133}\), and that work is inside \(T\).

Memory schedule. While the baby chain runs, the only large buffer is the Jacobian triple array, \(96 \cdot 2^{127} = 3\cdot 2^{132}\) bytes. After the chain, \(Y\) is dropped; the batch prefix reuses that space. The mergesort then holds two record arrays, another \(96\cdot 2^{127}\) bytes, and does not keep the Jacobian array. Giant blocks hold \(2^{12}\) points, under \(2^{20}\) bytes, beside the finished table of \(48\cdot 2^{127}\) bytes. Code and registers are under \(2^{20}\) bytes. Peak storage is therefore at most \(3\cdot 2^{132} + 2^{20} < 4\cdot 2^{132} = 2^{134}\) bytes.

Advice is 0 bytes. The modulus, the generator, \(R^2\) and \(n_0'\) are profile constants. Nothing is searched and embedded.

The procedure is deterministic. Coins are not used. Every \(x \in [0,r)\) has a proof above, so the success probability on a uniform challenge is 1. There is no heuristic, no abort, and no restart.

## What the experiment checks

`experiment.py` exhausts a prime-order toy curve, checks the odd partition at widths 3, 4 and 5, checks that the mergesort comparison count matches the recurrence, checks four-limb CIOS against \(ab\cdot 2^{-256} \bmod p\) on 49 products, and checks the integer inequalities \(61480 < 64000\), preprocessing coefficient \(< 32000\), and \(3\cdot 2^{132}+2^{20} < 2^{134}\). It does not measure AArch64 and it does not run the production loop.
