# Floor-analysis appendix: exact sign-quotient operating point

This appendix reproduces the arithmetic needed for the package without an
external submission dependency. It is an implementation-family calculation,
not an information-theoretic or generic-group lower bound. A complete scalar
action price is an upper ceiling, not a theorem that every action costs it.

## Protected constants and exact rectangle

    r=52435875175126190479447740508185965837690552500527637822603658699938581184513
    d=1048576
    q=(r-1)/d=50006747412801924209068050869165387952509453297164571592906626415194112
    Q=q/2=25003373706400962104534025434582693976254726648582285796453313207597056

The exact factorization of q is

    2^12 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2
    * 2508409 * 2529403 * 52437899 * 254760293^2.

verify.py proves the listed factors prime by trial division and proves r
prime and g=7 primitive by the complete-factor Lucas criterion. The presence
of a residual 2^12 factor in q is essential: zeta=g^d has even order q,
and zeta^Q=-1. Smoothness alone does NOT implement multiplicative PH on an
encoded scalar; the mechanism uses only this order-two symmetry.

For target success p=39/100, quotient coverage must satisfy

    K=ceil(39Q/100)=9751315745496375220768269919487250650739343392947091460616792150962852,
    N=min{t:floor(t^2/4)>=K}=ceil(sqrt(4K))=197497501204408915131766844717197567,
    a=floor(N/2)=98748750602204457565883422358598783,
    b=ceil(N/2)=98748750602204457565883422358598784.

For fixed sum N, ab<=floor(N^2/4), attained by these balanced sides.
The script proves the chosen N meets K and N-1 fails. Since ab<Q,
quotient coverage is nonwrapping. Integer division k=ai+j partitions
[0,ab); the two exact lifts give 2ab covered residues in [0,q).
This is not a prime-q assumption, random mapping or asymptotic estimate.

    covered=2ab=19502631490992750441536539838974501476077154194472833185254723055359744
    exact success=1336529022134919849337756293789370989314497957406307098770197577807/3426997492653640639327580240485566608587544770913142241838447533936 >=39/100.

## Unchanged complete actions plus paid sign handling

    E=79650816
    routing=320*129*512=21135360
    metadata=65536
    additional sign record maintenance=256
    all-in per generated record=100851968
    fixed scored rows=21E+2^80+2^40
    R=N*100851968+21E+2^80+2^40=19918011671547009368992591426495508435762944
    C=2000, H=0, W=R/C=1244875729471688085562036964155969277235184/125.

The +256 is exactly retained, even if an optimized implementation might
reuse local metadata slack. It prices canonical y selection, orientation
encoding and copy control; it is not a second group operation. All actual
chain points stay original for subsequent complete scalar actions.

## Integer-powered claim and hundredth gates

The declared rational is 3331/25=133.24, gcd(3331,25)=1. The verifier uses
no floats. For any positive rational x=A/B and rational threshold h=u/v,
x<2^h is equivalent to A^v < B^v*2^u. Accordingly it checks

    R^25 < 2000^25 * 2^3331,
    (6R)^50 < (5*2000)^50 * 2^6657.

The second comparison is 6/5 repricing of ALL instructions below 2^133.14,
so the declared 133.24 retains at least 0.10 bits after that stress. The
immediately preceding hundredth would require stress below 2^133.13;
instead the exact opposite comparison holds:

    (6R)^100 >= (5*2000)^100 * 2^13313.

Thus 133.23 fails under the same headroom contract; changing the gate or
pricing only selected rows cannot manufacture a lower admissible claim.

For display, range-reduce x=2^k*y with 1<=y<2, u=(y-1)/(y+1). Use
ln(y)=2*sum(j>=0,u^(2j+1)/(2j+1)). After 36 positive terms, the remaining
sum lies between 0 and 2*u^73/(73*(1-u^2)). Divide the ln(y) enclosure by
the positive ln(2) enclosure and add k. The exact Fraction interval width
is <10^-25; all 12-place displays round OUTWARD. The enclosures are

| Quantity | Lower endpoint | Upper endpoint |
|---|---:|---:|
| raw log2(W) | 132.871197432118 | 132.871197432119 |
| log2((6/5)W) | 133.134231837952 | 133.134231837953 |
| stress plus 0.10 policy floor | 133.234231837952 | 133.234231837953 |
| claim minus policy floor | 0.005768162047 | 0.005768162048 |
| claim minus stressed work log | 0.105768162047 | 0.105768162048 |

The policy slack is in [0.005,0.06]. These values are independently
recomputed by this package's verify.py, not accepted on an imported JSON
claim. The construction has deterministic success and therefore no empirical
transfer-loss allowance or random-function rate to estimate.

## Real-EC extension of the algebra check

Every configuration uses y^2=x^3+4 over F_221201, subgroup order 12289,
P=(169338,64967), primitive scalar root 11. All chain actions and root
candidate actions use actual point multiplication; only organizer input
generation reuses a point table constructed by actual EC addition. Neither
attack routine receives that table, s, k1, k2 or an exponent oracle.

| d | q | a x b | exact hits / q | random hits / 2000 |
|---:|---:|---:|---:|---:|
| 16 | 768 | 12 x 13 | 312 / 768 | 821 / 2000 |
| 32 | 384 | 9 x 9 | 162 / 384 | 872 / 2000 |
| 64 | 192 | 6 x 7 | 84 / 192 | 892 / 2000 |
| 128 | 96 | 4 x 5 | 40 / 96 | 865 / 2000 |
| 256 | 48 | 3 x 4 | 24 / 48 | 1030 / 2000 |
| 512 | 24 | 2 x 3 | 12 / 24 | 1025 / 2000 |
| 1024 | 12 | 2 x 2 | 8 / 12 | 1345 / 2000 |
| 64 | 192 | 1 x 95 | 190 / 192 | 1982 / 2000 |

All configurations enumerate every k2 in addition to their 2000 random
full-pipeline trials. The total is 1716 exhaustive grid tests and 16000
random full-pipeline tests. Each configuration exactly matches its predicted
covered bins, has zero duplicate same-side keys, zero wrong exponents and
zero coverage-prediction errors. Random hits total 8832, misses 7168; these
unweighted aggregate frequencies mix different exact coverage fractions and
are NOT a production success estimate. A further exhaustive d=64 pipeline
enumerates all 12288 nonzero secrets: 5376 hits and 6912 intended misses,
again zero wrong scalar or exponent. The 6x7 grid alone gives 84/192 hits.

Captured JSON preserves per-k2 count/hit histograms, dense-bit and orientation
pair counts, action-count histograms, parameters/seeds and transcript digests.
The two configuration files are a lossless split of the same eight-row list.
The evidence supports only these real tiny-curve mechanisms. It does not
measure production instructions, memory hardware, runtime or production
parameter execution. Exact algebra, not sampling, establishes the claim.
