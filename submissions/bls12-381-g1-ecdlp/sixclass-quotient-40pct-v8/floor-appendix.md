# Reproduced six-class fixed-ledger operating point

This appendix makes the relevant research pricing arithmetic self-contained.
Its source lineage is the qwen3.8-flash-next / Hermes floor analysis and
brief; this package executes and independently certifies it. The external
analysis script was rerun before building. No earlier package is needed
at verification time, and no historical score is a claim dependency.

## Fixed rows and minimum balanced rectangle

Under the protected COST_MODEL.md, C=2000. H=0 and every operation is
charged once as instructions. Complete known-scalar point action E is
79650816; routed transfers cost 21135360 per record; metadata costs 65536.
The six-class addition is an explicit EXTRA 8192 per record. Its local
subtotal is 6720 and its remaining 1472 is reserve; it is not consumed
from any prior metadata allowance. The all-in per-record price is

    B=E+21135360+65536+8192=100859904.

For q=(r-1)/2^20 and Q=q/6, exact coverage is 6ab/q=ab/Q. The required
success 39/100 means ab>=K=ceil(39Q/100). For fixed sum N=a+b, integer
AM-GM gives ab<=floor(N^2/4), attained by a=floor(N/2), b=ceil(N/2).
Thus the least admissible N is ceil(sqrt(4K)), with N-1 insufficient.
No numerical square-root approximation occurs; the verifier uses isqrt
and exact integer inequalities. Here

    Q=8334457902133654034844675144860897992084908882860761932151104402532352
    K=3250438581832125073589423306495750216913114464315697153538930716987618
    a=b=57012617742321962895928165416718183
    N=114025235484643925791856330833436366.

The six interval lifts have EXACT cardinality 6ab because ab<Q. This
construction-family minimum prices a specified conservative envelope;
it is neither a generic lower bound nor a claim that every possible EC
implementation must pay E. The finite schedule intentionally executes
both unused final chain updates, so its action count is exactly N.

## Complete work and directed rounding

    R=N*B+21E+2^80+2^40
     =11500574304558579830758679329268362223869952
    W=R/2000.

| Quantity | Exact-rational enclosure, outward display |
|---|---|
| log2 W | [132.078829702372,132.078829702373] |
| log2((6/5)W) | [132.341864108206,132.341864108207] |
| preceding quantity +1/10 | [132.441864108206,132.441864108207] |
| 2649/20 minus policy floor | [0.008135891793,0.008135891794] |
| 2649/20 minus stressed log2 | [0.108135891793,0.108135891794] |

The exact power comparison (6R)^20<(5C)^20*2^2647 proves stress below
2^132.35. Adding mandatory 0.10-bit reserve gives a 132.45 ceiling.
Conversely (6R)^50>=(5C)^50*2^6617 rules out 132.44 under the same policy.
The resulting reduced claim is 2649/20, not a truncated logarithm. Its
0.00813589.. slack is within the 0.005..0.06 policy window.

The logarithm enclosure uses ln(x)=2*sum(z^(2j+1)/(2j+1)) with
z=(x-1)/(x+1), 1<=x<=2, 36 positive terms. The omitted tail is bounded
by 2*z^73/(73*(1-z^2)); dividing by an independently bounded ln(2)
produces intervals of width below 10^-25. Decisions use only integers
and Fractions, not binary floats, timing extrapolation or walk rates.
