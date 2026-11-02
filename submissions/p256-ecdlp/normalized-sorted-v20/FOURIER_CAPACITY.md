# Fixed-bit capacity lemma

Let p be the protected P-256 field prime, r its prime curve order, R0=2^256, and rho=(R0 mod p)^(-1) mod p. For canonical Montgomery keyk, the number of affine curvepoints is1+chi((rho*k)^3-3rho*k+b). Define g(k) to be this quadratic character for0<=k<p and0 forp<=k<R0.

Primary reference for the mixed character bound: Todd Cochrane, Exponential sums modulo prime powers, Acta Arithmetica101.2(2002), introduction equation(1.2), author PDF https://www.math.k-state.edu/~cochrane/research/expprpowa.pdf . A nonzero additive frequency has polynomial degree1 and the multiplicative cubic degree3, characterorder2, giving3sqrt(p). At zero frequency Hasse gives2sqrt(p). The cubic is nonsingular and not a square.

An arbitrary realfrequency exponential over0..p-1 can be expanded in the p additive characters. The geometric coefficients satisfy min(1,1/(2p||theta-t/p||)). Counting at mosttwo frequencygridpoints per distanceannulus bounds their normalizedL1 sum by2+H_floor(p/2)<3+lnp; the deliberately looser4+2lnp is also safe. Hence every R0-DFT coefficient of zero-extendedg has magnitude at most3sqrt(p)(4+2lnp).

For an interval of lengthL in Z/q, its normalized FourierL1 norm is at most2+lnL: use|D_L(t/q)|<=min(L,q/(2|t|)) and the decreasing-sum integral. Fixing a bitblock[a,b) is a periodic interval of length2^a and period2^b; periodic lifting preserves normalizedL1. The lowest-bit fixedblock has norm exactly1. The complete bitmask is a POINTWISE PRODUCT of these indicators, so the Fourier convolution inequality bounds its normalizedL1 by the product of these norms. Do not use the false normalizedL1 product claim for pointwise Fourier products of AP convolutions.

For SIDE35 and d3=26, bucket bits are[0,35),[64,99),[128,154). There are2^160 total integerkeys in the cell. Its indicator has norm at most(2+64ln2)(2+128ln2). Pairing againstg gives pointcount <=2^160+3sqrt(p)(4+2lnp)(2+64ln2)(2+128ln2). Zero extension avoids wrap or proper-GAP assumptions. The count of canonical integers in the cell is at most2^160.

Safe rational constants: ln2<7/10, lnp<180, sqrt(p)<2^128,2+64ln2<47,2+128ln2<92. Error<3*364*47*92*2^128=4721808*2^128<2^151. Thus use KB96=(513/512)*2^160/(r-1). For either35-bit slot, norm<=47, error<3*364*47*2^128<2^144, far below2^221/512; the same513/512 coefficient is safe. Retain crude two-point fiber bounds for all other patterns, particularly144and182bits. This explicit lemma is used by proof.md; its argument and cited character-sum theorem are material evidence.

Conditional leaves with one free independent generator are uniform or uniform excluding one group point; denominatorr-1 is retained. Recheck every dependency/variance step and every positive capgap when applying the new KBs. Do not merelysquare a negative gap.
