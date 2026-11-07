# Explicit cold field-basis conversion and its budget

The wrapper's map (no precomputed isomorphism is free; needs h and g irreducible, [0]/[0b]). In F[Y], F the protected
representation of F_(p^12), factor h(Y) = Y^6+Y^3+Y^2-Y-1 and take a root y (h splits over F), then a root x of
X^2+|u|y (it splits in F): Y -> y, X -> x is a field isomorphism, the 12 vectors y^i x^j an Fp basis; invert the 12 x 12
coordinate matrix by Gaussian elimination. One map serves both inputs; all is charged once per attack.

**Bounded factorization.** A completely splitting squarefree f of degree d <= 6 over F splits by gcd(R, f) or
gcd(R^((p^12-1)/2) - 1, f), R of degree < d uniform, repeated on unresolved factors: the values at the roots are jointly
uniform (Vandermonde), and two roots get opposite quadratic characters with probability ((|F|-1)/|F|)^2/2 > 1/4,
independently of earlier trials. Five splits serve h, one the quadratic; capped at 128 trials each, failure <= 6
(3/4)^128 < 2^-50 (counted). Uniform F elements use twelve Fp coordinates by 381-bit rejection sampling (<= 128
attempts, failure < 2^-128 each; < 2^16 coordinates in all, failure < 2^-112); a 2^-48 allowance covers both, and
random bits, rejections and failed trials are charged.

**Cost.** Modulo a degree-<= 6 polynomial a product takes < 2^10 field operations, powering with < 4572 bits < 2^24,
gcds and inversions < 2^24: 2^26 Fp12 operations per trial, < 2^36 over <= 6 x 128 < 2^10 trials, at < 2^31
instructions each below 2^67 instructions. Inversion, validation, pairing, scalar output, coordinate generation and
control fit with large slack in the 2^80-instruction wrapper cap. Maps and polynomials take < 2^30 bytes.

## Idempotent projection and protected anchor

In the common STNFS tower letN=p^12-1,h=N/r andk=h^(-1)modr. Exact arithmetic giveshmodr=13860579304294609274647100502445671956830243711357126442998 andk=24831094541369581385436675797109464480232629635954898292486708404843532461745. Usepi(z)=z^(hk). Its image lies inGT, it fixesGT, pi^2=pi, andh is a multiple ofp^6-1, so it kills nonzero Fp6 multipliers. The checker and independent Wolfram code/raw verify these exact identities.

For alpha0=1+x onf, Norm(alpha0)=1399 prime, and its simple unramified ideal isI=(1399,y+40,x+1). Its sigma conjugate hasx-1 and is distinct. Setomega=pi(eval_f(alpha0)). Exact tower arithmetic verifiesomega!=1,omega^r=1,sigma(omega)omega=1; its x-coefficient0 is1823330539042453828028069116870546032896926197478728071363869126167756345638209516204631600716902807865027915471103. It therefore generatesGT. Constructing and checkingomega, including norms and primes, uses fewer than2^14 tower products at the conservative2^31instructions each plus control, below2^50instructions inside the wrapper2^80 cap. Its ONE Schirokauer map is charged in map103.

Define the theoretical true logL_K(alpha)=log_omega(pi(eval_K(alpha))). It is already sigma-odd; if an odd projection is written, use(L-Lsigma)/2. In this normalizationL_f(alpha0)=1. Orient the protected orbit byI itself, overriding any default sigma representative, and protectI and all21 Schirokauer coordinates through filtering. The class/unit prerequisites in proof.md give the exact identityv_I+ell_f.SM_f(alpha0)=1. For any extracted nonzero kernel vectorw=c*v, eta=w_I+w_SMf.SM_f(alpha0)=c!=0; dividew byeta. This computes the true normalized advice with charged field operations, without a class-number, unit-regulator or DLOG oracle. The proof and live checker use the same sign/half conventions. For the actual bootT^e with uniformnonzeroe, divide its computed log bye; the Fp6 lift multiplier addszero. Additive blinding byomega^b, if used, would addexactlyb, but is not the chosen boot. Both logarithms share this normalization; their ratio recovers the supplied-Q exponent.
