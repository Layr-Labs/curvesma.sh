# Exact Fp4 trace key for full twelve-element GT orbits

This is an algebraic representation change, not a projection heuristic. It uses the full protected p,r,raw GT states and one protected reference multiplication per transition. The raw one-sided independent-increment walk and its ideal collision probability are unchanged. Hash values and actual trajectories change with the exact key and NH fingerprint; their ideal coupling is proved explicitly. No auxiliary powers, searched advice or unchecked relation are introduced.

## Trace polynomial and exact three-element fibers

Let sigma(X)=X^p, xi=p modr of order12. For X in GT define

 c(X)=X+sigma^4(X)+sigma^8(X) in Fp4.

The norm over Fp4 is X^(1+p^4+p^8)=1, since 1+xi^4+xi^8=0. Also sigma^6(X)=X^-1. Therefore the three conjugates are precisely the roots of

 F_X(T)=T³-c(X)T²+c(X)^(p²)T-1.

Indeed the quadratic coefficient is their sum, the constant their norm, and the pair-product sum is their inverse sum, equal to sigma²(c(X)) because {xi²,xi⁶,xi¹⁰}={-xi⁸,-1,-xi⁴}. No probability statement is used.

For nonidentity X the three conjugates are distinct: a sigma^4 fixed point in the prime order r subgroup would imply r divides p^4-1, contrary to embedding degree12. If Y in GT has the same c, it is a root of the SAME cubic and hence is one of X,sigma^4X,sigma^8X. Conversely those three have the same trace. Thus this trace has exactly the desired three-element fibers on nonidentity GT.

This is the characteristic-polynomial mechanism of XTR. Material external source: A.K. Lenstra and E.R. Verheul, The XTR public key system, CRYPTO2000, Lemma2.2.1, https://www.cs.ru.nl/E.Verheul/papers/crypto2000/crypto2000.pdf . Their usual notation uses a degree-six extension; here the field-size parameter is p². The identities above, native coefficient projection below, full twelve-orbit key and complete work bound are derived explicitly for this protected profile. No XTR performance theorem, novelty or scientific qualification is claimed from the citation.

## Native coordinate extraction, not a free trace primitive

Use the pinned native tower Fp2=Fp[u]/(u²+1), Fp6=Fp2[v]/(v³-1-u), Fp12=Fp6[w]/(w²-v). Write

 X=a0+a1*w²+a2*w⁴+b0*w+b1*w³+b2*w⁵.

Under sigma^4 the coefficient field Fp2 is fixed and w is multiplied by a primitive cube root. Thus the sum over powers0,4,8 vanishes on degrees1,2,4,5 and is3 on degrees0,3. Consequently

 c(X)/3 = a0+b1*w³ = S(X).

Here3 is invertible and need not be multiplied or divided in the implementation. We use S itself as the key precursor. In the actual576-byte native state, a0 is at0..95 and b1 at384..479. The current implementation borrows a0 directly and copies b1 into candidate zero with six LDP/six STP, fusing its identity OR reduction. Full pointer/copy/test work is paid in canonical-parent.md. There is no call to an unpriced trace, exponentiation or full Frobenius routine. The pinned tower definition is at https://raw.githubusercontent.com/supranational/blst/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fp12_tower.c .

For nonidentity GT, b1 cannot be zero. Otherwise c lies in Fp2, so c^(p²)=c and F_X factors as

 (T-1)(T²+(1-c)T+1).

A nonidentity root then has degree at most2 over Fp2, so belongs to Fp4; but GT intersects Fp4 only at1. Contradiction. Thus b1=0 identifies the identity on GT. Its visit is the same charged rare identity event as before, not an additional exception.

## Four exact Frobenius images identify the full twelve-orbit

Set t=w³. Then t²=1+u, Fp4=Fp2[t], and sigma²(t)=-t. With A=a0,B=b1, S=A+Bt has four images

 (A,B), (conj(A),K(B)), (A,-B), (conj(A),-K(B)),

where K(B)=C*conj(B), C=(1+u)^((p-1)/2). Since p=19 mod24, C=(s,s), with s nonzero in Fp and2s²=-1. Therefore K(b_re,b_im)=(s(b_re+b_im),s(b_re-b_im)), two native Fp products and two add/sub operations. These constants are constructed cold.

B!=0 implies all four S-images are distinct: sigma² changes B to-B, and a fixed point of sigma or sigma³ would also be fixed by sigma². Choose the unique minimum in fixed unsigned native Montgomery coefficient order (A_real,A_imag,B_real,B_imag), and let k in0..3 be its applied exponent. Then k(sigma X)=k(X)-1 MOD4. The minimum S is invariant on the full twelve-element orbit of X. Conversely equal canonical S implies equality of traces after one of four twists; the cubic fiber theorem then implies equality of raw points after an additional twist4m. Thus canonical S keys are EQUAL IFF raw GT points lie in the same full Frobenius12 orbit.

No claim that k transforms modulo12 is made or needed. Before the first orbit repeat, the fresh jump index is uniform independently of k and the group past. k is used only in terminal twist recovery. The raw independent-increment Markov law follows from fresh labels, not orientation equivariance. The algorithm stores EVERY visited orbit and stops at its first repeat, so no quotient-cycle or post-repeat equivariance assumption enters.

## Exact terminal twist recovery and G2 output

At a matching canonical key, saved/current k give the three possible raw twists

 ell=k_saved-k_current+4m mod12, m=0,1,2.

Reconstruct the saved raw point X_saved=h^A_saved*z^B_saved using its already stored coefficients, with at most2048 paid dense products and fixed interfaces. Test these three twists against current full raw X; exactly one matches because nonidentity GT has free twelve-orbits. The needed full Frobenius constants are initialized once cold,without constructing rotated jump records. This is terminal fixed cold work, not a missing per-step snapshot or uncharged decompression oracle. It is far below the existing complete cold allowance; finishing, pairings, validation and verification still fit the prior2^120 term.

For this true group-determined ell, form D=B_current-xi^ell*B_saved. On D0 fail; otherwise a=(xi^ell*A_saved-A_current)/D, verify h^a=z, form canonical Q=[a]Q0, check curve/subgroup and pairing(Q,P)=z, then output. The reconstruction depends on coefficient values only through their invariant group point, so the selected true twist is group-path determined and the finite Frobenius interval/kernel singular-match bound in probability.md applies. Scalar-only or unchecked output is never accepted.


Complete work and all statuses are specified once in proof.md and scalar-audit.md;success.md gives the two-term finite certificate.
