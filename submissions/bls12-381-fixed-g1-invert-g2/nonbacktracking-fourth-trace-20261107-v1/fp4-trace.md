# Exact Fp4 trace key for full twelve-element GT orbits

This is an algebraic representation change, not a projection heuristic. It uses the full protected p,r,raw GT states and one protected reference multiplication per transition. The raw nonbacktracking walk and its ideal collision probability are unchanged. Hash values and actual trajectories change because the exact key and tabulation positions change. No auxiliary powers, searched advice or unchecked relation are introduced.

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

Here3 is invertible and need not be multiplied or divided in the implementation. We use S itself as the key precursor. In the actual576-byte native state, a0 is at0..95 and b1 at384..479. Extracting these two96-byte blocks is twelve LDP/twelve STP plus paid pointer setup, at most27 scalar instructions. There is no call to an unpriced trace, exponentiation or full Frobenius routine. The pinned tower definition is at https://raw.githubusercontent.com/supranational/blst/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fp12_tower.c .

For nonidentity GT, b1 cannot be zero. Otherwise c lies in Fp2, so c^(p²)=c and F_X factors as

 (T-1)(T²+(1-c)T+1).

A nonidentity root then has degree at most2 over Fp2, so belongs to Fp4; but GT intersects Fp4 only at1. Contradiction. Thus b1=0 identifies the identity on GT. Its visit is the same charged rare identity event as before, not an additional exception.

## Four exact Frobenius images identify the full twelve-orbit

Set t=w³. Then t²=1+u, Fp4=Fp2[t], and sigma²(t)=-t. With A=a0,B=b1, S=A+Bt has four images

 (A,B), (conj(A),K(B)), (A,-B), (conj(A),-K(B)),

where K(B)=C*conj(B), C=(1+u)^((p-1)/2). Since p=19 mod24, C=(s,s), with s nonzero in Fp and2s²=-1. Therefore K(b_re,b_im)=(s(b_re+b_im),s(b_re-b_im)), two native Fp products and two add/sub operations. These constants are constructed cold.

B!=0 implies all four S-images are distinct: sigma² changes B to-B, and a fixed point of sigma or sigma³ would also be fixed by sigma². Choose the unique minimum in fixed unsigned native Montgomery coefficient order (A_real,A_imag,B_real,B_imag), and let k in0..3 be its applied exponent. Then k(sigma X)=k(X)-1 MOD4. The minimum S is invariant on the full twelve-element orbit of X. Conversely equal canonical S implies equality of traces after one of four twists; the cubic fiber theorem then implies equality of raw points after an additional twist4m. Thus canonical S keys are EQUAL IFF raw GT points lie in the same full Frobenius12 orbit.

No claim that k transforms modulo12 is made or needed. Before the first orbit repeat, the fresh twelve-way phase label minus any deterministic group-only k is still uniform modulo12. The raw nonbacktracking Markov law is exactly the same. The algorithm stores EVERY visited orbit and stops at its first repeat, so no quotient-cycle or post-repeat equivariance assumption enters.

## Exact terminal twist recovery and G2 output

At a matching canonical key, saved/current k give the three possible raw twists

 ell=k_saved-k_current+4m mod12, m=0,1,2.

Reconstruct the saved raw point X_saved=h^A_saved*z^B_saved using its already stored coefficients, with at most2048 paid dense products and fixed interfaces. Test these three twists against current full raw X; exactly one matches because nonidentity GT has free twelve-orbits. The needed full Frobenius constants were constructed with the jump variants. This is terminal fixed cold work, not a missing per-step snapshot or uncharged decompression oracle. It is far below the existing complete cold allowance; finishing, pairings, validation and verification still fit the prior2^120 term.

For this true group-determined ell, form D=B_current-xi^ell*B_saved. On D0 fail; otherwise a=(xi^ell*A_saved-A_current)/D, verify h^a=z, form canonical Q=[a]Q0, check curve/subgroup and pairing(Q,P)=z, then output. The reconstruction depends on coefficient values only through their invariant group point, so the selected true twist is group-path determined and the previous uniform-denominator/fruitless proof applies unchanged. Scalar-only or unchecked output is never accepted.

## Complete key/hash/storage changes

The key is canonical S in192 initialized native bytes: four canonical Fp coefficients,24 little-endian64-bit words. Tabulation uses48 independent32-bit digit tables and one constant-zero digit table T48, still320-bit entries. Construct all(48*2^32+2)entries cold, including unused T48[1]. TABLE message/tag/part convention is unchanged, except positions are0..48. Distinct full12-orbit keys differ in a real digit; pairwise fingerprint collision is still2^-320 under the SAME standard BLAKE3 bridge H. No new trace-distribution, independence or compression heuristic is used.

The paired tabulator has24 PAIR groups, no key/base rebase instructions,4 constant-accumulator loads,5REV,3output stores andRET:24*22+4+5+3+1=541 scalar instructions. All digit offsets<=184 and table-base offsets<=368 fit their ordinary paired encodings. The U40 big-endian serialization and full standard JUMP hash are unchanged. The conditional phase/bucket proof remains valid: fresh labels depend only on the exact full-orbit key, k is group-only, and actual execution stops at first repeat.

Payloads remain512 bytes with existing A/B/k/U/next offsets. Store the192-byte canonical S at0..191;192..303 is initialized zero padding. A304,B336,k368 (now0..3),U384,next432 are unchanged. Compare24 canonical words at the once-only U confirmation; mismatched key aborts as before. Reconstruct saved raw X only at terminal match. Payloadwrite224 and head/chain bounds remain conservative with fewer copied key bytes. Shared tables are smaller and still below the previous shared2^45 capacity. Physical payload/head/jump capacity remains<2^136 bytes, zero advice.

## Deterministic canonical instruction bound

Always compute K(B),-B,-K(B) and-A_imag:2native products and7field operations (twoadd/sub,fivenegations). No cost counter, orientation average or canonical-budget failure is used. Native product384 and fieldleaf42 remain fully audited. One fixed cold-built field program has9calls; interpreter/entry/interfaces<=58, parent preservation<=16, top interface<=16, total90<128 wrappers. Actual helper ABIs are encoded per descriptor. All constants/programs and working buffers are initialized cold.

trace-min.S supplies the complete minimum leaf. If A_imag!=0, its comparison to its negative chooses parity, then only B versus-B in that parity is compared. If A_imag=0, all A prefixes coincide and all four B candidates are compared. Zero/nonzero internal paths count131/88 including all result copies andRET. Add16 interface and18 initial B-copy instructions, maximum165; retain176. All paired offsets are aligned/in range, borrow carry is used correctly, x18 untouched, caller registers only, no hidden frame or spill. Output192 bytes is fully initialized and disjoint from inputs.

Remaining control: native A/B extraction27, B identity test18, fixed program/buffer setup<=10, parent context<=16, final dispatch<=12,total83<128. Minimum/interface/Bcopy are separate in176; field-program/ABI work separate in128. Total canonical<=2*384+7*42+176+128+128=1494 per transition, plus its cold initial evaluation.

The old product counter is removed. trace-hot-loop.S makes the complete ordinary controller explicit:77 main instructions on its longest ordinary branch and17 in the called jump selector,94 total. Retain384 for the controller with290 spare instructions; all helper bodies are separately charged. Retry interface/rechecks are included in the separately capped extra-attempt work. The final exhausted check is fixed cold work. trace-controller.md gives ABI, context offsets, exact phase arithmetic and bank-selector details.

Full ordinary ledger1494canonical+541tabulation+797jointlabelhash+46scalar+384controller+192head+224payload+32cappedtraversals+1retry=3711. Retain3776K,65 slack. Complete coldD=540*2^119 includes all head zeros, table/generator/scalar/variant construction, attempts, new constants/code, initial state/key, terminal raw reconstruction/three-twist search, failures, pairings and exact canonical G2 output. The prior largest construction terms<2^108 and bounded terminal terms<2^82 still sum far below2^120, so the same complete D certificate holds.

With one protected reference per transition,

 W<=(23776K+540*2^119)/20000.

Exact protected-parameter integer certificate(23776K+540*2^119)^1000<20000^1000*2^125879 establishes proposed125879/1000=125.879,above diagnostic125.87866896466751. All previous table/mixing/collision/fruitless/phase/retry/dictionary failure losses are retained. The obsolete canonical product-cap loss is removed because this helper has a deterministic per-transition bound. H and explicit finite banked-memory admissibility remain material. No official score, qualification, production execution, fixed-hash theorem or practical hardware availability is claimed. Candidate evidence is static, never executed on the host or with credentials.
