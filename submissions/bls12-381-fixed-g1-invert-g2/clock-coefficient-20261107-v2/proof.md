# Unrotated independent-increment walk with full twelve-orbit detection

Proposed time_log2=125843/1000=125.843: complete computation<2^(125843/1000) protected dense-Fp12 equivalents at C20000,success>=0.39 under the disclosed full standard BLAKE3 bridge H and explicit finite banked-memory admissibility. Peak memory<2^136 bytes,zero advice. All per-target cold construction,scalar attempts,initialization,transfers,full hashes,dictionary operations,failures,verified extraction and canonical G2 output are charged. These are conditional analytic bounds, not qualification,promotion,an official score,a fixed-hash theorem,practical attack execution or hardware installation.

## Exact supplied problem

Use the full protected BLS12-381 p,r,Q0 and the exact supplied canonical nonidentity P. Validate P and GT z. If z=1 output G2 infinity. Otherwise h=pairing(Q0,P) generates GT for every allowed P by prime-order nondegeneracy. Solve z=h^a,verify h^a=z,return canonical Q=[a]Q0 only after curve/subgroup and pairing(Q,P)=z checks. No replacement P,hidden exponent,auxiliary powers,scalar-only answer or unchecked result is used.

Maintain raw X=h^A*z^B in full576-byte native GT plus canonical Fr histories A,B. fp4-trace.md proves that S=Tr_(Fp12/Fp4)(X)/3=a0+b1*w^3 has exact fibers{X,sigma4X,sigma8X}; canonicalizing its four Frobenius images identifies exactly the full12 orbit. The key is192 native bytes,not a free group operand or unproved projection. Canonical orientation k0..3 is stored for terminal recovery only. Every visited orbit is stored globally and first repeat detected directly; no post-repeat quotient-map law is assumed.

## Paid cold construction and standard hashes

Let W=2^76. Sample independent uniform Fr v_j and beta0. Construct one raw jump M_j=h*z^v_j and store it in a640-byte record with canonical fixed u=1 and v_j. Initialize X0=h*z^beta0,A0=1,B0=beta0. Each accepted transition gives A_t=t+1,not a random history sum. Only B needs modular Fr addition. Since z!=1 is a prime-order generator,the raw jumps and initial raw X are independently uniform GT;mixing.md proves this explicitly. No Frobenius jump variants are stored. Terminal Frobenius constants are constructed once cold. Every power,record,coefficient,write,rejected scalar attempt and constant is target-specific paid work,not advice or amortization.

B3 is standard unkeyed BLAKE3-256,all seven rounds. INIT,SCAL,JUMP are exact16-byte ASCII FAPI-RHO-INIT-v1,FAPI-RHO-SCAL-v1,FAPI-RHO-JUMP-v1. seed=B3(INIT||P||z),P canonical x||y each48-byte big-endian,z twelve canonical48-byte Fp coefficients in pinned order. INIT is688 bytes: eleven blocks,first flags1,intermediate0,final48 flags10,counter0. SCAL uses a global unsigned128-bit counter encoded16-byte big-endian,advanced each attempt. Its64-byte message SCAL||seed||counter has len64,flags11,counter0. Interpret digest big-endian,take low255 bits,accept iff<r,cap1024 tries per scalar. AllW+1<2^77 scalars and<2^88 attempts fit the counter; r>7*2^252 bounds scalar-cap loss<2^-970.

nh-fingerprint.md specifies65 new full standard B3 seed calls in separate NHK/NHC domains. They initialize a1040-byte record for five independent NH128 lanes and independent odd-multiplier compression to five64-bit outputs. Its full scalar implementation tabulate.S costs460. Serialize the resulting U320 as40 big-endian bytes. Fixed-pair collision<2^-312; independent ideal-path union gives coupling failure<2^-61,retained2^-55. Actual values are computed/charged. No simple-tabulation table or constant-zero digit remains.

JUMP=B3(JUMP||seed||zero16||U40) is104 bytes. Prefix64 len64 flags1 is compressed cold;prepared-hash.md proves its eleven-instruction fixed-prefix preparation preserves the full standard digest. Final40 uses counter0,len40,flags10 and all56 G/seven rounds,776 scalar instructions. Read digest Y little-endian. Set j=Y mod2^76,bucket=(Y>>76)mod2^128. These disjoint bits are exactly independent uniform under H. No phase,orientation correction,previous edge,inverse guard,fictitious generator,retry hash,bias allowance or per-step rejection exists.

Initialize the raw state and histories as above,compute canonical key,k,U and initial label/bucket,aborting on identity;insert its initialized payload. For at most K=ceil(1243R/1250) accepted transitions,R=min{R:12R^2>=r},multiply raw X by one selected M[j]. The accepted-step counter t is passed as two64-bit words to scalar-and-copy.S:it writes canonical A=t+1 as four limbs,upper two zero because K+1<2^126<r,and updates all four B limbs modulo full r. Recompute key,k,U,label/bucket and query the global dictionary. Every new orbit is inserted;the first repeat is detected immediately. All source ABIs,initializations and transfers are paid.

## Extraction and finite success

On equal U,compare all24 exact native key words; mismatch is fingerprint failure. On equal key,reconstruct saved raw Xs=h^As*z^Bs in paid terminal cold work. Test ell=ks-kc+4m mod12,m0..2 against current full X. Exactly one matches because nonidentity GT has free12 orbits. Form D=Bcur-xi^ell*Bs. D0 fails; otherwise a=(xi^ell*As-Acur)/D,verify h^a=z,form canonical Q and verify exact pairing/curve/subgroup before output. Reconstruction returns its group invariant,so the true twist is group-path determined,not selected by scalar freedom.

mixing.md proves stationary independent-increment Fourier kernels for W unrotated random jumps,with delta2^-33 and good-table failure<2^-98. success.md proves finite even-sixth Bonferroni repeat probability>0.3900538335845625,including all short correlated gaps and rotated orbit constraints. probability.md proves a complete full-size Frobenius interval lemma using U=0xd201000000010000,r=U^4-U^2+1 and xi=-U modr. Among1<=s<t<=K+1,the only equality t=xi^ell*s is ell7,t=U*s. A singular extraction denominator would force such an equality between known A values. There are at most floor((K+1)/U) exceptional time pairs,all at lags>=U-1>T128. Pointwise good-table mixing and a union bound give singular-match probability<2^-180. Every other genuine match has nonzero D deterministically. No uniform-denominator,coefficient-count tail,formal norm,random-alpha or relation heuristic is used.

Subtract fingerprint2^-55,bad table2^-98,scalar cap2^-970,identity2^-100,dictionary cap2^-100,and singular2^-180. Complete success>0.39005383358456>0.39. Losses use a union bound,not independence. H remains an explicit unproved cryptographic bridge at full query scale. No random-mapping,Poisson,pair-independence,trace-distribution or heuristic mixing assumption is added.

## Complete computation and physical capacity

Scalar ordinary work<=1408canonical+460fingerprint+776hash+30Frupdates+72controller+192head+224payload+32cappedchain=3194,retained3208 with14 spare instructions. Controller has50 main+14 selector instructions,retained72; every helper body/ABI is separately charged. There are exactly K protected reference multiplications charged once. Complete per-target coldD=540*2^119 includes all2^128 initialized head zeros,bank/index setup,group/scalar/variant/seed/code/constants/program construction,initial state/hash/payload,terminal raw reconstruction/three-twist search,failures,inversion,canonical G2 output and exact verification. Every allocated raw jump is initialized/paid.

Work <=K+(3208K+540*2^119)/20000=(23208K+540*2^119)/20000.
Exact integer certificate(23208K+540*2^119)^1000<20000^1000*2^125843 gives proposed125.843,above diagnostic125.84292322111403. Reference prices and full parameters are unchanged.

Physical capacity is2^126*512 payload+2^128*16 heads+2^46*2^40 jump banks+shared<2^45,below2^136 bytes including all padding. storage.md gives97-bit logical bank selection by ordinary64-bit latches/commit/barriers and native48-bit window addresses. All selector arithmetic/transfers are charged; no flat136-bit address or free search/allocator exists. If this explicit finite banked-memory model is inadmissible,the claim must be rejected. It is not installed on the VM. Analytic evidence is allowed and local intake does not allocate production attack memory.

## Sources and status

Published23d3d4cd-768d-4cac-a9a5-872fc010580d supplies the supplied-P reduction and finite Fourier/Bonferroni/cold-accounting framework. Published13f0490c-4535-4ee8-8112-6a72a5254024 supplies the full BLAKE3/helper/Fr contracts; its AI-valid127.756 was unpromoted and does not qualify this revision. Direct material local predecessor201ab5a6-9f73-4789-bea0-eff2855e3b43 supplies the exact native field ledger,NH/trace/hash implementation and earlier finite moment framework. Its unpublished formal link is deferred,not its material attribution omitted. Other local graph/hash/trace predecessors and external NH/XTR mathematical sources are credited in the public note and dedicated evidence. No novelty,independence or qualification is inferred from operational review failures.

Pinned primary definitions:
https://github.com/BLAKE3-team/BLAKE3/blob/f55849f89cd85083c1c9daa2c5de20766f309c00/reference_impl/reference_impl.rs
https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/src/fp12_tower.c

Candidate evidence is untrusted static text,never compiled/executed on the host or with credentials. Shared llm-verifier owns scientific execution/review/scoring; experiments are empty. Root submissions remain add-only. All six production tracks,one computation metric,full sizes,success0.39 and reference prices are preserved.
