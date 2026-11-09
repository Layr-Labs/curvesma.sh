# Constant-A inversion with a quartic trace projection

Proposed time_log2=12621/100=126.21; conditional success>0.390003; memory<2^138 bytes; cold preprocessing<2^124 scalar instructions INCLUDED; advice0. Full protected p,r, supplied nonidentity P, threshold0.39 and C=20000 are unchanged. This is an analytic bound, not an execution or official score. Material conditions: the full-standard-hash event bridge H below and ordinary bounded word-message storage. Neither is established hardware/production behavior.

Material published sources: e750ade2-0c61-4061-87a7-9f66de49eaae supplies the supplied-P reduction, constant-A eleven-phase collision argument, trace-orbit key, full hash domains and index radix; authentic proof SHA256=1d034090feb99988f518b36b46f84aeabd610f318d1c64ac20fc9fb1806cf5a3. f00694ca-67df-4231-b803-8825a49ef896 supplies finite word-message routing/controller initialization, not BSGS or its score; authentic proof SHA256=eaee7faf9dc4f134f2ea3c65debd0663b61771fab87acc9fc0f6b80348ef2863 at recorded commit6c8029b96e3f9bde62bfcfc31ac1483ff25c5149. Exact source/recipient passages are declared. No qualification/ancestry certification is inherited. New: nine-Fp2-product quartic projection with charged cached sums, native arithmetic and word transfers. Unpublished development is disclosed in the public note.

## Complete algorithm

Decode and validate the full protected input. For z=1 return verified G2 infinity. Otherwise h=pairing(Q0,P) has order r; write z=h^a,a!=0 for analysis. W=2^92,R=min integer with11R^2>=r,L=ceil(10305R/10000)=71148527125680446472998557496245104460. Generate full-Fr v_j,w_k for j,k<W by SCAL below. Cold construct T0[j]=h*z^v_j,T1[k]=z^w_k. A table record960 contains native GT576,zero32,beta32@608,six48-byte component sums@640..927,zero32. Six canonical Fp additions compute the sums of real/imaginary parts of the six Fp2 coordinates. These are computed target-dependent data, not supplied powers/advice.

For t<L,PAIR yields j,k. Copy both full records to coordinator buffers,add beta=v_j+w_k mod r. project_quartic.S computes ONLY a0,b1 of conceptual X=T0*T1=h*z^beta. Abort if b1=0. Compute exact S192 and fingerprint U32 below. Immutable payload256=(S192,beta32,U32),metadata64=(U32,t128,zero16). Two stable80-bit radix passes sort metadata by all160 fingerprint bits. Enumerate equal-fingerprint unordered pairs,aborting BEFORE more than floor(L/128) confirmations. Fetch full payloads by exact ordinals. Reject unequal S and equal beta. For FIRST equal S/different beta pair,try eleven lambda=(p mod r)^ell,ell1..11. For D=beta_t-lambda*beta_s nonzero,test a=(lambda-1)/D by h^a=z. True phase has1-lambda+aD=0 and D!=0. Equal beta would give equal raw X;different beta excludes phase0. Return Q=[a]Q0 only after full canonical G2/curve/subgroup and pairing(Q,P)=z verification for the actual P. At most eleven terminal attempts TOTAL. No walk or mixing premise.

## Actual hashes and explicit heuristic

Full seven-round BLAKE3-256,primary immutable definition https://github.com/BLAKE3-team/BLAKE3/blob/f55849f89cd85083c1c9daa2c5de20766f309c00/reference_impl/reference_impl.rs .16-byte ASCII tags FAPI-RHO-INIT-v1,FAPI-RHO-SCAL-v1,FAPI-RHO-PAIR-v1,FAPI-RHO-FING-v1. seed=B3(INIT||P96||z576),canonical48-byte BE coordinates in protected order:688 bytes,eleven blocks,first64 flags1,middle0,last48 flags10,zero padding,counters0. SCAL hashes tag||seed32||counter128_BE16,length64,flags11,counter0;digestBE low255 accepted iff<r. Increment global counter on EVERY attempt,cap1024 per scalar.2W accepts need<=2^103 attempts. PAIR hashes tag||seed32||t128_BE16,length64,flags11,counter0;digestLE low92 j,next92 k. FING same64-byte form,separate increment-on-every-attempt counter,LE low61 accepted iff<m=2^61-1,1024 cap each,81 coefficients cold. All actual calls are computed.

compression.S implements standard CV8||IVfirst4||counter2||length||flags and all seven message schedules. G is16 instructions:two message loads plus ADD/XOR/ROR,rotations16,12,8,7.17 header+896 rounds+15 final=928 throughRET,929 including BL;retain960 inclusive. Intermediate outputs chain CV;ROOT only on final block. Encoding/control is separately main/cold paid.

In the ideal experiment domain outputs are independent uniform,consistent on repeats;SCAL/FING rejection gives uniform accepted values. H is the heuristic that the ACTUAL full standard domains,averaged over uniform full-production P,z,support the finite collision moments,rejection bounds and unequal-key expectation below. Full rounds/domain separation support plausibility;no small test proves this enormous-scale extrapolation. This is NOT literal joint IID or total-variation equivalence for a fixed256-bit seed. If H is unsupported,the conditional success claim fails. No coin device is assumed;otherwise no fresh coins. This version uses only the moment event bridge,not the stronger no-collision bridge of our unpublished dependency draft.

## Exact orbit key

Fp2=Fp[u]/(u^2+1),Fp6=Fp2[v]/(v^3-xi),xi=1+u,Fp12=Fp6[w]/(w^2-v);canonical native Montgomery scale2^384. Protected polynomial conversion u=w^6-1 is cold paid. X=a0+a1w^2+a2w^4+b0w+b1w^3+b2w^5. X,X^(p^4),X^(p^8) have product1,trace c=3(a0+b1w^3),second symmetric coefficient c^(p^2);they are exactly the roots of Y^3-cY^2+c^(p^2)Y-1. Thus the trace identifies exactly three points (XTR,primary https://www.cs.ru.nl/E.Verheul/papers/crypto2000/crypto2000.pdf ). p^2 fixes a0 and negates b1. Choose between b1,-b1 by first nonzero canonical Fp component:negate BOTH iff that component>=ceil(p/2). Copy a0,b1 without multiplying by3. For nonidentity GT,b1!=0:otherwise the cubic factors(Y-1)(Y^2+(1-c)Y+1) over Fp2,so a nonidentity root is in Fp4,contradicting embedding degree12. F=(a0,canon_pm b1) has six even-Frobenius fibers.

Cold compute ACTUAL C0=xi^((p-1)/2)=(s,s),with its actual sign and2s^2=-1. C0^2=-u and p mod8=3 exclude opposite components since2 is nonsquare. Odd Frobenius conjugates a0 and sends B=Br+Bi*u to(s(Br+Bi),s(Br-Bi)),TWO Fp products. Canonicalize again. S=min(F(X),F(X^p)) has exactly12 fibers;odd/even classes distinct,p swaps them. Full24-word key,not truncation. Budget:projection/sign192+conjugation128+sum/diff256+2 products800+second sign/interface256+full comparison128+copy24+other pointers/frame/constants256=2040;retain2304. First projection/sign helper97+wrapper<=52+BL fits192;full six-limb adds/subs include one p correction,loads/stores/carries and fit coarse128. Identity/zero paths are included.

## Nine-product projection and native bounds

Set t=w^3,t^2=xi. Over Fp4=Fp2[t],U0=a0+b1t,U1=b0+a2t,U2=a1+b2t;V0=c0+d1t,V1=d0+c2t,V2=c1+d2t. EXACTLY

    Z=U0V0+t(U1V2+U2V1)=a0(XY)+b1(XY)t.

Each quartic product uses ac=A*C,bd=B*D,cross=(A+B)*(C+D),real=ac+xi*bd,imag=cross-ac-bd. Thus NINE Fp2 products. Six hot Fp2 sums cost12 Fp adds;their real+imaginary sums cost6. Three recombinations cost24 Fp ops;final sum/xi/output costs10;total52 Fp ops. xi*(x+y*u)=(x-y)+(x+y)u. All raw-coordinate component sums are cold paid. This is ring algebra for arbitrary full-width coefficients,not a statistical approximation. Only a0@0,b1@384 are written in output576;holes are NEVER read as GT. Conceptual X is GT since table inputs are GT.

fp381_add.S/sub.S have31 instructions each throughRET,canonical full381-bit six-limb results. Add:6 input LDP,6 carry adds,3 modulus LDP,6 trial subtraction,6 CSEL,3 STP,RET;2p<2^384. Sub:6 LDP,6 subtract,3 modulus LDP,6 select p/0 under ORIGINAL borrow,6 adds,3 STP,RET;all selections precede flag changes. No callee/x18/SP change,reads before writes. Four arguments+BL<=36;retain40.

Scalar multiply:EXACT pinned blst mul_mont_384 wrapper29+straight-line core341=370,including internalBL/bothRETs. Primary https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/build/elf/mul_mont_384-armv8.S ,full SHA256=a8234680738f79d8f438c6b27832b209d14fe610cc83bfbd37bef9f74213a2b7. Canonical six-limb inputs give ab/2^384 modp,canonical same native scale;standard Montgomery T<2p/carry bound and one correction apply. Callee registers restored,x18 untouched. Caller<=10;retain400 for key products. Protected383 profiling corroborates,not a changed primitive price.

native_fp2_pre.S accepts actual canonical real+imaginary sums;3 pinned products+3 canonical subtracts implement Karatsuba.224-byte aligned frame:80 saved registers,3*48 temporaries. Own52 includes all internal calls/setups;52+3*370+3*31+7 argument instructions+BL=1263;retain1280 inclusive. All input/sum reads finish before output stores matter.

project_quartic.S frame1808:80 save,864 products,576 hot sums,288 sum components. x24=sp+80,x25=sp+944,x26=sp+1520;raw coordinate q uses component sum input+640+(q/96)*48. Three groups of(ac,bd,cross) produce imag in cross FIRST,then real in ac using dead hot-sum scratch for xi*bd. Final Z uses those real/imag slots;no live operand is overwritten. Native source gives353 own instructions,including9*8+52*5=332 setup/BL already covered by inclusive budgets. Remaining21+outer<=10 fit256.9*1280+52*40+256=13856;retain13952. Full loads/stores/frames/corrections paid. Cold/terminal full operations retain C=20000. Native evidence is never host-executed.

## Finite success, fingerprint and aborts

Fix a!=0. In ideal analysis C_j=1+a*v_j,E_k=a*w_k independent uniform Fr,x_t=C_Jt+E_Kt. Distinct types give uniform independent raw pairs even sharing one index. Let B_i=binom(L,i),C count distinct-type nonzero raw pairs related by11 nonidentity phases. Then

    mu=B_2*(1-W^-2)*11*(r-1)/r^2,
    E binom(C,2)<=U2=3*121*(B_3*(r^-2+1/(rW^2))
                                     +B_4*(r^-2+4/(rW^4))).

Two rows in labeled C/E have coefficients1,-lambda,1-lambda,all nonzero. Rank2 costs<=r^-2. Rank1 requires matching support in BOTH families:three times force other endpoints identical (W^-2);four times permit two crossed matchings/family (4/W^4). Singleton support forces all four equal and is covered. Three shared/disjoint shapes and121 phase pairs give U2. No empty support because phase0 excluded. Pointwise1[C>0]>=C-binom(C,2) gives probability>=mu-U2,without independent collisions or Poisson approximation.

Fingerprint:prime m=2^61-1 (Lucas-Lehmer59 steps from4 gives0,exponent61 prime;primary https://isa-afp.org/entries/Mersenne_Primes.html ). Entire S192 injects into25 consecutive60-bit LE chunks and final36-bit chunk x_i<m. Three independent rows of27 coefficients compute y=b+sum a_i*x_i modm. Distinct keys give rank2 onto F_m^2,so two evaluations uniform independent per row. Low54,54,52 bits pack160;remaining12 U bytes zero. Bucket size<=2^(61-b) gives unequal-key probability<=delta=2^23/m^3. Pairwise expectation suffices,conditional on all keys. Standard universal hashing credit: https://research.ibm.com/publications/universal-classes-of-hash-functions--1 .

mersenne_fingerprint.S:3 deferred128-bit dot products;26 shared chunk extractions3+3*(load,MUL,UMULH,ADDS,ADC)=18 instructions/term. Sums<2^126. Split61/61/6 bits,sum<=2m+15,fold once<=m+2,one SUBS/CSEL.7 init+468 terms+33 reductions/masks+4 packing+2STP+RET=515;three args+BL519;retain640. Exact160-bit packing,word23@184 last load,no overread. Caller-saved only,x18/SP/callees untouched,all stores after reads.

All equal-S pairs have expectation<=B_2*(W^-2+12/r),unequal-S fingerprint pairs<=B_2*delta. Every confirmation is one of these. Cap loss<=B_2*(W^-2+12/r+delta)/(floor(L/128)+1). Identity abort<=L/r. SCAL per-attempt rejection<1/8:2^93 capped draws fail<2^-2979;FING81 caps fail<=81*2^-62464;combined<2^-970. No failure independence needed. Exact rational evaluation gives

    mu-U2-caploss-L/r-2^-970 >390003/1000000 >0.39.

Diagnostic0.3900031399012319 is not the certificate. H transfers these finite bounds to actual average uniform P,z;identity target is solved directly.

## Ordinary word storage and complete work

Tree fanouts2^27,2^35,2^35,depth3,2^97 leaves,<2^98 processors. Each leaf1-TiB data RAM at2^40,code/private below2^40,parent port2^42,child j port2^43+(j<<11);all pointers<2^47.97-bit bank TWO64-bit words;child IDs bank>>70 (high>>6),low35(bank>>35)(EXTR/AND),low35bank. No global137-bit pointer or unit-cost wide arithmetic.

Link endpoints2KiB,request/reply<=1KiB. Four-word header banklo/hi,localoffset40,op/length;max120-word payload gives992 bytes leaving32 control. Reliable blocking ready-word protocol,no busy spins. Links only move fixed words,not DMA/search/arithmetic/entropy. Leaf count/next use paid LDP,ADDS/ADC,STP. One RPC globally in flight;sum ALL processors. Communication admissibility is material,not a physical certificate.

Roundtrip3 request+3 reply hops. Paired senderSTP+receiverLDP costs2 instructions/two words/hop,thus payload3N;headers8*3=24. Each directed hop has64 OTHER sender+receiver instructions:8 save,8 decode,8 route/address,8 dispatch,12 handshake/barriers,12 restore/control,8 spare. Fixed unrolled N={0,2,8,32,120}. R(N)=3N+408. Root construction64 replaces old32 selector per RPC;endpoint load/store separate.

Banks/offsets:metadata buffer b=0,1,index i<2^127:bank b*2^93+(i>>34),offset low34(i)*64. Payload:bank2^94+(i>>32),offset low32(i)*256. Table b,j<W:bank2^96+b*2^62+(j>>30),offset low30(j)*960. Count/next b,d<2^80:bank2^96+2^64+b*2^44+(d>>36),offset low36(d)*16. Disjoint bank ranges<2^97,records fit1TiB,local data pointers<2^41. All used bytes initialized;unused capacity supplies no data/advice.

Radix:zero all2^80 two-limb counts,scan metadata counting selected80-bit digit,prefix ALL digits into separate next array,stable scan/scatter ENTIRE64-byte records;swap buffers,low80 then high80. Digits use64+16 bits at0/8 or10/18,unaligned fixed-word load permitted.128-bit ordinals/counters always two limbs. Payload never sorted. Count/scatter endpoint96+176 per record/pass=544L;scan128L. Full24-word S and4-word beta comparisons and all losing paths charged.

Generation key2304+projection13952+finger640+PAIR960+Fr-add256+main640=18752. Main covers rawT0 copy72/scalar4,args48,PAIR encoding32,fingerprint pack48,payload/meta construction/copies128,indices80,frames64,remaining root bookkeeping/margin. Add radix544,scan128,confirm32=19456;retain19680 with224 spare. Old32 selectors included and replaced below.

Generation two READ120,WRITE32,WRITE8:Rsum2472+four selector replacements128+extra T1 raw endpoint80+extra T0 padding4+two added40-word portions endpoint80=2764L. Per radix pass count READ8+INC_ACK0:840+64 selector replacement+2 digit extraction=906;scatter READ8+FETCH/INC_NEXT2+WRITE8:1278+96=1374. Two passes4560L. Scan READ8:432+32+8 endpoint margin=472L. Network total7796L;19680+7796=27476,retain27520. Every loaded/stored cached sum/padding byte is paid.

Confirmation<=4096:512 fullS/beta/ordinal/control+1872 transport(twoREAD8,twoREAD32)+128 selector replacements+EXTRA256 endpoint/padding=2768. Cap floor(L/128) gives<=32L already included. First useful pair triggers only11 terminal attempts.

Cold init<2^98 processors*2^24 instructions/node=2^122,including delivery/bootstrap of<=1MiB code/private/page-table image.1GiB block mappings need1024 descriptors/leaf RAM and<=65536/internal ports,fit image. Separately all4KiB/link flags/endpoints initialized within2^118. Each2^93 table entry<=512 dense references*20000+1024 SCAL attempts*8192+65536 scalar/power/record/six sum additions/WRITE120/control=18694144<2^25;all tables<2^118. ALL zero/prefix RPCs<2^96. Seed/FING,decoding/conversion/subgroups,h pairing,C0/11phases and ALL11 scalar/G2/pairing/output attempts<2^92;fixed exponents<4608bits give this loose cap. Sum setup2^122+links2^118+tables2^118+radix2^96+fixed2^92<2^124,no amortization.

Peak leaf capacity2^137+private/page/code<2^118+link endpoints<2^110+root<2^40 gives<2^138 bytes. Advice0. N=27520L+2^124;EXACT integer N^100<20000^100*2^12621 proves submitted126.21. No production run,empirical success,physical installation or qualification is claimed. Experiments empty for this analytic variant;shared llm-verifier alone owns execution/review/scoring.
