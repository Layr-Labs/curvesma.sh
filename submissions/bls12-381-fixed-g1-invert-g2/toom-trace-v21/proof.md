# Twenty-one scalar products through exact interpolation

Proposed time_log2=3139/25=125.56; conditional success>0.390003; memory<2^137 bytes; cold preprocessing<2^120 scalar instructions INCLUDED; advice0. Full protected p,r, supplied nonidentity P, threshold0.39 and C=20000 are unchanged. This is an analytic bound, not an execution or official score. Material conditions: the full-standard-hash event bridge H below and finite bank-window storage B. Neither is established hardware/production behavior.

Material published sources: e750ade2-0c61-4061-87a7-9f66de49eaae supplies the supplied-P reduction, constant-A eleven-phase collision argument, trace-orbit key, full hash domains and index radix; authentic proof SHA256=1d034090feb99988f518b36b46f84aeabd610f318d1c64ac20fc9fb1806cf5a3. The source finite97-bit bank-window model and32-instruction selectors are reused under explicit B;the discarded word-message routing architecture supplies no cost or qualification. Exact source/recipient passages are declared. No qualification/ancestry certification is inherited. New: exact degree3 evaluation/interpolation for the trace projection,21 scalar products,all weighted evaluations cold paid,with explicit native recombination. Standard interpolation is not claimed as invented. Unpublished development is disclosed in the public note.

## Complete algorithm

Decode and validate the full protected input. For z=1 return verified G2 infinity. Otherwise h=pairing(Q0,P) has order r; write z=h^a,a!=0 for analysis. W=2^92,R=min integer with11R^2>=r,L=ceil(10305R/10000)=71148527125680446472998557496245104460. Generate full-Fr v_j,w_k for j,k<W by SCAL below. Cold construct T0[j]=h*z^v_j,T1[k]=z^w_k. For EACH GT point compute the21 full-Fp native evaluations specified below. A record1088 contains21*48=1008 evaluation bytes,canonical scalar32@1008,zero48. T0 holds unweighted evaluations,T1 weighted evaluations. GT points need not remain stored after this exact transformation. Every evaluation/inverse-weight multiplication/byte is cold paid,target-dependent data,not supplied powers or advice.

For t<L,PAIR yields j,k. Copy T0 evaluations1008 and scalar32 into shared scratch (130 CPUcopy instructions);select T1 and keep its bank committed through native projection,then add beta=v_j+w_k mod r. project_toom.S computes ONLY a0,b1 of conceptual X=T0*T1=h*z^beta. Abort if b1=0. Compute exact S192 and fingerprint U32 below. Immutable payload256=(S192,beta32,U32),metadata64=(U32,t128,zero16). Two stable80-bit radix passes sort metadata by all160 fingerprint bits. Enumerate equal-fingerprint unordered pairs,aborting BEFORE more than floor(L/128) confirmations. Fetch full payloads by exact ordinals. Reject unequal S and equal beta. For FIRST equal S/different beta pair,try eleven lambda=(p mod r)^ell,ell1..11. For D=beta_t-lambda*beta_s nonzero,test a=(lambda-1)/D by h^a=z. True phase has1-lambda+aD=0 and D!=0. Equal beta would give equal raw X;different beta excludes phase0. Return Q=[a]Q0 only after full canonical G2/curve/subgroup and pairing(Q,P)=z verification for the actual P. At most eleven terminal attempts TOTAL. No walk or mixing premise.

## Actual hashes and explicit heuristic

Full seven-round BLAKE3-256,primary immutable definition https://github.com/BLAKE3-team/BLAKE3/blob/f55849f89cd85083c1c9daa2c5de20766f309c00/reference_impl/reference_impl.rs .16-byte ASCII tags FAPI-RHO-INIT-v1,FAPI-RHO-SCAL-v1,FAPI-RHO-PAIR-v1,FAPI-RHO-FING-v1. seed=B3(INIT||P96||z576),canonical48-byte BE coordinates in protected order:688 bytes,eleven blocks,first64 flags1,middle0,last48 flags10,zero padding,counters0. SCAL hashes tag||seed32||counter128_BE16,length64,flags11,counter0;digestBE low255 accepted iff<r. Increment global counter on EVERY attempt,cap1024 per scalar.2W accepts need<=2^103 attempts. PAIR hashes tag||seed32||t128_BE16,length64,flags11,counter0;digestLE low92 j,next92 k. FING same64-byte form,separate increment-on-every-attempt counter,LE low61 accepted iff<m=2^61-1,1024 cap each,81 coefficients cold. All actual calls are computed.

compression.S implements standard CV8||IVfirst4||counter2||length||flags and all seven message schedules. G is16 instructions:two message loads plus ADD/XOR/ROR,rotations16,12,8,7.17 header+896 rounds+15 final=928 throughRET,929 including BL;retain960 inclusive. Intermediate outputs chain CV;ROOT only on final block. Encoding/control is separately main/cold paid.

In the ideal experiment domain outputs are independent uniform,consistent on repeats;SCAL/FING rejection gives uniform accepted values. H is the heuristic that the ACTUAL full standard domains,averaged over uniform full-production P,z,support the finite collision moments,rejection bounds and unequal-key expectation below. Full rounds/domain separation support plausibility;no small test proves this enormous-scale extrapolation. This is NOT literal joint IID or total-variation equivalence for a fixed256-bit seed. If H is unsupported,the conditional success claim fails. No coin device is assumed;otherwise no fresh coins. This version uses only the moment event bridge,not the stronger no-collision bridge of our unpublished dependency draft.

## Exact orbit key

Fp2=Fp[u]/(u^2+1),Fp6=Fp2[v]/(v^3-xi),xi=1+u,Fp12=Fp6[w]/(w^2-v);canonical native Montgomery scale2^384. Protected polynomial conversion u=w^6-1 is cold paid. X=a0+a1w^2+a2w^4+b0w+b1w^3+b2w^5. X,X^(p^4),X^(p^8) have product1,trace c=3(a0+b1w^3),second symmetric coefficient c^(p^2);they are exactly the roots of Y^3-cY^2+c^(p^2)Y-1. Thus the trace identifies exactly three points (XTR,primary https://www.cs.ru.nl/E.Verheul/papers/crypto2000/crypto2000.pdf ). p^2 fixes a0 and negates b1. Choose between b1,-b1 by first nonzero canonical Fp component:negate BOTH iff that component>=ceil(p/2). Copy a0,b1 without multiplying by3. For nonidentity GT,b1!=0:otherwise the cubic factors(Y-1)(Y^2+(1-c)Y+1) over Fp2,so a nonidentity root is in Fp4,contradicting embedding degree12. F=(a0,canon_pm b1) has six even-Frobenius fibers.

Cold compute ACTUAL C0=xi^((p-1)/2)=(s,s),with its actual sign and2s^2=-1. C0^2=-u and p mod8=3 exclude opposite components since2 is nonsquare. Odd Frobenius conjugates a0 and sends B=Br+Bi*u to(s(Br+Bi),s(Br-Bi)),TWO Fp products. Canonicalize again. S=min(F(X),F(X^p)) has exactly12 fibers;odd/even classes distinct,p swaps them. Full24-word key,not truncation. Budget:projection/sign192+conjugation128+sum/diff256+2 products800+second sign/interface256+full comparison128+copy24+other pointers/frame/constants256=2040;retain2304. First projection/sign helper97+wrapper<=52+BL fits192;full six-limb adds/subs include one p correction,loads/stores/carries and fit coarse128. Identity/zero paths are included.

## Exact evaluation/interpolation projection

Let t=w^3,t^2=xi=1+u,so q(t)=t^4-2t^2+2=0. Over Fp4=Fp2[t],U0=a0+b1t,U1=b0+a2t,U2=a1+b2t,and V analogously. EXACTLY

    Z=U0V0+t(U1V2+U2V1)=a0(XY)+b1(XY)t.

Convert A+B*t,with A=Ar+Ai*u,B=Br+Bi*u,into the degree3 polynomial(Ar-Ai)+(Br-Bi)t+Ai*t^2+Bi*t^3. Thus three products of degree<=3 polynomials suffice. For six finite points x=(0,1,-1,2,-2,3),let V(T)=product(T-x_i),N_i(T)=V(T)/(T-x_i),d_i=N_i(x_i)=(-12,12,24,-24,-120,120). For ANY polynomial f degree<=6,

    f(T)=sum_{i=0}^5 f(x_i)*N_i(T)/d_i + f_6*V(T).

Both sides have the same leading coefficient and all six roots of their difference;the remaining degree<=5 difference is zero. All d_i invertible at protected full p. The infinity evaluation is simply the leading coefficient. Reduce N_i,V modulo q,then map a reduced vector(c0,c1,c2,c3) to(a0r,a0i,b1r,b1i)=(c0+c2,c2,c1+c3,c3). For t times the reduced vector,the map is(c1,c1+2c3,c0+c2,c2). This derives the following exact columns,not empirical fitting. Point order is0,1,-1,2,-2,3,infinity;d_infinity=1.

    plain = [(3,9,-1,-3),(8,4,5,-5),(16,8,-13,1),
             (1,-1,-1,-5),(5,-5,-1,7),(0,0,-1,-3),(2,-4,3,9)]
    g     = [1,2,2,1,1,2,2]
    shift = [(2,-4,3,9),(5,0,4,2),(-7,-6,8,4),
             (4,-6,1,-1),(-8,6,5,-5),(1,-2,0,0),(-3,6,1,-2)]

For T0 store seven unweighted evaluations of EACH U0,U1,U2. For T1 store V0(x_i)/d_i,and g_i*V1(x_i)/d_i,g_i*V2(x_i)/d_i;at infinity use leading coefficients with d=1. All quantities are full canonical Fp in the same native Montgomery scale2^384. Cold scalar multiplication by public rational weights is fully charged. Define P_i=U0(x_i)V0(x_i)/d_i and Q_i=g_i*(U1(x_i)V2(x_i)+U2(x_i)V1(x_i))/d_i,with leading products at infinity. Then EXACTLY

    (a0r,a0i,b1r,b1i)=sum plain_i*P_i+sum shift_i*Q_i modp.

The division of shifted columns by g is compensated by the paid T1 weights. No full GT multiplication is performed hot. Twenty-one scalar products produce P0,P1,P2;seven field additions combine the latter two into Q. The displayed integer columns require61 positive-multiple additions (binary addition chains with shared powers per variable) and47 accumulator adds/subtracts. Four initial field copies cost7 instructions each. All signs,wraps and normalization use full canonical p arithmetic. This identity holds for arbitrary full-width coefficients;for computed GT inputs it gives the SAME conceptual-product a0,b1 and hence the SAME S and success events. Output576 writes ONLY a0@0,b1@384;holes are NEVER read as GT.

## Complete native implementation

fp381_add.S/sub.S each have31 scalar instructions throughRET,canonical six-limb full381-bit arithmetic. Add uses6inputLDP,6carry adds,3modulusLDP,6trial subtract,6CSEL,3STP,RET;2p<2^384. Sub uses6inputLDP,6subtract,3modulusLDP,6p/0selections under ORIGINAL borrow,6adds,3STP,RET. All selections precede flag changes. Reads precede stores;no SP/x18/callee changes. Five caller setup/BL instructions make36 inclusive for each call here.

pinned_mont_core.S is the EXACT341-instruction __mul_mont_384 core at blst commitf7a88067ed98385ebc6cf41ccecf21ead706a5b7,primary https://github.com/supranational/blst/blob/f7a88067ed98385ebc6cf41ccecf21ead706a5b7/build/elf/mul_mont_384-armv8.S ,full original SHA256=a8234680738f79d8f438c6b27832b209d14fe610cc83bfbd37bef9f74213a2b7. Function instructions unchanged,global binding added;Apache licence included. It computes canonical ab/2^384 modp in the same native scale. Standard Montgomery full carry/T<2p and one-correction bounds apply. The ordinary29-instruction scalar wrapper is also included;key products29+341+caller<=400,not a changed C or primitive price.

The private core expects A inx11..16,B0inx17,Bpointerx2,pinx5..10,n0inx4,and context n0/out atx29+96/104. It changes neither SP,x18 norx29. Modulus registersx5..10 are NEVER written;its final x4 write reloads n0. Therefore load p/n0 ONCE before the21 uninterrupted calls. Before EACH call load the original input registers and write its context slots;all341 instructions execute. This hoisting follows the actual pinned body,not a public ABI promise or free multiplication oracle.

project_toom.S frame4128:callee savesx19..28,x29/LR@0..95;core n0/out@96/104;metadata(out,A,B,p,n0)@128..167;products21*48@192..1199;61 multiple temporaries@1200..4127. SP adjusts by4096+32 using TWO legal immediate instructions;restores similarly. All used slots initialized before reading. Products are ordered U0V0,U1V2,U2V1,seven each;A offsets48i,B offsets48i,48(14+i),48(7+i),respectively. Afterwards Q overwrites the second product block;remaining products and fresh multiple slots feed the explicit straight-line field circuit. All input loads finish before output writes;alias inputs are safe. The actual algorithm preserves table values and uses shared output scratch.

Macro FMUL expands13 setup/BL/store instructions;21 calls. FADD/FSUB expand5;115 calls (7Q+61multiples+47output sums). FCOPY expands7;four calls. Prologue12,one p/n0 load4,pointer reload4,epilogue9 give905 own instructions. Nested21*341+115*31 give11631 throughRET;outer arguments/BL<=10 give11641. Retain11776 inclusive. Frame/register/memory effects and EVERY load/store/call/carry correction are included. No runtime interpolation loop or unspecified vector instruction is assumed. Scalar key products and cold/full GT operations retain their original costs.

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

## Finite bank-window model B and complete work

B is EXACTLY the finite storage-address mechanism of public source e750ade2,whose admissibility is material. It is not an installed VM,free search/DMA/arithmetic/entropy oracle or physical hardware certification. A noncacheable2^40-byte window starts2^46,with code/shared scratch below2^46 and all native pointers<2^47. A97-bit bank uses TWO ordinary64-bit latches plus commit. Each switch performs DSB/low/high/commit/DSB;all two-limb mask/shift/offset/interface/barrier work fits32 scalar instructions,as explicitly retained per selector. The controller ONLY selects passive RAM;CPU instructions still perform all loads/stores,arithmetic and copies. No unit-cost137-bit pointer,atomic128 or processor swarm is assumed. If B is inadmissible,the conditional bound fails.

Keep source disjoint namespaces,using1088-byte evaluation records:
metadata buffer b=0,1,index i<2^127:bank b*2^93+(i>>34),offset low34(i)*64;
payload:bank2^94+(i>>32),offset low32(i)*256;
table b=0,1,j<2^92:bank2^96+b*2^63+(j>>29),offset low29(j)*1088;
count/next b=0,1,d<2^80:bank2^96+2^64+b*2^44+(d>>36),offset low36(d)*16.
All bank ranges<2^97 are disjoint and each record fits the window;1088*2^29<2^40. Computing1088*j as(j<<10)+(j<<6) is bounded ordinary word work. Every used byte is initialized by a paid write;unused capacity supplies no data/advice.

Generation preserves T0's21 evaluations and scalar in1040 shared bytes:130 paired load/store copy instructions;zero padding is never read. T1 remains in the committed window throughout all21 scalar products,so its field loads are already included in905/341 counts. Complete scalar extraction/Fr addition before leaving its bank. Output a0,b1 and beta then reside in shared scratch before payload/metadata switches. No stale pointer,bulk DMA or hidden second operand copy exists.

Radix:zero all2^80 two-limb counts,scan metadata counting selected80-bit digit,prefix ALL digits into separate next array,stable scan/scatter ENTIRE64-byte records;swap buffers,low80 then high80. Digits use64+16 bits at0/8 or10/18,unaligned fixed-word load permitted.128-bit ordinals/counters always two limbs. Payload never sorted. Count/scatter endpoint96+176 per record/pass=544L;scan128L. Full24-word S and4-word beta comparisons and all losing paths charged.

Generation key2304+projection11776+finger640+PAIR960+Fr-add256+main704=16640. Main704 includes two table selectors64,T0 evaluation/scalar copies130,projection arguments48,PAIR encoding32,fingerprint packing48,payload selector/full256-byte copy64,metadata construction/selector/full64-byte copy64,two-limb indices/control80,frame/state64,remaining bookkeeping10;sum704. Inclusive kernel counts pay direct T1 operand loads. No extra network is present in B. Radix544,scan128,confirmation32 bring work17344L;retain17568 with224 spare,then17600 as complete hot bound. Every selector is paid at32 and all used bytes/copies are accounted.

Each confirmation<=4096:512 includes four32-instruction selectors,metadata/ordinal copies,fullS/beta comparisons,counters/frames/control. An EXTRA256 covers full endpoint/padding copies,still768<4096. All false fingerprints/repeats and losing paths paid. With cap floor(L/128),confirmations<=32L already in the base. FIRST useful pair triggers only11 cold terminal attempts.

Fixed<=2^92 instructions pays input decoding/conversion/subgroups,h pairing,C0/eleven phases,seed/FING including every capped rejection,all code/constants/shared scratch/controller bindings,ALL count zeros/prefix scans and ALL11 scalar/G2/pairing/output attempts. B uses a fixed bank decode/latch mechanism,not per-leaf processors or software page tables. Radix loops over all2^80 digits cost<2^90 including selectors/two-limb operations. Fixed exponents<4608 bits and conventional bounded field/curve routines fit this loose cap. All application bytes are separately written and paid;unused RAM supplies no data/advice.

Each2^93 target table entry costs<=512 dense references*20000+1024 full SCAL attempts*8192+2^18 evaluation/weight/scalar/record/index/selector/control=18890752<2^25. In particular six finite points times3 degree3 polynomials require54 Horner field products and54 adds;weighted T1 requires21 more products. Even pricing all75 at the coarse2048 price,adds/conversion at128 and every record write/control leaves room in2^18. Public d inverses/native constants,seed/FING,full input/subgroup checks,h pairing,C0/eleven phases,code/controller bindings,ALL count zeros/prefix scans and ALL11 terminal scalar/G2/pairing/output attempts fit the separate fixed2^92 bound. Thus table2^118+fixed2^92<2^120,complete cold and finish included,no amortization/advice.

Allocated namespaces:two metadata buffers together2^134 bytes,payload2^135,table bank padding<=2^104,count/next2^85,shared/code<2^46,total<2^136;retain peak<2^137. N=17600L+2^120;EXACT integer N^25<20000^25*2^3139 proves submitted125.56. Full protected p/r,supplied P,threshold0.39 and C20000 remain unchanged. No production run,physical bank approval or qualification is claimed. The optional shared-Docker audit interprets actual905 instructions,all21 pinned cores and115 field helpers;compares with independent dense convolution at full p;checks ABI,stack/aliases/initialization and11631 steps. It does not certify H/B or execute production inversion. No pass is claimed;all candidate code is NEVER host-executed or run with credentials. Shared llm-verifier alone owns execution/review/scoring.
