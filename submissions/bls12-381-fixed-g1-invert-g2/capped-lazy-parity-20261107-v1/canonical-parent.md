# Capped lazy Frobenius candidates with a finite branch-count proof

This changes native implementation,not the raw walk,key,hash or group/scalar invariant. First compute canonical negative ANEG of raw Aimag and candidates B,-B. Compare raw Aimag with ANEG as unsigned canonical six-limb Montgomery values. If Aimag<ANEG,the unchanged trace-min.S necessarily selects even k0/2 and reads only B,-B;skip K(B),-K(B). Otherwise compute all four candidates before invoking the same minimum. Equality means Aimag=0 because p is odd;this path computes every candidate needed by the minimum's zero-imaginary branch. Odd parity also computes every needed candidate. No uninitialized B1/B3 is read on a skipped path. The output S192,k0..3 and explicit STATUS are identical for every nonidentity input.

canonical-parent.S contains all thirteen comparison/branch instructions and a ten-instruction heavy counter/cap check. The counter increments only on a heavy branch,before its expensive body. Its32-byte shared control stores count128 and cap128 H=floor((1/2+2^-20)K)+1. All bytes are cold initialized. The cold initial-state canonical call is charged separately;reset count to zero after that call and before the ordinary main loop. A cap overflow returns STATUS0 and fails without an unchecked result. All ADRP/ADD,loads,128-bit increment/store,128-bit comparison/branch are present. K<2^126 prevents counter overflow. The new control/code bindings fit the existing shared image and are cold paid.

## The zero-imaginary set is small at the full protected order

Let sigma(X)=X^p,U=0xd201000000010000,r=U^4-U^2+1,xi=p modr=-U. For raw X in GT,the native coefficient A=a0 satisfies6A=sum_(i=0..5)X^(xi^(2i)),by the exact Fp12/Fp2 trace in the pinned tower. Thus Aimag=0 implies A^p=A and

 sum_(even ell)X^(xi^ell) -sum_(odd ell)X^(xi^ell)=0.

Use the full exact representatives for ell0..11:
1,-U,U^2,-U^3,U^2-1,-U^3+U,-1,U,-U^2,U^3,1-U^2,U^3-U.
Multiplying this Laurent polynomial by X^(U^3) yields a nonzero polynomial of degree2U^3. Its12 exponents are distinct;the extreme coefficients are -1 in the full odd-characteristic field. Hence at most2U^3 elements of GT have Aimag=0. Exact public arithmetic gives2U^3/r<2^-61. No uniform-coordinate assumption or heuristic distribution is used.

Frobenius sigma is a bijection of GT and sends Aimag to its negative,which in canonical native Montgomery form is p-Aimag when nonzero. It bijects the sets Aimag<ANEG and Aimag>ANEG. Therefore,for a uniform raw group element,the heavy indicator I=1[Aimag>=ANEG] has mean mu=1/2+Pr[Aimag=0]/2<1/2+2^-62. This uses an exact bijection and polynomial root bound,not an assumed half-cost orientation average.

## Finite dependent-path concentration and explicit failure allowance

Condition on any good generator table from mixing.md. The ideal extended raw walk is stationary uniform and its nontrivial Fourier eigenvalues have modulus<=rho=2^-42. For centered indicator g=I-mu,Plancherel and the exact transition eigenvalues give

 |Cov(I_t,I_(t+d))|<=rho^d E[g^2]<=rho^d/4.

This is valid for the actual possibly nonreversible Cayley kernel:use the modulus of each complex eigenvalue,not an independent-endpoint assumption. Thus for S=sum_(t=1..K)I_t,

 Var(S)<=K/4*(1+2rho/(1-rho))<K.

The cap H above exceeds E[S] by more than2^-21 K. Chebyshev therefore bounds Pr[S>H]<2^42/K<2^-80,using the exact full K>2^124. Prefix counts are nondecreasing,so any ordinary cap overflow before termination is included in this extended-walk event. No post-repeat behavior of the actual deterministic hash walk is assumed:the ideal extension is only for bounding the pre-repeat stopping event under the same H/NH coupling. Add2^-80 to the existing failure union. Exact final success remains>0.39005147408482>0.39 at the two-term capK=ceil5153R/5000.

## Complete deterministic work on every execution

The full source contains100 parent instructions conservatively counting both mutually exclusive STATUS assignments. Of these,22 argument/BL instructions and10 counter/check instructions occur only on a heavy branch;the base parent has68. All frame/address/copy/identity/minimum argument/status/restore/RET instructions and the13-instruction parity comparison are included. Parent frames remain48 aligned bytes preserving usedx19..x23 and LR;x18,x24..x28 are untouched. There is no hidden wrapper,normalization,orientation update or dictionary/hash change.

Common field bodies are ANEG and two B negations,3*24=72. The heavy field body has two375-instruction products,two28 add/sub leaves and two24 K-negations,total854. The unchanged minimum is conservatively charged its full131-instruction worst path on EVERY call,even when nonzero Aimag would use88. Thus base<=68+72+131=271,heavy extra<=22+10+854=886. For ordinary executions including cap failures,total canonical work<=271K+886H+10. Exact arithmetic bounds this by728K+896;the896 rounding/failing-check allowance is within unchanged coldD540*2^119. No average-cost discount survives beyond the explicit capped count.

Retain728 ordinary canonical instructions. Whole scalar subtotal728+460NH+776fullBLAKE3+26Bupdate+72controller+160dictionaryfixed+12cappedvisits=2234,retained2264 with30 spare. One protectedC20000 reference per transition and complete coldD give work<=(22264K+540*2^119)/20000. The exact certificate(22264K+540*2^119)^10000<20000^10000*2^1258347 supports reduced proposed1258347/10000=125.8347,above diagnostic125.83460209360135. All per-target initialization,full W94 construction,control reset/rounding allowance,failures and verified raw-twist/G2 finish remain cold paid. Full parameters,supplied P,threshold,reference prices,six tracks,one metric,zero advice and capacity<2^136 are unchanged. Standard full BLAKE3 H and explicit finite banked-memory admissibility remain material;this is conditional analytic evidence,not official scoring,qualification,promotion or practical installation. Candidate sources are static untrusted evidence and never host-executed.

## Native field and layout contracts

Workspace remains576 initialized bytes:Bvariants0..383,ANEG384..431,SUM432..479,DIFF480..527,native s528..575. The separate32-byte heavy-control record is shared and initialized cold. Actual C=(1+u)^((p-1)/2)=(s,s) determines the native s through the full pinned Fp2 tower. Raw A is borrowed read-only through x19; B0 is copied from raw384..479 with its fused all-twelve-limb identity test. The outputS192 is disjoint. ANEG and B2 are always produced;B1/B3 are produced precisely before any minimum path that could read them. SUM/DIFF precede products;K(B) precedes its negations. All sources/destinations in the direct calls are disjoint. The parent preserves x19..x23/LR in48 aligned bytes;all helper frames andRET are counted separately.

field-ops.S neg(input,output) costs24 worst-case,13 for zero. Three LDP load six limbs,five ORR/CBZ test zero. The nonzero path binds p in two instructions,loads it in three,subtracts a with six borrow-chain instructions,stores in three andRET. Thus3+5+1+2+3+6+3+1=24;the zero path stores all48 zero bytes andRET,for13. For0<a<p,the output p-a is canonical;Montgomery scaling is linear,so this is the exact native field negative. Only caller registers are used,x18/callee state preserved. No hidden correction or adapter exists.

Addition/subtraction use(input,input,output) and retain x2 as output pointer directly. Each loads a,b with six LDP. Addition computes sum in six carry instructions;full p<2^381 ensures sum<2^382,no384-bit overflow. It binds/loads p in five,computes sum-p in six preserving sum,and uses B.LO to select one of two explicit three-STP/RET blocks. Each path6+6+5+6+1+3+1=28. Subtraction computes a-b in six borrow instructions and B.HS skips correction. The no-borrow path costs17 including stores/RET;the borrow path binds/loads p in five and adds p modulo2^384 in six,for worst28. Since -p<a-b<0 in that case,the result a-b+p is canonical. All loads precede stores;caller scratch only,no hidden masks/wrappers. Native multiplication still uses output/input/input and costs375 including its constant-binding adapter and pinned primary body. All code/constants/address bindings are initialized/paid cold.
