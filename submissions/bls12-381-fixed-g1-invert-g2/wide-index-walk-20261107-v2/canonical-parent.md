# Direct field calls for exact trace canonicalization

canonical-parent.S supplies the same fapi_trace_canonical(raw960,outS192),exact projection/cubic fibers/four-image minimum and status/k contract as its predecessor. It replaces the descriptor interpreter and five negation wrappers with nine explicit direct BL calls. No field value,canonical key,orientation,hash message,raw trajectory,failure event or success law changes. Sources are static untrusted evidence,never compiled/executed on the host or with credentials.

## Fully initialized workspace and constants

The fixed shared816-byte aligned workspace is:TRAW0..191,ANEG192..239,Bvariants240..623,SUM624..671,DIFF672..719,native s720..767,zeroFp768..815. x21=base,x23=ANEG,x22=Bvariants,x24=zeroFp. The last48 bytes are an immutable canonical native field zero;every byte is explicitly initialized and paid cold. No descriptor/function-pointer table or interpreter is present. All code,address bindings,workspace and zero bytes are charged per target,not free program advice.

Compute actual C=(1+u)^((p-1)/2)=(s,s) in the pinned Fp2 tower cold and write canonical Montgomery s at the same exported fapi_trace_s720. This fixes the sign from the exact tower rather than an arbitrary-root oracle. The full protected p,Montgomery n0 and all native helper constants are initialized cold. Place code/data/helpers in the shared image below2^45,in direct branch/ADRP range,with no PLT,veneer,C wrapper or per-step constructor. All816 bytes are disjoint from main context,current raw state and NH seeds;the code overwrites every transient input before reading it.

## Complete direct call sequence and actual ABI

The nine calls,each with THREE explicit argument instructions and its BL,are:

1. sub(zeroFp,TRAW+48,ANEG) to negate Aimag.
2. add(TRAW+96,TRAW+144,SUM).
3. sub(TRAW+96,TRAW+144,DIFF).
4. native_fp_mul(output Bvariants+96,input SUM,input s).
5. native_fp_mul(output Bvariants+144,input DIFF,input s).
6. sub(zeroFp,TRAW+96,Bvariants+192).
7. sub(zeroFp,TRAW+144,Bvariants+240).
8. sub(zeroFp,Bvariants+96,Bvariants+288).
9. sub(zeroFp,Bvariants+144,Bvariants+336).

This is exactly the earlier two-product/seven-field-operation computation. Each negation is the exact field subtraction0-input modulo full p;its output remains canonical. The direct calls match the actual helper signatures:field add/sub(input,input,output),native multiplication(output,input,input). Every destination is disjoint from its sources or explicitly safe. SUM/DIFF precede their products;candidate1 precedes its negation;zeroFp and s are never written during the walk.

The field body cost is now2*375+36+6*37=1008. Five old negation wrappers of5 instructions each are not invoked or charged as if they ran;their argument work is instead explicit in the parent. The old field program's52 instructions and its3 parent setup/call instructions are replaced by9*4=36 direct argument/BL instructions. This is an implementation saving,not a changed reference price or unpriced interface. Pinned native multiplication still includes its29+341+5 body/constant adapter=375. Full field add/sub include modulus loads andRET at36/37. No SIMD or unbounded-width primitive is used.

## Projection,minimum and preservation

The parent saves x19..x24 and originalLR in its64-byte aligned frame. x19=raw,x20=outS,x21=workspace,x22=Bvariants,x23=ANEG,x24=zeroFp. All helpers preserve these callee registers,trace-min uses caller scratch,x25..x28 are untouched,and x18 is untouched. Every called leaf and final parent restore/RET are included. There is no hidden context/frame work or overwritten originalLR.

Copy native a0 at0..95 and b1 at384..479 into TRAW with12 LDP/12 STP=24. Test all12 B limbs in18 includingCBZ;B0 iff identity on GT returns STATUS0 and skips all later calls. Copy raw B96 to candidate0 in12 memory instructions. Then execute the nine direct calls and invoke trace-min(TRAW,ANEG,Bvariants,outS,raw+944). The minimum writes every192 output byte and k0..3;the parent explicitly returns STATUS1 even when k0. Scalar-A metadata576..607 is never accessed. Identity paths are cheaper and their uninitialized output is never consumed because the main aborts.

Every paired displacement is aligned and<=464;all workspace ADD offsets<=768 and the k pointerADD944 are encodable. Native pointer width8;no pointer/address oracle or atomic wide store is used. The .zero declarations describe publicly initialized data;cold construction pays every byte and computes actual s. No stale B candidate or scratch is read:all four B candidates,ANEG,SUM andDIFF are written before the minimum,while zero and s are cold-initialized constants.

## Complete worst-path computation

Direct parent source has115 instructions conservatively counting BOTH mutually exclusive status assignments. This includes its full frame,address binding,projection,B reduction,initial copy,all36 direct field-call argument/BL instructions,the minimum's6 argument/BL instructions,status branch,restore andRET. No separate interpreter or negation wrapper is executed. The longest minimum path is131 including every comparison/output/k store/RET;the other path is88. The field bodies total1008 as above.

Thus canonical cost<=115+1008+131=1254,retained1264 with10 spare instructions. This is deterministic for every nonidentity GT input;identity exits are cheaper. No orientation-average/cap failure is introduced. The initial canonical call and new768+48-byte image/constant setup are paid in the complete per-target cold bound. Old descriptor/interpreter/negation files are omitted because this implementation does not call them.

Whole ordinary scalar work<=1264canonical+460NH+776fullBLAKE3+26Bupdate+72controller+160dictionaryfixed+12cappedvisits=2770,retained2800 with30 extra spare instructions. One protected reference C20000 per transition plus coldD=540*2^119 gives work<=(22800K+540*2^119)/20000,K=ceil9943R/10000. Diagnostic log2=125.81719586241536;exact certificate(22800K+540*2^119)^2500<20000^2500*2^314543 supports reduced proposed314543/2500=125.8172. Full finite success remains>0.39001077987607>0.39,including all singular/fingerprint/dictionary/sampler/identity/good-table losses. All protected p/r,exact supplied P,six tracks,one metric,threshold0.39,reference prices,zero advice and capacity<2^136 remain unchanged. H and explicit finite banked-memory admissibility remain material;this is conditional analytic evidence,not an official score,qualification,promotion,practical execution or fixed-hash theorem.
