# Complete exact-trace canonical parent

canonical-parent.S supplies fapi_trace_canonical(raw960,outS192),which the ordinary main calls. It instantiates the same exact Fp4trace/cubic-fiber/four-image minimum from fp4-trace.md with the same nine field computations and same minimum leaf. It changes no group operation,canonical key,orientation k,hash input,randomness,trajectory or success law. All source remains untrusted static evidence,never compiled/executed on the host or with credentials.

## Workspace and cold initialization

One16-byte-aligned shared1056-byte workspace holds:

|offset|bytes|value|
|---|---:|---|
|0|192|TRAW=(native a0,native b1)|
|192|48|negative Aimag|
|240|384|four96-byte B candidates|
|624|48|Bre+Bim|
|672|48|Bre-Bim|
|720|48|native Montgomery s|
|768|288|nine32-byte descriptors|

The whole program/code/data image,zero bytes,function/operand addresses and constants are initialized and paid cold. The .zero directives reserve/publicly describe bytes;they do not supply free runtime initialization or advice. Cold compute C=(1+u)^((p-1)/2)=(s,s) in the pinned Fp2 tower and write its actual s as a canonical native Montgomery coefficient at exported fapi_trace_s. This fixes the sign from the exact tower,not an arbitrary square root. The same pmod24=19 identity2s^2=-1 holds. No secret or searched target advice is provided. All workspace/code/control initialization is within the existing complete coldD;new image size/copy work<2^20 is negligible against its explicit slack.

Place workspace/code/helper objects in the fixed shared image below2^45,disjoint from main context,current raw buffers,NH seeds and each other. Direct static branches and ADRP data references are placed in range,with no C wrapper,PLT,veneer,unknown address translation or per-step constructor. Every table pointer is an ordinary64-bit initialized address. The fixed native reference's narrowed Fp coefficients and all provided leaves use canonical Montgomery residues;the key order compares those full six-limb values,not a truncated representation.

## Exact actual helper ABIs

Every descriptor is(function,arg0,arg1,arg2),matching the actual leaf:

1. neg(input TRAW+48,output ANEG);third argument ignored.
2. add(TRAW+96,TRAW+144,SUM).
3. sub(TRAW+96,TRAW+144,DIFF).
4. native_fp_mul(output Bcandidate+96,input SUM,input s).
5. native_fp_mul(output Bcandidate+144,input DIFF,input s).
6. neg(TRAW+96,Bcandidate+192).
7. neg(TRAW+144,Bcandidate+240).
8. neg(Bcandidate+96,Bcandidate+288).
9. neg(Bcandidate+144,Bcandidate+336).

These are exactly two native Fp products,one addition,one subtraction and five negations:2*375+36+37+5*42=1033 body instructions,as certified by the primary excerpts/constant adapter in scalar-audit.md. The zero-field operand for negation,full p,n0 and all code/constants are initialized cold. Descriptor order computes every source before use;outputs are disjoint from their inputs or explicitly safe. fapi_run_field_program has4 entry instructions,5*9 loop/dispatch instructions and3 exit instructions=52,including all descriptor loads,BLR,counters,frame and RET. Native helper bodies preserve their actual callee state and use exactly the ABI encoded here.

## Complete parent flow

The parent saves x19..x24 and originalLR in a64-byte aligned frame. x19=raw record,x20=outS,x21=workspace,x23=ANEG,x22=Bvariants,x24=program. Other callee registers x25..x28 are untouched. Field-program/native helpers preserve the six live values;the minimum leaf uses only caller scratch. x18 is untouched. There is no hidden ordinary wrapper/spill or return to an overwritten LR.

Copy raw native a0 at0..95 and b1 at384..479 into TRAW using12 LDP/12 STP=24 instructions. These are exact trace/3 coordinates by fp4-trace.md. OR all12 B limbs in18 instructions including its final CBZ;B0 is exactly identity on GT,returnSTATUS0 immediately. Otherwise copy B96 to candidate0 with6 LDP/6 STP=12. No stale workspace value is read.

Invoke the fixed nine-call field program with three argument/call instructions. It fully initializes ANEG,SUM,DIFF and candidates1/2/3. Then six argument/call instructions invoke fapi_trace_min(TRAW,ANEG,Bvariants,outS,raw+944). It initializes all192 output bytes and stores k0..3 at raw944. The parent sets STATUS1 after that call:it must not return k as status because k0 is valid. Identity skips every unneeded field/minimum call;no uninitialized output is subsequently consumed because the main aborts on STATUS0.

All paired offsets in the parent are aligned and<=464;the minimum's raw944 k pointer is formed by an ordinary ADD and single STR in the leaf. The parent frame's largest load offset48 and workspace address additions fit the actual AArch64 encodings. All scalar-A metadata576..607 remains unused padding;this parent reads only raw group bytes and writes k944.

## Deterministic complete instruction certificate

Direct parent source has82 instructions even if BOTH mutually exclusive success/identity status assignments are counted. That includes every frame/load/store,ADRP/ADD,projection,B identity reduction,initial B copy,program/minimum argument and BL,status branch,restore andRET. A success path actually uses81;retaining82 is conservative.

The existing minimum leaf has longest internal paths131 for zero Aimag,88 otherwise,including every comparison,output copy,k store andRET. The131 count can be checked statically:entry/zero-test15,zero-path B initialization1,three28-instruction comparisons84,branch1,complete output30. This is a worst-case bound for every input,no orientation-average cost or cap failure. Its call interface is already in the82 parent count.

Therefore complete canonical cost<=1033field bodies+52field-program+131minimum+82parent=1298. Retain1312,with14 additional spare instructions. This includes ALL interfaces and data movement;no old separate176/96/96 allowances are also added. The initial cold canonical evaluation uses the same paid code and is inside completeD. Identity paths are cheaper and do not weaken the hard bound.

Whole ordinary scalar work<=1312canonical+460NH+776fullBLAKE3+26Bupdate+72controller+160dictionaryfixed+12cappedvisits=2818,retained2848 with30 extra spare instructions. One protected reference C20000 per transition plus complete per-target coldD=540*2^119 gives work<=(22848K+540*2^119)/20000,K=ceil1243R/1250. Diagnostic log2=125.82037425140932;exact certificate(22848K+540*2^119)^1000<20000^1000*2^125821 supports proposed125821/1000=125.821. The full-size finite success bound remains>0.39005383358456>0.39. No reference price,parameter,threshold,six-track structure,one computation metric,zero-advice condition or memory<2^136 changes. H and explicit finite banked-memory admissibility remain material;this is conditional analytic evidence,not official scoring,qualification,promotion,practical execution or a fixed-hash theorem.
