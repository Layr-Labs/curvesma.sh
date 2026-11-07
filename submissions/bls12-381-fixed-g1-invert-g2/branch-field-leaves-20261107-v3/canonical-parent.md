# Native canonical negation and complete direct parent

This revision uses a dedicated fapi_fp_neg(input48,output48) in field-ops.S for five exact negations, replacing general sub(zero,input,output). This material implementation change preserves every full-size field value, orbit key, orientation, hash input, raw walk and success condition. Candidate sources remain static untrusted analytic evidence, never compiled/executed on the host or with credentials.

## Exact leaf and deterministic count

Inputs are canonical ordinary six-limb representatives a of native Montgomery values, 0<=a<p. Negation in this representation is0 when a=0 and p-a otherwise; it represents the same field negative because Montgomery scaling is linear. The leaf loads all six input limbs with3 LDP, combines them with5 ORR and tests with CBZ. The zero path stores all48 zero bytes with3 STP andRET:13 instructions including the common prefix. The nonzero path binds the actual protected p with ADRP/ADD, loads its six limbs with3 LDP, performs full p-a through SUBS/four SBCS/final SBC, stores all48 bytes with3 STP,andRET. Total worst path=3+5+1+2+3+6+3+1=24 instructions. No loop, average, reduced modulus, unspecified normalization or hidden adapter exists. All helpers and constant/address bytes are initialized/paid cold.

For nonzero canonical a,0<p-a<p, so no conditional modular correction is required. For zero, the branch writes exactly zero. All arithmetic is over six native64-bit limbs, with borrow propagated explicitly. The leaf uses only x0..x14, preserves x18 and all callee registers, and never calls another function or changes SP/LR. x1 remains the output pointer; all input limbs are loaded before any store. Inputs/destinations are disjoint in the supplied parent. There is no unbounded-width instruction or memory-search oracle. The branch does not create an extra scientific failure event; the deterministic bound covers both paths.

The five unary calls require two argument instructions plus BL each, rather than three plus BL. Parent field dispatch therefore costs31 instead of36. Field bodies are2*375+28+28+5*24=926. Full native multiplication still includes its375-instruction body/constant adapter; add/sub now use the complete28-instruction branch leaves described below, including modulus loads andRET. No protected reference price changes.

## Parent workspace, calls and raw lifetimes

The fully cold-initialized shared workspace is576 bytes: Bvariants0..383,ANEG384..431,SUM432..479,DIFF480..527,native s528..575. x21=base,x22=base,x23=ANEG. There is no zero operand buffer. Compute the actual exact C=(1+u)^((p-1)/2)=(s,s) in the pinned tower cold and write canonical native s. Constants/code/data/helper bindings are initialized and charged per target. Direct BL/ADRP bindings fit the same shared image below2^45 without PLT,veneer or C wrapper.

The nine direct calls are:
1. neg(raw+48,ANEG).
2. add(Bvariants,Bvariants+48,SUM).
3. sub(Bvariants,Bvariants+48,DIFF).
4. native_fp_mul(output Bvariants+96,input SUM,input s).
5. native_fp_mul(output Bvariants+144,input DIFF,input s).
6. neg(Bvariants,Bvariants+192).
7. neg(Bvariants+48,Bvariants+240).
8. neg(Bvariants+96,Bvariants+288).
9. neg(Bvariants+144,Bvariants+336).

The parent borrows raw A0..95 directly read-only, copies raw B384..479 directly into candidate zero, and fuses its all-twelve-limb identity OR test into those six LDP/STP pairs. Only identity exits STATUS0. Every B/ANEG/SUM/DIFF transient is written before use. All destinations are disjoint from operands, raw and outputS192. It invokes trace-min(raw,ANEG,Bvariants,outS,raw+944), whose A reads stop at95. The minimum writes all192 key bytes and k0..3;the parent returns STATUS1 explicitly even when k0. Scalar-A metadata576..607 remains unused. No raw group bytes are overwritten during canonicalization. k944 is disjoint from every operand.

The reduced48-byte aligned frame saves x19/x20,x21/x22,x23/LR in three STP and restores them in three LDP. x24 is no longer used or modified; x18 and other callee registers are untouched. The frame is valid across every helper call and original LR survives. No implicit wrapper or once-per-step constructor is omitted. All pointer offsets and paired accesses are encodable.

## Complete work and success

The parent source has77 instructions conservatively counting both mutually exclusive status assignments: predecessor85 minus five unary-argument instructions, two frame instructions and one obsolete zero-pointer binding. Its full projection/identity test is24, field interfaces31, minimum interface6;every frame/address/status/RET is counted. The minimum remains131 worst-case instructions including every comparison,output andRET. Canonical work<=77+926+131=1134,retained1144 with10 spare. This is deterministic for every nonidentity GT input and adds no orientation average, cost cap or failure allowance.

Ordinary scalar1144+460+776+26+72+160+12=2650,retained2680 with30 spare. One protected C20000 reference per transition and full D540*2^119 give complete work<=(22680K+540*2^119)/20000,K=ceil9943R/10000. The exact integer certificate(22680K+540*2^119)^625<20000^625*2^78631 supports reduced proposed78631/625=125.8096,above diagnostic125.80958450596707. Finite success remains>0.39001077987607>0.39. All W94 construction,head clearing,code/constants,scalar/group work,first-state work and verified raw-twist/G2 finish remain cold paid. No cold bound is reduced merely because workspace shrinks. Full protected parameters,supplied P,six tracks,one metric,threshold,reference prices,zero advice and memory<2^136 remain. H for actual full standard BLAKE3 and explicit finite banked-memory admissibility remain material. This is conditional analytic evidence,not official scoring,qualification,promotion or practical memory installation.

## Full canonical addition and subtraction leaves

field-ops.S gives both leaves with the same actual ABI x0=a,x1=b,x2=output. Each loads a into x3..x8 and b into x9..x14 with six LDP, preserving x2 directly and avoiding three pointer moves. Both use caller registers only and preserve x18/callee state/SP/LR. All input loads precede output stores. Constants and code are explicitly initialized/paid cold, no veneer or helper is omitted.

Addition: six carry-chain instructions compute a+b in x3..x8. Because full protected p<2^381 and a,b<p, the sum is below2^382, so no384-bit overflow is lost. Two address instructions bind p,three LDP load it,six SUBS/SBCS compute sum-p into x9..x14 while preserving sum. B.LO selects the original sum store block exactly when sum<p; otherwise it stores the reduced limbs. Each branch has three STP andRET. Both paths cost6+6+2+3+6+1+3+1=28. Every result is in[0,p). No six CSEL instructions or hidden data copies remain.

Subtraction: six SUBS/SBCS compute a-b modulo2^384 in x3..x8. B.HS skips correction exactly when a>=b; that path has17 instructions including stores/RET. On borrow,two address instructions and three LDP bind/load p,then six carry-chain instructions add p modulo2^384. Since -p<a-b<0, this gives the canonical value a-b+p in(0,p). Three STP andRET finish. Worst case6+6+1+2+3+6+3+1=28. The correction branch has no AND masks or conditional-data operations left to charge.

Relative to the material predecessor,addition36->28 and subtraction37->28 save17 actual instructions with no changed field result or failure allowance. The full canonical count77+926+131=1134 is retained1144,with10 spare; the resulting retained ordinary bound2680 has30 additional spare instructions. These deterministic branch bounds do not assume secret-independent timing or an average input distribution.
