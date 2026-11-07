# Direct raw-A projection and fused B initialization

canonical-parent.S computes exactly the same full-size Fp4 trace key and orientation as its material predecessor wide-index-walk-20261107-v2. The native a0 is already contiguous in raw0..95. trace-min.S reads only these first96 bytes from its raw-A input, so pass the preserved raw pointer directly. No temporary A copy is needed. Copy raw b1 at384..479 directly to B candidate zero, combining the identity reduction with these six LDP/STP pairs. The eleven ORRs and CBZ cover all twelve B limbs. B=0 iff identity on GT; this path returns STATUS0 before using an output.

## Workspace and lifetimes

The shared624-byte aligned workspace is Bvariants0..383,ANEG384..431,SUM432..479,DIFF480..527,native s528..575,zeroFp576..623. x21=base,x22=base=Bvariants,x23=base+384,x24=base+576. Raw A is borrowed read-only through x19. The outputS192 in the main context is disjoint from raw and workspace. All624 workspace bytes and exact native s/zero constants are initialized and paid cold. Actual C=(1+u)^((p-1)/2)=(s,s) is computed in the pinned Fp2 tower and s written canonically Montgomery. There is no arbitrary-root, stale-data or free-initialization oracle. All code/data/helpers fit the same shared image below2^45 with direct BL/ADRP bindings, no veneers/PLT/wrappers.

The nine direct calls each have three explicit argument instructions and BL:
1. sub(zeroFp,raw+48,ANEG).
2. add(Bvariants,Bvariants+48,SUM).
3. sub(Bvariants,Bvariants+48,DIFF).
4. native_fp_mul(output Bvariants+96,input SUM,input s).
5. native_fp_mul(output Bvariants+144,input DIFF,input s).
6. sub(zeroFp,Bvariants,Bvariants+192).
7. sub(zeroFp,Bvariants+48,Bvariants+240).
8. sub(zeroFp,Bvariants+96,Bvariants+288).
9. sub(zeroFp,Bvariants+144,Bvariants+336).

Field add/sub use input/input/output and native multiplication output/input/input. All destinations are disjoint from their source operands. Candidate zero is initialized before every subsequent read; SUM/DIFF precede products; candidate one precedes its negations. Neither raw group bytes nor immutable zero/s are written. The minimum is called with(raw,ANEG,Bvariants,outS,raw+944), reads A only at raw0..95, writes all192 output bytes and k0..3, and does not overwrite raw A or any B candidate. The parent returns STATUS1 explicitly even for k0. Scalar-A metadata576..607 remains unused padding.

The parent preserves x19..x24 and LR in its64-byte aligned frame. Other callee registers and x18 are untouched. All helper frames and returns are counted separately. Every pointer offset is encodable; paired raw accesses384..464 are within576 group bytes, paired workspace accesses0..80 cover B0. No pointer or memory-search oracle is introduced. The smaller workspace is cold initialized, not recycled uninitialized storage.

## Complete count

The parent has85 instructions conservatively counting both mutually exclusive status assignments. Relative to its115-instruction predecessor, twelve instructions copying A disappear, twelve copying B to a second candidate disappear, and six redundant LDPs used only for the identity reduction disappear. Its projection/reduction is now12 copy instructions+11 ORRs+CBZ=24, replacing54. All36 field argument/BL instructions, six minimum argument/BL instructions, frame/address/status/restore/RET remain included.

The field bodies are unchanged:2*375+36+6*37=1008. The full minimum worst path is131 including comparisons, every output/k store andRET. Thus canonical cost<=85+1008+131=1224, retained1232 with8 spare. This is a deterministic bound for every nonidentity GT input. No orientation average, additional failure, hash change, group operation or changed walk is involved.

Ordinary1232+460+776+26+72+160+12=2738, retained2768 with30 spare. One protected reference C20000 per transition and full coldD540*2^119 give complete work<=(22768K+540*2^119)/20000, K=ceil9943R/10000. The exact integer certificate(22768K+540*2^119)^1250<20000^1250*2^157269 establishes reduced proposed157269/1250=125.8152, above diagnostic125.81517009049527. Full finite success remains>0.39001077987607>0.39, including all generator/NH/dictionary/sampler/identity/singular losses. W94 construction and addressing are fully paid in cold-construction.md. Protected parameters, supplied P, six tracks, one computation metric, success threshold, C20000, zero advice and memory<2^136 remain unchanged. Standard full BLAKE3 H and explicit finite banked-memory admissibility remain material. This is conditional analytic evidence, not an official score or qualification. Candidate evidence is static and never host-executed.
