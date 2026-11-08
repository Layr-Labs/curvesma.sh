# One-fold Mersenne Horner and complete EXTR extraction

This computes exactly the same five evaluations U as universal-proof.md. Let
b=2^61 and q=b-1. All arithmetic uses ordinary scalar64 instructions. The source
5039bd36's thirteen-instruction update folded twice. One fold already yields a
value below2q, permitting an exact single subtraction/select. No heuristic or
field-primality assumption beyond the reproduced q certificate is added.

For canonical h,k<q and digit c<2^60, compute P=h*k+c using a full128-bit
MUL/UMULH pair and carry. The following is exactly ten scalar instructions:

    MUL lo,h,k
    UMULH hi,h,k
    ADDS lo,lo,c
    ADC hi,hi,xzr
    LSR tmp,lo,#61
    ORR tmp,tmp,hi,LSL #3
    AND lo,lo,q
    ADD lo,lo,tmp
    SUBS tmp,lo,q
    CSEL h,tmp,lo,HS

P<2^122, so hi<<3 does not overflow64 bits and its low3 bits are zero;
lo>>61 fits3 bits. Their OR is exactly floor(P/2^61), not an approximate
high half. The other operand is P mod2^61. Thus their sum T is congruent to
P modulo q. For a precise bound, h,k<=b-2 and c<=b/2-1 give
P<=b^2-(7/2)b+3. Therefore floor(P/b)<=b-4, P modb<=b-1, and
0<=T<=2b-5=2q-3<2q<2^62. No64-bit addition overflow occurs. SUBS/CSEL
therefore returns T-q when T>=q and T otherwise, exactly the canonical
P modq. q itself reduces to zero. Every subsequent Horner input is canonical.
The 128-bit product addition cannot overflow128 bits because P<2^122.

Use the same five keys and39 digits c38,...,c0. Rather than execute the first
Horner step from zero, initialize all five accumulators directly to c38. This
is valid because c38<2^60<q. Then perform38 updates for each of five keys:
190 full ten-instruction updates,1900 instructions. The five initialization
MOVs replace the source's five zero MOVs and are in the initialization allowance.
This computes the identical polynomial at every key for every possible input.

The 37 initialized little-indexed input limbs are d0,...,d36, as completely
constructed and charged in preparation-audit.md. For each digit l=0,...,38,
put j=floor(60l/64) and r_l=60l mod64. Fully unroll:

    LDR lo,[input,#8*j]
    LDR hi,[input,#8*(j+1)]
    EXTR digit,hi,lo,#r_l
    AND digit,digit,#0x0fffffffffffffff

EXTR concatenates (hi:lo), shifts right by the immediate in0..63, and extracts
its low64 bits. Hence the result after AND is exactly (M>>(60l)) & (2^60-1).
At r_l=0 the instruction returns lo; loading hi remains charged. The largest
j is35 at l=38, so j+1=36 is initialized and no out-of-bounds read occurs.
All byte offsets are fixed and below296, using ordinary immediate loads.
The60-bit mask is an encodable logical immediate. Reuse lo/hi before their
values become dead at the next arithmetic update. There is no runtime division,
loop counter, table of offsets or free arbitrary-width extraction. Total39*4=156.
The general EXTR semantics are documented by Arm's Instruction Set Reference:
https://documentation-service.arm.com/static/6245c734b059dc5ff9a8bdab, D2.83.
The semantics needed for this exact schedule are reproduced above.

Five accumulators, five keys, q, digit, lo, hi, tmp, input/output bases and one
saved pointer fit within19 registers, as in the inherited leaf. Use x0..x17
and one saved x19 when x18 is ABI-reserved; the routine is a leaf, preserves LR
and uses no runtime calls. Allow128 instructions for loading five keys, five
initialization MOVs, q construction, input/output setup, saved-register frame,
five output REV/ST pairs, argument preservation, stack control and RET. A
concrete budget is: five key LDRs5, initialization MOVs5, q at most4, pointer
moves/address setup16, frame save/restore and stack allocation12, output REV/ST10,
RET1, and other control/moves32:85<128. Any fixed sign/digit input preparation
is already charged once in preparation-audit.md. No data-dependent loop occurs.

Complete compressor cost is1900+156+128=2184. Charge2208, with24 spare. The
separate512 preparation allowance remains unchanged, not silently reduced.
Other work is3296 SHA+128 sign+4992 outer+384 reference interface. Total other
scalar work is11008 per transition; the dense reference product is H=1.
No candidate code or assembly is host-executed. This is a static complete
instruction-equivalent schedule, not a measured production timing.
