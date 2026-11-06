# Exact integer dot products and final Mersenne reduction

The mathematical compressor is unchanged: q=2^61-1,39 digits c_ell<2^60, f_M(t)=sum from ell=0 to38 of c_ell*t^ell mod q, evaluated at five independently initialized keys. This implementation produces exactly the same40-byte U(M) for every289-byte input and every permitted key, including zero keys and leading zero digits. It changes no hash label, group transition or probability argument.

## Cold powers and exact sum bound

After sampling each key k_i, construct t_(i,ell)=k_i^ell mod q for ell=1,...,38; the constant term c_0 is used directly. Store these190 canonical residues as five keys per nonconstant digit in a fixed1520-byte table. Each next power can use the predecessor's bounded61-bit product/reduction routine; at most190 such operations, writes and address work are far below2^16 scalar instructions. This is charged inside the existing2^90 cold fixed-work budget for every target. The key-dependent table is internal working storage, not external advice; include its1520 bytes in the unchanged2^50 peak bound.

Initialize each low sum word to c_0 and each high word to zero, then add the38 nonconstant terms. The final unreduced integer is still S_i=sum_(ell=0)^38 c_ell*t_(i,ell), where t_(i,0)=1 mathematically but is neither stored nor multiplied. Every term is below2^60*q<2^121. Thus S_i<39*2^121<2^127. The entire sum fits in two ordinary64-bit words (lo,hi). It is not an arbitrary-width unit-cost integer. Each term uses:

LDR power,[table,#fixed_offset]
MUL prod_lo,c,power
UMULH prod_hi,c,power
ADDS sum_lo,sum_lo,prod_lo
ADC sum_hi,sum_hi,prod_hi

All additions are the ordinary two-word carry schedule. No high carry is lost since every partial sum is nonnegative and at most the final sum below2^127. Charge five scalar instructions per nonconstant term,190 terms total950. The five initial c_0 assignments replace five zero assignments in the same fixed initialization budget. All table offsets are fixed nonnegative multiples of8 below1520, within immediate LDR addressing. The same digit is extracted once and reused for five keys. Fully unroll all39 digits and five keys; no runtime loop counter, exponent computation or dynamic address work is concealed in these five instructions.

## Exact final reduction in twelve instructions

For a two-word S=lo+2^64*hi<2^127, define a0=lo & (2^61-1), a1=(lo>>61)|((hi & (2^58-1))<<3), a2=hi>>58. These are precisely the three base2^61 digits, and S=a0+2^61*a1+2^122*a2. Hence S mod q=(a0+a1+a2) mod q. Here a0,a1<=q and a2<=31. Their sum T<=2q+31<2^63, so the scalar additions below cannot overflow. Fold once with V=(T&q)+(T>>61). Then V<=q+2<2q, so one conditional subtraction yields the unique canonical residue. This is exact for every bounded S, including S=0 and values near a word boundary.

Use this complete scalar schedule, with q preloaded and scratch registers reused:

AND a0,lo,q
LSR a1,lo,#61
UBFX tmp,hi,#0,#58
ORR a1,a1,tmp,LSL #3
LSR a2,hi,#58
ADD a0,a0,a1
ADD a0,a0,a2
LSR tmp,a0,#61
AND a0,a0,q
ADD a0,a0,tmp
SUBS tmp,a0,q
CSEL residue,tmp,a0,HS

Exactly12 instructions per key,60 total for five keys. UBFX and shifted ORR are scalar64-bit instructions; neither is a wider arithmetic primitive. The no-borrow condition implements residue=V-q if V>=q and V otherwise.

## Extraction, register allocation and all overhead

Retain the aligned37-word representation d_0,...,d_36 of the injective289-byte integer M. For digit ell in0,...,38 let j=floor(60*ell/64), s=(60*ell) mod64. Load lo=d_j and hi=d_(j+1), use EXTR digit,hi,lo,#s, then AND digit,digit,#(2^60-1). This is exactly the low60 bits of M>>(60*ell). All digits use four scalar instructions: two LDR, one EXTR, one logical-immediate AND. The mask of60 consecutive low one bits is an encodable logical immediate. All offsets are fixed at compile time. The maximum j is35, so j+1<=36 and every load is inside the fully initialized296-byte buffer, even at s=0. EXTR's immediate is always0..63, so no shift by64 or undefined boundary case arises. Fully unroll39 extractions, cost156. c_0 initializes the sums; remaining digits are shared across all five keys. Arm's scalar EXTR definition is in its Instruction Set Reference Guide, sectionD2.83, linked in prefix-audit.md.

Ten registers hold the five two-word sums. One holds the shared digit; two hold a loaded power/product halves (the product can overwrite the consumed power after last use); pointers for input, powers and output and one q register use four more. One extraction scratch gives at most18 registers. LDR power may use a separate scratch reused before MUL/UMULH, or schedule UMULH before overwriting a source. For clarity, use two product registers and one power register: the shared digit and input-extraction scratch do not need to coexist after extraction, so reuse that scratch for power. This still fits18. The reduction reuses consumed sum and product registers sequentially after all terms. No calls or per-term spills occur.

Charge128 more instructions for five low-word c_0 assignments and five high-word zero initializations, q and pointer setup, caller BL and return, an optional preserved-register pair on ABIs reserving x18, five big-endian output reversals/stores and any fixed interface. This comfortably covers the listed fixed work. Total950+60+156+128=1294; charge1344, leaving50 slack. Every power-table load and output store is included. Code is fully unrolled; its ordinary code storage remains below the broad peak memory allowance.

The separate SHA1024 preparation allowance includes forming d from the576-byte full group state's b and sign. Read/reverse36 b words and form d_0=(B_35<<8)|bit, d_j=(B_(35-j)<<8)|(B_(36-j)>>56) for1<=j<=35, and d_36=B_0>>56. Reading/reversing costs at most72,35 shift/shift/OR/store constructions140, the two boundary constructions at most8, and moving the retained neighbor at most36. Including sign load, fixed pointers, optional saved registers and caller/return, this is below512. It is computed once per transition in the checkpoint pipeline, plus one fixed initialization; comparisons reuse d.

The remaining SHA preparation is bounded by128: copy eight32-bit chaining words with LDR/STR pairs(16); copy five64-bit compressor-output words into the40-byte payload(10); construct the three remaining64-bit padded words including0x80 and832-bit length (below16, counting immediate construction, byte reversal and stores); extract the40-bit label from the32-byte digest (below8); and at most64 for all pointer/control/interface instructions. Sum below114, within128. SHA's own input parsing, round/schedule loads, feed-forward and block function interface are charged separately in its128 fixed-block allowance. All payload bytes and zero padding are initialized. Therefore whole preparation is below512+128=640; charge1024, leaving384 slack. No required serialization, buffer copy, padding or control has been discarded.

Standard SHA-256 costs64*24+48*14+128+1024=3360 including that preparation. The final step budget is1344 compressor +3360 SHA +128 sign +4640 outer +384 interface =9856 outside one named dense reference multiplication. Whole cold work is(29856*K+2^90)/20000. The exact score inequality remains in horizon-audit.md. No optional experiment or measured instruction trace is claimed.
