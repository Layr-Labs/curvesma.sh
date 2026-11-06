# Exact integer dot products and final Mersenne reduction

The mathematical compressor is unchanged: q=2^61-1,39 digits c_ell<2^60, f_M(t)=sum from ell=0 to38 of c_ell*t^ell mod q, evaluated at five independently initialized keys. This implementation produces exactly the same40-byte U(M) for every289-byte input and every permitted key, including zero keys and leading zero digits. It changes no hash label, group transition or probability argument.

## Cold powers and exact sum bound

After sampling each key k_i, construct t_(i,ell)=k_i^ell mod q for ell=0,...,38, with t_(i,0)=1. Store these195 canonical residues as five groups per digit in a fixed1560-byte table. Each next power can use the predecessor's bounded61-bit product/reduction routine; at most195 such operations, writes and address work are far below2^16 scalar instructions. This is charged inside the existing2^90 cold fixed-work budget for every target. The key-dependent table is internal working storage, not external advice; include its1560 bytes in the unchanged2^50 peak bound.

For each i compute the unreduced nonnegative integer S_i=sum_(ell=0)^38 c_ell*t_(i,ell). Every term is below2^60*q<2^121. Thus S_i<39*2^121<2^127. The entire sum fits in two ordinary64-bit words (lo,hi). It is not an arbitrary-width unit-cost integer. Each term uses:

LDR power,[table,#fixed_offset]
MUL prod_lo,c,power
UMULH prod_hi,c,power
ADDS sum_lo,sum_lo,prod_lo
ADC sum_hi,sum_hi,prod_hi

All additions are the ordinary two-word carry schedule. No high carry is lost since every partial sum is nonnegative and at most the final sum below2^127. Charge five scalar instructions per term,195 terms total975. All table offsets are fixed nonnegative multiples of8 below1560, within immediate LDR addressing. The same digit is extracted once and reused for five keys. Fully unroll all39 digits and five keys; no runtime loop counter, exponent computation or dynamic address work is concealed in these five instructions.

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

Retain the aligned37-word representation d_0,...,d_36 of the injective289-byte integer M. Each60-bit digit reads at most two aligned64-bit limbs, shifts by fixed immediate amounts, combines and masks. Charge the full8 instructions per digit,312 for39 digits, as in the predecessor. All offsets are fixed below296 and no uninitialized padding is read. The final digit may be extracted more cheaply but the complete bound is retained.

Ten registers hold the five two-word sums. One holds the shared digit; two hold a loaded power/product halves (the product can overwrite the consumed power after last use); pointers for input, powers and output and one q register use four more. One extraction scratch gives at most18 registers. LDR power may use a separate scratch reused before MUL/UMULH, or schedule UMULH before overwriting a source. For clarity, use two product registers and one power register: the shared digit and input-extraction scratch do not need to coexist after extraction, so reuse that scratch for power. This still fits18. The reduction reuses consumed sum and product registers sequentially after all terms. No calls or per-term spills occur.

Charge128 more instructions for ten sum-zero initializations, q and pointer setup, caller BL and return, an optional preserved-register pair on ABIs reserving x18, five big-endian output reversals/stores and any fixed interface. This comfortably covers the listed fixed work. Total975+60+312+128=1475; charge1536, leaving61 slack. Every power-table load and output store is included. Code is fully unrolled; its ordinary code storage remains below the broad peak memory allowance.

The separate SHA2048 preparation allowance still includes forming d from the576-byte full group state's b and sign: read/reverse36 b words, form d_0=(B_35<<8)|bit, d_j=(B_(35-j)<<8)|(B_(36-j)>>56) for1<=j<=35, and d_36=B_0>>56. This is below512 scalar instructions. It is computed once per transition in the checkpoint pipeline, plus one fixed initialization; comparisons reuse the buffer. Copying the32-byte SHA chaining value and formatting the64-byte final block, padding,832-bit length and40-bit label extraction cost below512 more. The retained2048 leaves over1024 slack. No preparation is shifted into an uncharged phase.

Standard SHA-256 costs64*24+48*14+128+2048=4384 including that preparation. The final step budget is1536 compressor +4384 SHA +128 sign +4640 outer +384 interface =11072 outside one named dense reference multiplication. Whole cold work is(31072*K+2^90)/20000. The exact score inequality remains in horizon-audit.md. No optional experiment or measured instruction trace is claimed.
