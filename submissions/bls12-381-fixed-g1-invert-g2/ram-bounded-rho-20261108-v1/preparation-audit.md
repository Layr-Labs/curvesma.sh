# Complete 512-instruction preparation schedule

This replaces the predecessor's 2048 preparation allowance. No SHA round,
expansion, compressor arithmetic, sign computation or outer-control work is
removed. Those retain their separate allowances. The math message and output
are unchanged except that only its low 25 bits select the smaller jump array.

Let B_i, i=0,...,35, be the numerical unsigned words in the big-endian
serialization of the six canonical Montgomery Fp coefficients of b, coefficient
order as in compression-proof.md and highest limb first inside each coefficient.
These are already stored in the full raw GT state, with each Fp using six native
little-endian limbs. B_i is therefore available by ONE scalar LDR at the fixed
raw offset corresponding to coefficient floor(i/6), limb 5-(i mod 6). No REV
is needed to obtain its numerical value: a native LDR of the selected limb is
exactly the number represented by that big-endian serialized word. No b_bytes
buffer is created. All offsets are public and encodable ordinary immediates.

The integer M represented by b_bytes||sign has limbs
  d0=(B35<<8)|sign,
  dj=(B(35-j)<<8)|(B(36-j)>>56), 1<=j<=35,
  d36=B0>>56.
Here sign is the one-bit result from the separately charged 128-instruction
sign schedule. Load all 36 B words in reverse i order, one at a time. Keep the
previous and current word in two registers. For d0 execute LSL, ORR and STR;
for every middle word execute LSL, ORR with shifted-register LSR #56, and STR;
for d36 execute LSR and STR. Thus the complete transfer/arithmetic is
36 loads+36*3+2=146 scalar instructions. One LDRB obtains the initialized sign
if not already in a register. Unroll all operations, so no runtime loop counter,
index table or dynamic coefficient permutation exists. Each of all 37 destination
limbs is written before the compressor reads it. Allow 32 further instructions
for pointers, sign load, arguments and boundary control: total at most 178.
The compressor's own 3072 budget pays all its digit loads and Horner work.

After U emits its five big-endian 8-byte words to a disjoint 40-byte buffer,
copy those bytes to the 128-byte final SHA block using five LDR/STR pairs:
10 scalar instructions. The compressor's output REV/ST stores are already in
its 3072, and SHA's internal parsing REV/ST is already in its fixed 128.
Set final block word 5's memory bytes to 0x80 followed by seven zeros: MOV
of 0x80 and STR of that little-endian machine value, two instructions. Set
words 6,...,14 to zero with nine STR xzr instructions. Set word 15 to the
big-endian byte encoding of 1344 (the original 168-byte message length times
8): MOV of 1344, REV and STR, three instructions. Word 14 supplies the high
zero half of the required 128-bit length. All 128 final-block bytes are written.
This padding uses no uninitialized bytes and is exactly standard SHA-512.

Copy the immutable 64-byte instance prefix chaining value into a separate
mutable eight-word state with eight LDR/STR pairs: 16 instructions. The full
80-round/64-expansion SHA routine works on that state and returns all eight
feed-forward words. To obtain the low 25 digest bits, load its last native
state word and AND #0x1ffffff: two instructions. SHA digest serialization is
not required, because the digest's low bits are precisely that numerical word's
low bits. The mask is an encodable A64 logical immediate, not arbitrary-width
unit-cost arithmetic. The prepared block and copied state cost 10+2+9+3+16+2=42.

Add a separate allowance of 128 scalar instructions for the complete preparation
frame, call/return boundaries, preserved live input/state pointers, stack
alignment, padding/output addresses and any saved-register pairs. Between calls
all preparation temporaries are dead. The compressor and SHA internal ABI work
remain charged in their own budgets; even counting those interfaces again here
leaves the complete preparation at most 178+42+128=348. Charge 512, leaving
164 instructions of slack. The two preparation phases reuse scratch after their
inputs are dead, but no memory-aliasing saving is required. All working buffers
fit within the fixed 1 MiB scratch allowance.

The final SHA cost is 80*22+64*14+128+512=3296. Other per-step work is
3072 compressor +128 sign +4992 full-checkpoint/scalar/control +384 reference
interface, for total 11872 nonreference scalar instructions. One actual reference
dense Fp12 multiplication is counted separately as H=1. No cryptographic or
SIMD instructions, free oracle or host-executed candidate code are used.
