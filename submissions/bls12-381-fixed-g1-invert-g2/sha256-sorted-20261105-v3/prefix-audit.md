# Exact SHA-256 messages and fixed prefix

This revision uses SHA-256, not SHA-512, for all four hash domains.
INIT,SCAL,WALK,UKEY are the respective16-byte ASCII strings defined in proof.md
and universal-proof.md. seed=SHA256(INIT||canonical_P||canonical_z) is32 bytes.
P is96 bytes,z576 bytes,so this INIT message is688 bytes and padded to11 blocks.
SCAL||seed||counter_be64 is56 bytes and pads to two64-byte blocks. Take its
low255 bits for bounded rejection sampling modulo r, at most1024 attempts per
scalar. UKEY||seed||counter_be64 is also56 bytes and pads to two blocks; take
low61 bits and reject q for the independently keyed polynomial compressor.
The counter increases on every attempted SCAL/UKEY hash and is below2^36.
Their distinct16-byte domain tags make their input sets disjoint.

WALK||seed||zero16||U(X_bytes) has length104 bytes, with a fixed64-byte prefix
WALK||seed||zero16 and a40-byte state-dependent suffix. zero16 is exactly16
zero bytes. Precompute the compression of that fixed first block from the
standard SHA-256 IV once per cold attack. The second block contains the40-byte
suffix,0x80,15 zero bytes,then the eight-byte big-endian integer832, the bit
length of the complete104-byte message. It is one64-byte block, not a hash of
the suffix alone. Load the precomputed eight-word chaining state, compress this
one block, and take the low24 bits of the standard256-bit digest as the label.
These are ordinary SHA-256 evaluations with exact standard padding, not a new
custom primitive or a free random oracle. All four domains are fully computed.

Under the sole explicitly unproved SHA-256 random-function idealization their
disjoint full message inputs yield the same independent uniform labels/scalars
as in the finite idealized proof. No claim transfers a proved fixed-SHA-512
property to SHA-256; neither fixed construction has that theorem. The change
is an alternative disclosed cryptographic heuristic of the same form.

sha256-audit.md gives the complete scalar implementation cost:6144 instructions
for the one variable block and formatting, including aligned compressor-input
preparation. Charge4096 separately for the unchanged polynomial compressor.
Cold INIT and prefix computations and all two-block SCAL/UKEY hashes are paid
inside2^90 fixed work. Even2^36 attempted samples at65536 instructions each
cost below2^52;65536 covers both blocks and their preparation. Pairings remain
below2^79 together,so the total fixed-work budget retains ample slack.
