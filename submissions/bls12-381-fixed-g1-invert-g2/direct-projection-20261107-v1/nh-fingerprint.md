# Small-seed NH fingerprint with proved compression

This replaces simple tabulation, not the exact Fp4 trace key or raw transition algorithm. Canonical S is still192 native bytes/24 unsigned64-bit words. The output U is still40 bytes/five big-endian64-bit words, consumed by the same full standard BLAKE3 JUMP/retry functions. The digest distribution is coupled to the same ideal one-sided independent-increment chain until a fingerprint collision. Actual fingerprints and trajectories change; no same-seed trajectory equivalence to the tabulation predecessor is claimed.

## Published universal-hashing ingredient

Material source: J.Black, S.Halevi, H.Krawczyk, T.Krovetz and P.Rogaway, UMAC: Fast and Secure Message Authentication, full version September8,1999, Theorem4.2, https://www.cs.ucdavis.edu/~rogaway/umac/umac_full.pdf . For even-length w-bit word messages, independently add uniform w-bit keys modulo2^w to each word, multiply adjacent pairs to2w bits, and sum modulo2^(2w). NH's collision probability for two different messages of the same length is at most2^-w. Here use24 words and w=64, giving128-bit NH outputs with pairwise collision probability<=2^-64. The paper supplies this theorem and construction; the following compression, five-lane accounting and scalar AArch64 implementation are derived explicitly here. No UMAC performance, MAC security, SIMD assumption or qualification is inherited.

## Independent odd-multiplier compression lemma

For independent uniform odd A in Z/(2^128), define g_A(x)=high64(A*x mod2^128). For distinct128-bit x,y put d=x-y mod2^128=2^s*d0, d0 odd,0<=s<128. A*d is uniform on residues of valuation s: there are2^(127-s) such residues. Equal high halves of A*x and A*y require their cyclic difference to lie in the two intervals of residues with circular distance<2^64 from zero. If s>=64 this is impossible. If s<64, the two intervals together contain at most2^(64-s) residues of valuation s. Therefore Pr[g_A(x)=g_A(y)]<=2^-63. This is only a necessary-condition count, so no independence of A*x and A*d is presumed. It is a finite exact argument, not a heuristic truncation of NH.

For each of five lanes choose independent NH keys and an independent odd128-bit A. Let n_l(S) be its NH result and U_l=g_A_l(n_l). For two distinct S,S', first separate the event of equal NH results (<=2^-64) from unequal results followed by compression collision (<=2^-63). Each lane collides with probability<=3*2^-64. Lane key sets are independent, so all five collide with probability<=(3*2^-64)^5<2^-312. No independence of NH outputs across different message pairs is needed.

Extend an ideal raw walk with independent fresh labels, independent of these fingerprint keys, until its first full-orbit repeat. Its exact canonical S values are distinct until that event. There are K+1<2^126 values; a union bound over pairs and the preceding fixed-pair bound gives fingerprint coupling failure<2^251*2^-312=2^-61. Retain the conservative loss2^-55. This is the same ideal-path/dummy-extension argument as probability.md, not an invalid adaptive pairwise-universality union bound. On no fingerprint collision, all JUMP messages for distinct orbit keys are distinct, so the disclosed standard BLAKE3 bridge H supplies fresh independent labels as before. Conditional phase-bin bucket laws and dictionary-cap arguments are unchanged.

## Exact target-specific seed construction

Let NHK be ASCII FAPI-NH-K-v1 padded by four zero bytes to16; NHC is FAPI-NH-C-v1 with the same padding. These new domains are distinct from INIT,SCAL,JUMP and all initial/retry domains. For pair i0..11,lane l0..4, generate B3(NHK||seed||BE64(5i+l)||zero8), taking the first16 digest bytes as two little-endian64-bit keys. Each message is64 bytes, counter0,len64,flags11, standard unkeyed BLAKE3-256, all seven rounds. Store these60 pairs consecutively in shared seed-record bytes0..959. For lane l generate B3(NHC||seed||BE64(l)||zero8), take first16 bytes as little-endian128-bit A, force its low bit to1, and store at960+16l. Under H, all65 inputs are distinct and their outputs independently uniform; forcing the low bit samples exactly uniform odd multipliers. This uses no extra random oracle, hidden advice, table search or uncomputed seed.

Initialize every byte of the1040-byte record by paid cold hashes/writes, and store its native pointer at controller context48. The same outer fingerprint ABI (key pointer, seed pointer, out40) is retained. Seed derivation plus preparation/interfaces is below2^20 scalar instructions, replacing vastly larger previous tabulation construction. All new tag/code/constants bytes are initialized cold; complete coldD=540*2^119 remains conservative without advice or amortization. Shared capacity<2^45 remains more than sufficient. Every ordinary output byte is written; its initialized64-byte BLAKE3 block retains zero padding40..63.

## Scalar instruction certificate

The supplied tabulate.S is an unrolled leaf. It preserves x19..x22 in a32-byte aligned frame, leaves x18 untouched, and does not modify LR. x20 is the sequential seed cursor, x21 output and x22 key cursor. Accumulators are(x3,x4),(x5,x6),(x7,x8),(x9,x10),(x11,x12). Twelve input pairs are loaded by LDP with post-index16. Each pair uses its five independently keyed NH_LANE calls in order. ADD performs the required64-bit modular message/key additions. MUL/UMULH produce the full128-bit unsigned product; ADDS/ADC accumulate modulo2^128, including carry. No128-bit unit multiply or SIMD instruction is assumed.

For the first pair, write each product directly to its accumulator, with no uninitialized read or zeroing fiction. Each first lane costs5 instructions (seed LDP,twoADD,MUL,UMULH), plus one input LDP for the pair:26 total. Each of the remaining11 pairs costs1+5*7=36, total396. Callee preservation and three pointer copies cost5 on entry; restoring and RET cost3.

At the end, for NH accumulator lo+2^64*hi and multiplier a0+2^64*a1, the high half of the product modulo2^128 is

 high64(lo*a0)+hi*a0+lo*a1 mod2^64.

Each lane uses seed LDP, UMULH, two scalar MADD, REV, STR:6 instructions,30 across five lanes. MADD's64-bit result is exactly the required modular addition; a high128-bit multiplication primitive is not used. REV/STR serializes each result big-endian. Seed cursor consumes exactly960+80=1040 bytes; input exactly192; output exactly40. All offsets/post-indexes are aligned and encodable, caller scratch and preserved registers are disjoint, and all five accumulator pairs remain live until used.

Total complete fingerprint cost=5+26+11*36+5*6+3=460 instructions including memory, pointer work, frame, serialization and RET. Source is static untrusted evidence and has never been compiled/executed on the host or with credentials. Public checks of arithmetic or text are not attack experiments; shared llm-verifier owns isolated execution/review/scoring.

## Complete finite success and work

The fingerprint loss is conservatively2^-55. probability.md proves the other current losses:bad table2^-66,scalar cap2^-970,identity2^-100,dictionary cap2^-100,singular-match2^-180. Raw jumps h*z^v_j and initial h*z^beta0 are independently uniform for every nonidentity z;the finite repeat proof is unchanged. Success after all losses>0.39001077987607>0.39. The first history exponent is t+1;there is no uniform-denominator or count/norm-tail assumption.

The ordinary scalar subtotal is1232canonical+460fingerprint+776fullBLAKE3+26Bupdate+72controller+160dictionaryfixed+12cappedvisits=2738. Retain2768, with30 additional spare instructions. Complete work, including exactly one protected C20000 reference multiplication per accepted transition and complete per-target coldD, is <=(22768K+540*2^119)/20000. Diagnostic log2=125.81517009049527. Exact integer certificate(22768K+540*2^119)^1250<20000^1250*2^157269 establishes proposed157269/1250=125.8152. Peak physical capacity<2^136, zero advice, six production tracks and one computation metric are preserved. Full standard BLAKE3 H and explicit finite banked-memory admissibility remain material. This is conditional analytic evidence, not a measured score, qualification, promotion, practical production execution or fixed-hash theorem.

Material local sourcee17a155a-6758-4650-bd1f-47e087a733ef supplies the prepared full-round hash, exact trace key and finite probability/ledger framework. It is acknowledged; formal linkage is deferred while unpublished. NH's published theorem is independently credited above, not claimed novel or supplied by that predecessor.
