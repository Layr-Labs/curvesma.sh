# Machine model and instruction ledger

All work is charged in scalar AArch64 instruction equivalents R and divided by C = 20,000
(`tracks/docs/COST_MODEL.md`, BN254 section). No BN254 operation is charged as one unit;
the dense Fp12 multiplication and the cyclotomic squaring are priced through their
implementation, down to the base-field kernels listed below.

## 1. Representation

Fp elements are four 64-bit limbs in Montgomery form x~ = x * 2^256 mod p, fully reduced
(0 <= x~ < p). The same kernels serve Z/r (exponent arithmetic) with the constants of r.
p and r are both below 2^254, so a + b < 2^255 never carries out of 256 bits, and the
Montgomery product of reduced inputs needs one final conditional subtraction.

Tower (the one the profile maps to, checked in `experiment.py` part A):
Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 9 + u. Grouped for the product as
Fp6 = Fp2[v]/(v^3 - xi), v = w^2, Fp12 = Fp6[w]/(w^2 - v). The profile basis
(w^12 = 18 w^6 - 82, u = w^6 - 9) maps to tower coefficients by
c_m(tower) = (c_m + 9 c_{m+6}) + c_{m+6} u; this conversion is done once for z at setup.

## 2. Leaf kernels

The three leaves are straight-line listings that read operands through x1, x2, write the
result through x0, use only caller-saved registers x0..x17 (no spills, no frame) and end
in `ret`. Their constants (m0..m3 and n0 = -m^-1 mod 2^64) sit at label K of a block per
modulus. `experiment.py` part K executes every listing on an instruction-level simulator of
the AArch64 subset used (ldr/ldp/stp/adr, mul/umulh, adds/adcs/adc, subs/sbcs, csel,
csetm, and), on 12 x 12 edge pairs plus random and all-limbs-large inputs for both p and
r, and checks the result against Python integers. Declared run: 519 cases per modulus,
0 failures. The executed instruction count is the listing length (straight line, `ret`
included).

| Leaf | Listing length | Call site (BL + up to 3 pointer set-ups) | Charged per call |
| --- | ---: | ---: | ---: |
| fp_mul (Montgomery product, product scanning then 4 word REDC rounds, 1 conditional subtraction) | 158 | 4 | **162** |
| fp_add (add, trial subtract, 4 csel) | 22 | 4 | **26** |
| fp_sub (subtract, borrow mask, masked add of m) | 23 | 4 | **27** |

Every non-leaf routine call (Fp2 multiplication, Fp2 squaring, Fp4 squaring, Fp6
multiplication, the Fp12 routines) is charged **32** instructions: its call site, a frame
with four callee-saved registers, pointer bookkeeping for its temporaries on the stack, the
epilogue and slack. Fp2 additions, subtractions and multiplications by xi are inline pairs
of leaf calls and are charged as those leaf calls.

Multiplication by xi = 9 + u: (a0 + a1 u)(9 + u) = (9 a0 - a1) + (9 a1 + a0) u, with
9t = 8t + t by three doublings and one addition: 8 fp_add + 1 fp_add + 1 fp_sub per Fp2
element (10 leaf calls). No multiplication is used.

## 3. Counted Fp12 routines

`experiment.py` implements the two routines in a counting tower (every leaf call and every
routine call increments a counter) and checks their results against an independent
schoolbook product on random GT elements of the real BN254 group (part A, `mul_ok`,
`sq_ok`).

| Routine | routine calls | fp_mul | fp_add | fp_sub | Instructions |
| --- | ---: | ---: | ---: | ---: | ---: |
| dense Fp12 product (Fp6 Karatsuba over Fp2 Karatsuba: 3 Fp6 products of 6 Fp2 products) | 22 | 54 | 171 | 109 | 22*32 + 54*162 + 171*26 + 109*27 = **16,841** |
| Granger-Scott cyclotomic squaring (3 Fp4 squarings of 3 Fp2 squarings) | 13 | 18 | 96 | 31 | 13*32 + 18*162 + 96*26 + 31*27 = **6,665** |
| 12 Frobenius keys (section 4 of proof.md) | 0 | 4 | 1 | 8 | 4*162 + 26 + 8*27 = **890** |

Fp2 multiplication is Karatsuba (3 fp_mul, 2 fp_add, 3 fp_sub); Fp2 squaring is the
complex method (2 fp_mul, 2 fp_add, 1 fp_sub). No lazy reduction is assumed: every leaf
output is fully reduced, so each listing's input condition holds at every call.

## 4. Per-step bookkeeping

The walk step itself is the Fp12 routine of its branch plus the following scalar work. The
figures are itemised instruction estimates for straightforward code, then charged with a
margin.

Every step:

| Work | Instructions |
| --- | ---: |
| load limbs 0 and 1 of stored T = Re(c_0); branch bits (and), table index h (ubfx), DP test on (T >> 13) mod 2^60 across two limbs (lsr, orr with masked limb 1, cbz) | 12 |
| table address ((12h + t) * 384 + base: madd, lsl, add) | 4 |
| exponent update mod r: two fp_add leaf calls (alpha += A[h][t], beta += B[h][t], or alpha += alpha, beta += beta) | 52 |
| 2-cycle test: four limbs of T against T two steps back (2 ldp each, cmp, 3 ccmp, b.eq) | 10 |
| checkpoint test: four limbs of T against the checkpoint T | 10 |
| history: store T (2 stp) and rotate the three state-buffer pointers (y, alpha, beta) | 12 |
| per-walk step counter, used only modulo 64 for the checkpoint period (add, and #63, cbnz; a wrap of the 64-bit word is harmless because 2^64 is a multiple of 64), amortised checkpoint store | 6 |
| stuck-walk counter, two 64-bit words (the threshold 20/theta = 20 x 2^60 = 2^64.32 does not fit one word): ldp, adds #1, adc, stp, cmp high, ccmp low, b.hi | 7 |
| global step counter and step cap, three 64-bit words (cap 2^131.22 steps): ldp + ldr, adds #1, adcs, adc, stp + str, cmp / sbcs / sbcs against the three cap words, b.hs | 11 |
| work meter, three 64-bit words counting scalar instructions (cap 2^145.37 instructions): csel of the branch's step price, ldp + ldr, adds, adcs, adc, stp + str, cmp / sbcs / sbcs against the three cap words, b.hs | 12 |
| branch dispatch and call of the Fp12 routine (its frame is inside its own 32) | 6 |
| walk loop: next walk-state pointer, index increment, wrap test, branch | 8 |
| **itemised** | **150** |
| **charged (COMMON)** | **300** |

Multiply steps in addition, beyond the 890 of the key arithmetic:

| Work | Instructions |
| --- | ---: |
| pointers to the w^1 coefficient (a, b) and the twelve key slots | 4 |
| minimum by top limb: 12 loads, 11 x (cmp, csel value, csel index) | 45 |
| uniqueness of the top-limb minimum: 12 x (cmp, cinc), branch | 26 |
| rare full-limb comparison when two top limbs tie (taken with probability below 2^-55; bounded) | 60 |
| twist index t = (12 - kappa) mod 12 | 3 |
| **itemised, tie path included** | **138** |
| **charged (KEYMIN)** | **160** |

So the keys cost 890 + 160 = **1,050** instructions per multiply step, and

| Step | Fp12 routine | keys | bookkeeping | Instructions | Units |
| --- | ---: | ---: | ---: | ---: | ---: |
| multiply | 16,841 | 1,050 | 300 | **18,191** | **0.90955** |
| squaring | 6,665 | 0 | 300 | **6,965** | **0.34825** |

With squaring probability s = 1/2 the average step costs
c = (0.90955 + 0.34825) / 2 = **0.6289 units**.

Counter widths. The work meter counts scalar instructions in a 192-bit accumulator. Its cap is
ceil(20,000 x W_meter) < 2^145.4 instructions (proof.md section 3.4), every single addition to it is below
2^48 (a step price, a report, one verification, or one setup item), and no addition is made once the
accumulator has reached the cap, so its value never exceeds 2^145.4 + 2^48 < 2^192: it cannot wrap. A
128-bit accumulator would wrap after 2^128 instructions and is not used. The global step counter has 192
bits for a cap below 2^131.3. The stuck-walk counter has 128 bits for a threshold of 2^64.32 and is reset
at every distinguished point and restart. The report counter has 128 bits for a cap below 2^73.7; reports
are stored under two-word record indices. The restart counter (cap 2^52), the false-candidate counter (cap
2^20) and the walk index (2^20 walks) fit one word. Non-walk events add their charge to the same meter
with one three-word addition and comparison, inside their own charges (section 5). `gate.py` implements
these counters on 64-bit words with add-with-carry and word-by-word comparison and checks them against big
integers at and around their caps. `experiment.py` still names the earlier itemised figure 129 in a
comment; the charged constant COMMON = 300 that it uses is unchanged and covers the 150 itemised here.

## 5. Distinguished-point reports and setup

A report (2^65.62 expected, at most 2^73.62; proof.md sections 3.4 and 5) computes the keys if the step did
not (1,050), multiplies alpha and beta by the stored constant p^kappa * 2^256 mod r (2 x 162),
extracts the 128-bit fingerprint, hashes it into the chained multimap and writes an 80-byte
record (under 600 instructions including memory traffic): under 2,000 instructions,
**charged 10,000** (0.5 unit). Setup work (pairing-free: g is a profile constant, z is
given): conversion of z, the check z^r = 1, the table of 1,024 multipliers g^a z^b with
their 11 Frobenius images, 2^20 walk starts g^alpha z^beta, each exponentiation at most 254
squarings and 254 products (under 400 units), all together under 2^31 units.

## 6. Kernel listings (as executed by part K)

### fp_mul (158 instructions)

```
ldp x3, x4, [x2, #0]
ldp x5, x6, [x2, #16]
ldr x14, [x1, #0]
mul x7, x14, x3
mul x8, x14, x4
mul x9, x14, x5
mul x10, x14, x6
umulh x15, x14, x3
umulh x16, x14, x4
umulh x17, x14, x5
umulh x2, x14, x6
adds x8, x8, x15
adcs x9, x9, x16
adcs x10, x10, x17
adc x11, x2, xzr
ldr x14, [x1, #8]
mul x15, x14, x3
mul x16, x14, x4
mul x17, x14, x5
mul x2, x14, x6
adds x8, x8, x15
adcs x9, x9, x16
adcs x10, x10, x17
adcs x11, x11, x2
adc x12, xzr, xzr
umulh x15, x14, x3
umulh x16, x14, x4
umulh x17, x14, x5
umulh x2, x14, x6
adds x9, x9, x15
adcs x10, x10, x16
adcs x11, x11, x17
adc x12, x12, x2
ldr x14, [x1, #16]
mul x15, x14, x3
mul x16, x14, x4
mul x17, x14, x5
mul x2, x14, x6
adds x9, x9, x15
adcs x10, x10, x16
adcs x11, x11, x17
adcs x12, x12, x2
adc x13, xzr, xzr
umulh x15, x14, x3
umulh x16, x14, x4
umulh x17, x14, x5
umulh x2, x14, x6
adds x10, x10, x15
adcs x11, x11, x16
adcs x12, x12, x17
adc x13, x13, x2
ldr x1, [x1, #24]
mul x15, x1, x3
mul x16, x1, x4
mul x17, x1, x5
mul x2, x1, x6
adds x10, x10, x15
adcs x11, x11, x16
adcs x12, x12, x17
adcs x13, x13, x2
adc x14, xzr, xzr
umulh x15, x1, x3
umulh x16, x1, x4
umulh x17, x1, x5
umulh x2, x1, x6
adds x11, x11, x15
adcs x12, x12, x16
adcs x13, x13, x17
adc x14, x14, x2
adr x2, K
ldp x3, x4, [x2, #0]
ldp x5, x6, [x2, #16]
ldr x1, [x2, #32]
subs xzr, x7, #1
mul x7, x7, x1
mul x16, x7, x4
mul x17, x7, x5
mul x2, x7, x6
adcs x8, x8, x16
adcs x9, x9, x17
adcs x10, x10, x2
umulh x15, x7, x3
umulh x16, x7, x4
umulh x17, x7, x5
umulh x2, x7, x6
adc x2, x2, xzr
adds x8, x8, x15
adcs x9, x9, x16
adcs x10, x10, x17
adcs x11, x11, x2
adcs x12, x12, xzr
adcs x13, x13, xzr
adcs x14, x14, xzr
subs xzr, x8, #1
mul x8, x8, x1
mul x16, x8, x4
mul x17, x8, x5
mul x2, x8, x6
adcs x9, x9, x16
adcs x10, x10, x17
adcs x11, x11, x2
umulh x15, x8, x3
umulh x16, x8, x4
umulh x17, x8, x5
umulh x2, x8, x6
adc x2, x2, xzr
adds x9, x9, x15
adcs x10, x10, x16
adcs x11, x11, x17
adcs x12, x12, x2
adcs x13, x13, xzr
adcs x14, x14, xzr
subs xzr, x9, #1
mul x9, x9, x1
mul x16, x9, x4
mul x17, x9, x5
mul x2, x9, x6
adcs x10, x10, x16
adcs x11, x11, x17
adcs x12, x12, x2
umulh x15, x9, x3
umulh x16, x9, x4
umulh x17, x9, x5
umulh x2, x9, x6
adc x2, x2, xzr
adds x10, x10, x15
adcs x11, x11, x16
adcs x12, x12, x17
adcs x13, x13, x2
adcs x14, x14, xzr
subs xzr, x10, #1
mul x10, x10, x1
mul x16, x10, x4
mul x17, x10, x5
mul x2, x10, x6
adcs x11, x11, x16
adcs x12, x12, x17
adcs x13, x13, x2
umulh x15, x10, x3
umulh x16, x10, x4
umulh x17, x10, x5
umulh x2, x10, x6
adc x2, x2, xzr
adds x11, x11, x15
adcs x12, x12, x16
adcs x13, x13, x17
adcs x14, x14, x2
subs x15, x11, x3
sbcs x16, x12, x4
sbcs x17, x13, x5
sbcs x2, x14, x6
csel x11, x15, x11, cs
csel x12, x16, x12, cs
csel x13, x17, x13, cs
csel x14, x2, x14, cs
stp x11, x12, [x0, #0]
stp x13, x14, [x0, #16]
ret
```

### fp_add (22 instructions)

```
ldp x3, x4, [x1, #0]
ldp x5, x6, [x1, #16]
ldp x7, x8, [x2, #0]
ldp x9, x10, [x2, #16]
adds x3, x3, x7
adcs x4, x4, x8
adcs x5, x5, x9
adc x6, x6, x10
adr x11, K
ldp x7, x8, [x11, #0]
ldp x9, x10, [x11, #16]
subs x7, x3, x7
sbcs x8, x4, x8
sbcs x9, x5, x9
sbcs x10, x6, x10
csel x3, x7, x3, cs
csel x4, x8, x4, cs
csel x5, x9, x5, cs
csel x6, x10, x6, cs
stp x3, x4, [x0, #0]
stp x5, x6, [x0, #16]
ret
```

### fp_sub (23 instructions)

```
ldp x3, x4, [x1, #0]
ldp x5, x6, [x1, #16]
ldp x7, x8, [x2, #0]
ldp x9, x10, [x2, #16]
subs x3, x3, x7
sbcs x4, x4, x8
sbcs x5, x5, x9
sbcs x6, x6, x10
csetm x11, cc
adr x12, K
ldp x7, x8, [x12, #0]
ldp x9, x10, [x12, #16]
and x7, x7, x11
and x8, x8, x11
and x9, x9, x11
and x10, x10, x11
adds x3, x3, x7
adcs x4, x4, x8
adcs x5, x5, x9
adc x6, x6, x10
stp x3, x4, [x0, #0]
stp x5, x6, [x0, #16]
ret
```
