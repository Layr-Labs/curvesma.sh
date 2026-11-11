# Scalar implementation and instruction bounds

All figures are analytical upper bounds in scalar AArch64 instruction
equivalents (64-bit general registers; every load, store, address
computation, branch and call counted). They are not timings and are not
derived from Python. One processor is used; parallelism would only sum to
the same total. W = R/20,000.

## Base field: M=256, A=32

Residues mod p (and mod r) are four 64-bit limbs in Montgomery form,
canonical in [0,m). Since m<2^254, sums of two residues fit in 256 bits.

Montgomery multiplication (operand scanning, fully unrolled, limbs held in
registers; AArch64 has 31 general registers):

| part | instructions |
|---|---:|
| load 8 operand limbs, load 4 modulus limbs and n0 | 13 |
| 4 schoolbook rows: 4 MUL + 4 UMULH + at most 10 ADDS/ADCS/ADC each | 72 |
| 4 reduction rounds: 1 MUL (q) + 4 MUL + 4 UMULH + at most 10 adds each | 76 |
| final SUBS/SBCS x4, borrow word, CSEL x4 | 9 |
| store 4 limbs | 4 |
| call, return, saving and restoring up to 12 callee-saved registers | 28 |
| total | 202 <= **M=256** |

Correctness is the standard REDC argument: the 8-limb product is below
m*2^256, each round adds q*m with q chosen to clear the low limb, the
result after four rounds is below 2m<2^255, so one conditional subtraction
gives the canonical residue and no ninth limb is needed. A 64x64->128
product is MUL plus UMULH; no wide unit-cost arithmetic is assumed.

Modular add: 8 loads, ADDS+3 ADCS, load modulus (4), SUBS+3 SBCS, 4 CSEL,
4 stores = 28. Subtract/negate/double are the same shape (subtract with
borrow, then add the modulus masked by the borrow). All <= **A=32**.

## Tower arithmetic

Internal representation: Fp2=Fp[u]/(u^2+1); Fp6=Fp2[v]/(v^3-xi), xi=9+u;
Fp12=Fp6[w]/(w^2-v). This is the protected field (w^6=xi, u=w^6-9); the
input coefficient list c_0..c_11 converts by a_i=(c_i+9c_{i+6})+c_{i+6}u
for the coefficient of w^i, i<6 (done once, in F), and v^k w^e = w^(2k+e).

* Fp2 mul (Karatsuba): v0=a0b0, v1=a1b1, (a0+a1)(b0+b1)-v0-v1, v0-v1:
  3M+5A.
* xi*(a+bu) = (9a-b)+(a+9b)u; 9a as three doublings and one add: 10A.
* Fp6 mul: v0=a0b0, v1=a1b1, v2=a2b2;
  c0=v0+xi((a1+a2)(b1+b2)-v1-v2); c1=(a0+a1)(b0+b1)-v0-v1+xi*v2;
  c2=(a0+a2)(b0+b2)-v0-v2+v1. Six Fp2 mul (18M+30A), 15 Fp2 add (30A),
  two xi (20A): **18M+80A**.
* Fp12 mul: v0=a0b0, v1=a1b1; c0=v0+v*v1; c1=(a0+a1)(b0+b1)-v0-v1; with
  v*(x0,x1,x2)=(xi*x2,x0,x1) (10A). Three Fp6 mul (54M+240A), five Fp6
  add (30A), one v (10A): **54M+280A = 13,824+8,960 = 22,784**.
* Fp12 square (complex method): m=a0a1; c0=(a0+a1)(a0+v*a1)-m-v*m; c1=2m.
  Two Fp6 mul (36M+160A), five Fp6 add/double (30A), two v (20A):
  **36M+210A = 9,216+6,720 = 15,936**.
* Frobenius phi^j: phi^j(sum a_i w^i)=sum conj^j(a_i)*gamma_{j,i} w^i with
  gamma_{j,i}=xi^(i(p^j-1)/6) in Fp2 (because w^6=xi and p=1 mod 6, so
  w^(p^j)=w*xi^((p^j-1)/6); conj is the Fp2 conjugation, applied iff j is
  odd). At most 6 negations and 5 Fp2 multiplications by stored constants
  (gamma_{j,0}=1): **15M+31A = 3,840+992 = 4,832**. The 11x5 constants
  (64 bytes each) are computed in setup.

Sub-results are written to fixed scratch (under 4 KiB); the per-kernel call
overhead is inside M and A, and kernel-level calls are inside the control
allowance below.

## Hash, selector, distinguished test

The class invariant t=Re(a_0) is one stored Montgomery residue (canonical,
hence an injective encoding of t). The message is the 16-byte salt followed
by the four limbs of t, taken directly as 32-bit words: 48 bytes, so one
padded SHA-256 block. Standard SHA-256 compression on W-registers: per
round Sigma1 (ROR + 2 EOR with rotated operand = 3), Ch (3), Sigma0 (3),
Maj (4), 7 additions, loads of K and W (2): 22; message schedule for 48
rounds: two small sigmas (6), 3 additions, 4 loads/stores: 13. We allow 44
per round: 64*44=2,816. Message assembly (12 word moves, padding and
length words), IV load (8), feed-forward (8), digest extraction: under
384. **Hash <= 3,200.** The Python reference hashes a different fixed
injective encoding of the same t (big-endian standard residue); the
argument does not depend on which injective encoding is used.

Selector: first digest word, UDIV+MSUB by 11, index a 11-entry table of j
and a pointer table to the constants gamma_{j,.} and c_j: under 16.
Distinguished test: mask and compare two 64-bit words: under 8.

## Remaining per-step work

C <- C*c_j mod r: one Montgomery multiplication mod r by the stored
constant: 256. Global 128-bit step counter increment and comparison with
S, 128-bit walk-length counter and comparison with L: under 24. Copying
the 48-limb result back to the state (or swapping pointers): under 100.
Four kernel calls with argument set-up, loop branch, spills: under 400.
Selector and distinguished test: under 24. Allowance **992**.

    I = 22,784 + 15,936 + 4,832 + 3,200 + 256 + 992 = 48,000.

## Fixed work F (everything outside the step loop)

* Parsing the bounded 12x64-hex input, conversion, z^r=1 (at most 508
  Fp12 multiplications): under 2^25.
* Frobenius constants: 55 Fp2 exponentiations with exponents below 2^3048
  (or iterated from gamma_{1,1} by Frobenius), c_j, Montgomery constants:
  under 2^32.
* Per walk start: 512 fresh bits with rejection sampling (at most 1024
  draws per scalar, charged 64 instructions per bit), two Fp12
  exponentiations by 254-bit scalars (at most 1,016 Fp12 multiplications,
  under 2^25 instructions) and one multiplication: under 2^26. At most
  2^20 walks: under 2^46.
* Per distinguished point: a plain linear scan of the at most 2^20 stored
  keys t (16 instructions per four-limb comparison, under 2^24 per scan;
  no hash table is needed at this size); twelve Frobenius images and
  48-limb comparisons only for a stored key equal to t (under 2^17 per
  such key); one 512-byte record copy. Under 2^26 each; at most 2^20
  points: under 2^46.
* Recovery: a handful of multiplications and one inversion mod r (r-2
  exponentiation, under 2^18), verification g^x=z (under 2^25).
* Output Q=[x]Q0: at most 254 affine doublings and 254 affine additions
  on y^2=x^3+3/(9+u) over Fp2, each one Fp2 inversion (norm, then one Fp
  inversion by exponent p-2: under 2^18), a few Fp2 multiplications and
  the usual infinity/equal-x case tests: under 2^19 each, 2^28 total;
  conversion out of Montgomery form and hex output: under 2^20.
* Failure paths end inside the same caps.

Sum < 2^50. We charge **F=2^64**.

## Storage

Record: Y (384 bytes), C, a, b (96 bytes), padded to 512. At most 2^20
records: 2^29 bytes; constants, code, stacks and
scratch under 2^20 bytes. Peak below 2^30 bytes. Advice: none.

## Total

    S = floor(sqrt(3(r-1)/16)) + 1 = 64062824933574886057232921861673675758
    R_total <= S*48,000 + 2^64
    W = R_total/20,000,  log2(W) <= 126.853862
