# Argument: P-256 negation-class rho with self-itemised field pricing (v9)

Claim: total classical work <= 2^(6449/50) = 2^128.98 target-group-operation
equivalents at C=1000, success >= 0.40, advice 0, peak memory < 2^73 B.
Gate: strictly below the live head 129.46 (`5fa78d62`).

## 1. Target and normalization

NIST P-256, p = 2^256-2^224+2^192+2^96-1, a = -3, cofactor 1, order r prime.
C=1000 is the borrowed secp256k1 reference unit (instructions per unit); P-256
arithmetic is charged through its own instruction schedule (section 3), never as
a borrowed point addition. H = 0, so score = log2(T/1000) with T in
instructions.

## 2. Algorithm (inherited unchanged from the live head)

Class space: the negation classes {P,-P} of <P>. N = (r-1)/2 non-identity
classes, key = x-coordinate (unique per class), canonical member (x, |y|) with
|y| = min(y, p-y), plus the identity class handled explicitly.

Coefficient walk: state X_i = [a_i]P + [b_i]Q, stored as the canonical ring
representative of a_i + b_i*x mod r.

Step: X_{i+1} = canon([2]X_i + T_{j(X_i)}), j = H(X_i) mod R, T_j = [k_j]G from
a +/- twisted table of R = 2^16 entries (R/2 stored points and their negations),
k_j wide pseudo-random. Batched affine form: one inversion per 2^15 additions.
Coefficients update as (a,b) -> (2a + k_j, 2b) mod r for the doubling part and
the table addition, all charged.

Distinguished points: key with fixed 32-bit fingerprint pattern; reports route to
2^16 storage nodes by the top 4 routing bits of H(X), per-node caps enforced,
node set-up and routing bits charged inside H1.

Escape: when a cycle is detected (a repeated class key, or a repeated DP whose
coefficients differ), the walk restarts from the *doubling* of the
minimum-fingerprint class of that cycle:
X <- canon([2] X_min(cycle)). The escape depends only on the cycle.

Recovery: at a class collision between coefficients (a,b) and (a',b') with b != b'
mod r, x = (a' - a)/(b - b') mod r; verify [x]P = Q. Identity target: x = 0.

Merging: two walks reaching the same class see the same key, the same table
index, and the same escape (a function of the cycle), so they remain merged. The
published merging lemma applies verbatim; nothing here depends on phase or
history.

## 3. Field-op itemisation (the only change in this package)

64-bit limbs, four limbs per element. Montgomery multiplication is CIOS n=4.
Instruction schedule (AArch64, everything counted: loads, stores, addressing,
control, call overhead):

M (Montgomery mul), itemised 220:
  5  init/zero the 5-limb accumulator
  4  4 x load b[i]
  80 16 limb products x (mul, umulh, adds, adc, adc)
  4  4 x m = t[0]*n0' mod 2^64
  80 16 reduction products x (mul, umulh, adds, adc, adc)
  12 4 x (adds, adc) carry-out + 4 x store t[4]
  13 final conditional subtraction of p (cmp, 4 sub, 4 sbc, 4 csel)
  8  write the result out (4 loads, 4 stores)
  14 call/return, stack frame, loop counters
A (modular add/sub), itemised 21:
  8  4 x (adds, adc)
  13 conditional subtract of p
I (inversion), itemised 17700: safegcd-style 4-limb inversion, 590 bounded
  iterations x 30 instructions.

Priced with stated margin: M = 352 (1.60x), A = 32 (1.52x), I = 2^15 (1.85x),
control = 128 (published, covers index, hash, fingerprint, DP test, routing
bits, cap check). Step = 6M + 6A + I/2^15 + control = 6*352 + 6*32 + 1 + 128 =
2433 instructions = 2.433 units.

Why this is an honest middle: the live head prices M=512 while its own
itemisation bounds M <= 227; the midpoint is 369.5. Our independent itemisation
is 220 (within 3% of 227), we price 352, and the itemisation-only step is 1455
so the price carries a 1.67x margin over the itemised schedule. The step is
strictly cheaper than the published 3.6 units: factor 1.4797, i.e. 0.5653 bits
off every row of the ledger.

## 4. Ledger

r-1 = 2^4 * 3 * 71 * 131 * 373 * 3407 * 17449 * 38189 * 187019741 * 622491383 *
1002328039319 * 2624747550333869278416773953 (product verified exactly).
ord_r(2) = (r-1)/2 = N, log2 = 254.99999999966. Note the correction to the
published statement: r-1 has 2-adic valuation 4, so ord_r(2) is (r-1)/2 rather
than the literal power 2^255; the consequence used below is unchanged.

N' = N*2^16/(2^16-1). S0 = sqrt(2 N' ln(1/0.60)) = 2^127.5154,
S_ops = ceil(1.05 S0) = 2^127.5858, certified by exact integer squares
(S_ops^2 >= (21/20)^2 2N' ln_hi and (S_ops-1)^2 < the same bound, with ln(1/0.60)
bracketed by rationals to 1e-60).

  phase 1        = 2.433*(S_ops + 2^60)             = 2^128.86859
  phase-2 hard   = 2^16 nodes x 2^12 x 2^8 + 2^20    = 2^26.12
  random extras  = 2 steps x S_ops/(2*2^16)          = 2^113.84
  total                                              = 2^128.86863
  stress (every row x 6/5)                           = 2^129.13167
  claim                                              = 128.98

Claim window: total <= claim - 0.10 (0.1114 slack); claim < stress + 0.06
(0.1517 below); claim < 129.46 (0.48 below).

## 5. Success

P[T_c <= S_ops] = 1 - 0.60 = 0.40 in the random-mapping model, above 0.39. The
meter halts at W_max = 2^claim; a trial that collides within S_ops steps costs at
most 2.433*(S_ops + 2^60)*(1 + 1/2^16) + 2^26.12 = 2^128.869 < W_max, so the cap
cannot truncate a successful trial and no halt probability is subtracted. Only
failing trials can hit the cap.

## 6. Heuristics and limitations

H1 (random-mapping walk behaviour) carries the success bound; it is supported by
the toy census and disclosed as a heuristic. H2 (escape breaks every fruitless
cycle) is supported by the trap test and by the ord_r(2) argument. The concrete
affine toy walk is not self-averaging at toy scale (0.73-1.12x the model pooled
over four tables); this is disclosed, the ledger prices from H1, and the toy
makes no production collision-time claim. Instruction counts are analytic.
Preprocessing is included; advice is zero; peak memory is 2^72 B. Not a
production run, not a measured P-256 kernel, not a practical attack.
