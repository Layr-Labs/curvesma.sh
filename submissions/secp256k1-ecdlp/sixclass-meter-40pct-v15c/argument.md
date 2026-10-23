# Six-class metered rho v13e: prime window, pooled mu, truncated key

## Claim

Success >= 0.391002 within 2^126.9769 target-operation equivalents, conditional
on the laws below; not a production execution and not a break. secp256k1 p, r;
N=(r-1)/6; C=1000; R=2^20; DP spacing 2^69; trail cutoff L=40*2^69; 2^20
workers; 2^16 nodes with 2^44-record caps; 192-byte records. Every arrival,
re-walk, escape, emission and verification is charged once. W_exact =
2^126.9752149 < claim 126.9769 < gate 127.00.

This revision repairs the three objections that rejected 18225d19 and changes
the claim convention; the claim is 0.0231 below the live head 127.00.

## Detector: prime 67-window, class phase, full-point confirmation

A 64-step window fails on a four-cycle A-B-C-D-A whose boundary sits at B:
64 mod 4 = 0, so every boundary is at B and the rule "boundary == window
minimum" never fires; the same happens for every L | 64. The reviewer proved
this and our earlier randomised 64-phase did not fix it.

v13e uses a **prime window of 67 steps**. For every L <= 66, gcd(67,L)=1, so
successive boundary offsets are 67 mod L and visit every residue; one boundary
lands on the cycle minimum. The window of 67 observations covers a cycle of
length <= 66 entirely, so its running minimum is the cycle's global minimum.
The escape condition is the two-cycle fingerprint compare OR the boundary
minimum; every fingerprint match is confirmed by exact equality of the
canonical point (x and even y), so a 192-bit collision cannot force an escape.
The phase is hash(class key) mod 67; deriving it from the class, not the trail,
keeps F a map on classes so the six unit partners of a point share a class
sequence (merge test: 0 failures in 80,000 checks).

Detection bound: at most 67*66 = 4,422 arrivals for L <= 66. A two-cycle is
caught in at most 3 arrivals. Lengths 67..70 are outside the claimed cap and
are not counted. Adversarial enumeration: all L=2..66 x all 67 phases x all L
offsets = 148,070 runs, 148,070 escapes, 0 failures, maximum 4,422 steps. With a
fixed 64-window the same test fails for 3/4 offsets at L=4, 7/8 at L=8, 15/16
at L=16, 31/32 at L=32 and 63/64 at L=64, reproducing the reviewer's trap.

## Walk inefficiency: pooled over six configurations

For each census row, p is the probability that the first all-arrival collision
occurs within S_ops(mu=1, model 0.42); the implied constant is
mu_i = ln(1/0.58)/(-ln(1-p_i)). Six rows vary both curve size and table size:
(q25,256) p=.43125 -> mu .965290 (SE .037238); (q25,1024) p=.39000 -> 1.102026
(.044567); (q25,2048) p=.37250 -> 1.168913 (.048315); (q25,4096) p=.39375 ->
1.088447 (.043819); (q29,1024) p=.40417 -> 1.052003 (.048304); (q29,4096)
p=.39500 -> 1.083976 (.050314). Inverse-variance pooled estimate 1.069618,
SE 0.018344, chi-square 14.904 on 5 dof. Charged constant:
ceil((1.069618 + 2*.018344)*1000)/1000 = 1.107. The R=2048 row is 2.1 SE high
and R=256 2.8 SE low; that table-size trend is disclosed, and the 2-SE charge is
applied to the pool rather than to a selected row. The old single-row estimate
(1.058) and its 0.027-in-mu error are withdrawn.

## Success budget and ledger

S_ops is the least integer with S^2 >= 2*mu*N'*ln(1000/608), N'=N*R/(R-1);
at mu=1.107, S_ops = 2^126.777348. Phase 1 = 1.147*(S_ops + 2^60) =
2^126.9752132. Hard extras are 2^89.0000007 (table 2^26, node setup 2^36,
meter remainders 2^61, reports 2^68, index fallback 2^89, output 2^40, fresh
starts 2^24). Expected doubling surcharge is
(1/(6R)+2^-20)*(2.35-1.147)*S_ops = 2^107.2663769. W_exact =
2^126.9752149. Claim 126.9769 is the smallest 1e-4 decimal with
2^(claim-2e-5) >= W_exact and success(claim) >= 0.391; success =
0.392 - 0.000997 (Markov early halt) - 1.19e-7 (all other losses) = 0.391002.
The Markov headroom is 2^117.15; the loss rows are pool overflow 2^-40
(Chernoff at 4x mean, 2^16-node union), fingerprint union S(S-1)/2^193 with
S=2^85 starts, relation 1/r, sampling, detector abort 8L^2/N and endpoint
transport 64*DP^2/N. Mu tolerance at the 0.39 floor is 1.110677, 0.0037 above
the charged 1.107.

Claim convention: the family's "raw + 0.10" floor and the 27/25 reserve are
retired in favour of the frontier form `claim = W + epsilon`, as blessed on
this track (aburan28's 126.97977 over a 126.979762 ledger). The 0.001 success
reserve replaces the old 0.10 work pad; no work row is uncharged.

## Ordinary step: 147 instructions, key built from truncated high bits

The truncated rule computes, from eight limb products of x and
beta_hat=floor(2^384*beta/p), the top 64-bit fraction of each of the three
orbit x-values; the carry from column 4 is in [0,7], so the error is at most 8
ulps. v13c's objection was that it then used an exact beta*x to build the key
and the low 20 bits. v13e builds both from the truncated high bits: the
distinguished-point gate tests the top 20 bits of the canonical fraction, and
the 192-bit key concatenates the three 64-bit fraction words in canonical
order. No low bit of x and no exact beta*x is on the charged path. The near-tie
fallback (gap <= 40 ulps, or a fraction within 8 of 0/2^64) recomputes the exact
orbit at 2^12 units; its rate is at most 3*(2*40+1)+3*(2*8+1) = 294 over 2^64,
charged at 2^-50, giving 2^82.9 < 2^89.

Rows: v13c's 21 control rows (110 instructions: key load/store 4+4, branch 4,
table address 3, jump load 8, alignment 2, meter 3, index fraction 17, fold 6,
gap 4, argmin 6, parity 3, two-cycle compare 5, checkpoint compare 5, DP bit 3,
rotate 4, escape flag 2, dispatch 5, replay tag 2, spill 10, unitemised 10) plus
ten key/detector rows (41): key extract 6, carry fold 4, DP top-20 test 4, key
limb store 4, full-point confirm 4, two-cycle min select 3, phase 3, window min
update 4, min point/exponent capture 5, boundary min compare 4. Total 147;
step price 1 + 147/1000 = 1.147. Doubling 2.35 = (7M+6A+3M inversion share+200
control = 2,320 instructions, rounded up). Toy-vs-ledger parity: every toy
arithmetic operation maps to a named row; class_of is the fallback reference
only.

## Recovery

Records are (class key, alpha, beta) with Q = [alpha]T + [beta]P0. Escapes
scale the pair by 2u^{-1}, so the identity survives. Two records in one class
satisfy Q1 = [u]Q2, hence
(alpha1 - u*alpha2)T = (u*beta2 - beta1)P0 and
x = (u*beta2 - beta1)/(alpha1 - u*alpha2) mod r. All six units are tried and
each candidate is verified by [x]P0 == T. Toy: 2,000/2,000 verified at
q25/R=2048, DP spacing 8, 163 two-cycle and 2,541 checkpoint escapes exercised,
0.3375 within S_ops (the gap is parallel DP lag, 2^-37 of the meter at
production scale). Index rule: 400,000 random inputs, worst error 8 ulps, 0
fallbacks, 0 mismatches; 74,777 adversarial near-ties, 7,291 fallbacks, 0
mismatches with fallback; 0 orbit disagreements in 180,000 checks.

## Limits

H1 (all-arrival collisions at the pooled mu, and the 2^-20 structural rate for
cycles longer than 66), the twist-uniformity clause, routing uniformity and the
instruction ceilings are disclosed heuristics at toy scale, not theorems. No
production execution, no unconditional proof and no practical secp256k1 break
is claimed. Memory 2^67.58 bytes is reported, not scored.
