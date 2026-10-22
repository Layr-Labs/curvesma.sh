# bls12-381-fixed-g2-invert-g1: Pollard rho in GT on Frobenius classes

Claim: total work at most 2^126.19 dense-Fp12-multiplication equivalents (C = 20,000), success probability at least 0.39 over uniform nonidentity Q, uniform z in GT and fresh coins. Heuristics: random-mapping behaviour of the class walk, and the absence of long fruitless cycles. Both are named in section 6.

## 1. Reduction

Let r be the subgroup order, p the field modulus, P0 the G1 generator. Given (Q, z):
1. If z = 1 return P = infinity (valid exactly for z = 1). This has probability 1/r and costs nothing.
2. Compute g = e(P0, Q) with one reduced pairing. Q is not infinity and the pairing is non-degenerate on G1 x G2, so g != 1 and g generates GT (prime order r).
3. Find x with g^x = z (the DLP in GT, section 2 to 5).
4. Return P = [x]P0. By bilinearity e(P, Q) = g^x = z, for this exact Q.

## 2. Frobenius classes of GT

GT has prime order r and r divides Phi_12(p) = p^4 - p^2 + 1. Checked numerically: the multiplicative order of p modulo r is exactly 12, p^6 = -1 mod r, and p^4 - p^2 + 1 = 0 mod r. For y = g^e in GT, y^p = g^(ep), so Frobenius acts on exponents as multiplication by p mod r. The group H = {+-p^i : i = 0..5} has 12 elements (p^(i+6) = -p^i) and y^(p^6) = y^(-1) is Fp12 conjugation. For e != 0, p^i e = e mod r forces r | p^i - 1, impossible for 0 < i < 12. So every nonidentity class has exactly 12 elements and the number of classes is N' = (r-1)/12 + 1, log2 N' = 251.2721 (log2 r = 254.8571). A collision of two class members gives, with c1 = g^a1 z^b1 = c2 = g^a2 z^b2, the equation a1 + x b1 = a2 + x b2 mod r, so x = (a1 - a2)/(b2 - b1) when b1 != b2. The exponents b are close to uniform mod r, so b1 = b2 has probability about 1/r and such a collision is discarded.

## 3. Walk with conjugated jump tables

Work in the tower Fp2[w]/(w^6 - xi), xi = 1 + u. This is the profile's basis: w^6 = 1 + u, and coefficients convert by x_j = d_j + d_(j+6), y_j = d_(j+6), six additions. Write y = sum c_j w^j. Then y^(p^i) has coefficients conj^i(c_j) * gamma(i,j), with gamma(i,j) = xi^(j (p^i - 1)/6) in Fp2 (checked numerically against direct exponentiation), and y^(p^6) negates the odd coefficients.

Class key. For the 12 images y^(+-p^i), the key is the Fp-component Re(c_1) of the image's c_1. For i = 1..5 it is Re(gamma(i,1) conj^i(c_1)), two Fp mults each; the images with p^(i+6) have key minus that value mod p. The canonical image c = y^lam, lam = +-p^i, is the one with the least key (ties broken by the index (i, sign); a tie needs two distinct images with equal Re(c_1), probability below 2^-370 per step). The key depends only on the class.

Walk. Fix s = 2^20 random exponent pairs (a_j, b_j) and M_j = g^(a_j) z^(b_j). Store for every j and every lam in H the dense element M_j^(1/lam) and the pair (a_j, b_j)/lam mod r. The state is y with known exponents, y = g^a z^b. One step:
1. Compute the 12 keys, the least key and lam; j = bits 0 to 19 of the key read as an integer.
2. y' = y * M_j^(1/lam) (one dense Fp12 multiplication); (a, b) += (a_j, b_j)/lam mod r.
Then y'^lam = c * M_j, the canonical jump from the class of y. The class sequence is exactly that of the standard class walk: the next class is [c * M_j], j taken from the canonical key. The canonical image c is never built during a step.

Distinguished points. A class is distinguished when 80 further bits of its key vanish (theta = 2^-80): with the key read as an integer in [0, p), bits 0 to 19 select the jump index j and bits 20 to 99 all zero mark a distinguished point (the minimum over 12 images skews only the high bits; the low 100 bits are uniform up to 2^-280). Leaving Montgomery form costs one more Fp multiplication per step, inside the charge below. At a distinguished point, build the full image c = y^lam (at most 5 Fp2 constant multiplications), scale (a, b) by lam, hash c to 256 bits and insert (hash, a, b), 96 bytes, into one global store. A repeated hash with a different b yields x. A found x is checked by g^x = z before use.

Fruitless cycles. The map F from a class to its successor is a function of the class. A trail that enters a cycle of F never leaves it, so the rule below makes every walker leave a cycle through the same point.

Detection (fixed period alpha = 64). Each trail counts its own steps t = 0, 1, 2, ... At every step, before anything else, the walker tests its distinguished-point condition. Then: if a checkpoint key C is set and the current key equals C, a cycle is detected; else if t is a multiple of 64, it sets C := current key and W := current key; else it sets W := min(W, current key). The window of comparisons after a checkpoint is 64 steps (the equality test comes before the checkpoint is overwritten). If a trail is inside a cycle of length L <= 64, the first checkpoint after it entered lies on the cycle (at most 63 steps after entry), and the cycle returns to that point L <= 64 steps later, inside the window. So every cycle of length at most 64 is detected, and a detection is a true repeat of a class (keys of distinct classes agree with probability below 2^-370 per comparison). When detected, the trail sits at the checkpoint class and W is the least key over the whole cycle, because the window then covers exactly the L members.

Escape (a function of the cycle only). The walker keeps stepping along the cycle, with detection off, until its key equals W (at most L - 1 steps; if W is the checkpoint class, zero). At that member m* (the least-key class of the cycle) it replaces the jump index by j + 1 mod s, with the same lam and exponent update (a_(j+1), b_(j+1))/lam, and then clears C and W. The class reached is [c(m*) M_(j(m*)+1)], which depends on m* only. What depends on the walker is the detection time (the phase of t mod 64) and the entry point into the cycle; neither affects where the trail leaves. Two trails that have met are at the same class, so they follow the same F orbit: either both reach the same distinguished point on it, or both enter the same F-cycle, both pass every member of it (each walker travels a whole turn before it can detect, and the distinguished test runs at every step, so a cycle that holds a distinguished point stops both trails there), and both leave from the same m* at the same class. They stay merged. A second escape (new cycle after the escape step) is handled in the same way, since the walker keeps no memory except C, W and t.

Rates and cost per event. A 2-cycle needs the successor of X to have the same jump index (probability 1/s) and the canonicalising element at the successor to equal inversion (probability 1/12, since the 12 images of a fresh class are equally likely to hold the least key); other canonicalising elements close the cycle only if e (mu - lam) = m (lam + 1) holds for the class exponent e, probability about 1/r. So the 2-cycle rate is 1/(12 s) = 2^-23.585 per step. A cycle of length L >= 3 closes with probability above 1/r only when the exponent relation holds for every jump value as an identity, which forces each distinct jump index in the cycle to occur at least twice: at most D <= floor(L/2) distinct indices, probability at most S(L, D) s^-(L-D) with S counting partitions into blocks of size >= 2 (the Frobenius relation 1 - p^2 + p^4 = 0 gives the 3-cycle with D = 1, rate at most s^-2). Summing L = 3 to 64 gives at most 5 s^-2 = 2^-37.7 per step, which is 2^-14.1 of the 2-cycle rate. For L > 64 the same count gives at most 2^-450 per step. Total rate of cycle events rho <= 2^-23.585 (1 + 2^-14).
Cost per event, derived: from entry to the first checkpoint at most 63 steps, then L steps to the repeat, then at most L - 1 steps to m*, then the escape step, a total of 63 + L + (L - 1) + 1 <= 191 steps for L <= 64, of which the first L were fresh class visits. I charge 192 steps per event. Amortised waste per step is 192 rho <= 2^-16 (1 + 2^-14) = 2^-15.9999, charged as a factor (1 + 2^-15). The window test, the minimum and the checkpoint compare add about 50 instructions per step, inside the 8,000 charge. This follows Bos, Kaihara, Kleinjung, Lenstra, Montgomery (2010) and Bernstein, Lange, Schwabe (2011). A trail without distinguished point after 2^84 steps (16 times the mean trail length 2^80, probability e^-16) is dropped.

Starting points. Each trail begins at y = g^a z^b with fresh uniform (a, b), by a double exponentiation of at most 500 dense multiplications.

## 4. Success probability

Run until T steps have completed, with T = 1.05 k + 2^104, then finish all trails in flight (at most P 2^84 steps for P <= 2^20 processors, at most 2^104). The run is a fixed budget that halts, not an expected runtime. Here k = sqrt(-2 ln(0.61) N') = sqrt(0.98859 N'), log2 k = 125.6278.

Class visits of distinct trails behave as independent uniform draws on N' classes up to a first collision (a collision lies within one trail with probability about 2^84/k, negligible). The probability of at least one inter-trail collision among k visits is 1 - exp(-k(k-1)/(2N')), at least 0.39 for k = sqrt(0.98859 N'). Two trails that meet share all later points, so both reach the same distinguished point: van Oorschot and Wiener (1999). Trails of at most 2^84 steps have no intra-trail collision except through cycles, handled above; a random self-collision of a trail (an undetected cycle of length above 64) has total probability at most T 2^84 / (2 N') = 2^-42.6 over the run, and costs at most 2^84 steps. Non-randomness: Teske (2001) shows r-adding walks with at least 16 branches behave near-randomly; s = 2^20 branches with a hash-selected index make the deviation of the collision rate from 1/N' of order 1/s. I charge a 5% increase of k (0.07 bit) instead. Success at least 0.39 follows when the heuristic holds. The factor 1.05 covers table non-uniformity, key-hash bias, cycle losses (about 2^-16 of the steps) and the loss of a collision whose trail is dropped (2 e^-16).

## 5. Cost ledger (C = 20,000; one dense Fp12 multiplication is 1)

Per step: the multiplication 1. Key computation: 10 Fp mults at 383 instructions (COST_MODEL, blst_fp_mul) = 3,830, plus one Fp mult to leave Montgomery form = 383; 5 modular additions about 200; 6 modular negations about 150; 11 six-limb comparisons with index tracking about 450; exponent update (two 256-bit modular additions plus table loads) about 200; hash and distinguished test about 30; checkpoint compare, window minimum and loop about 70; constant loads about 200. Sum about 5,550 instructions. I charge 8,000 = 0.40 units. Step cost 1.40 units.

| Item | Charge |
| --- | --- |
| Steps T = 1.05 k + 2^104, log2 T = 125.698 | T x 1.40 x (1 + 2^-15 cycle handling) x 1.001 (capped trails and rounding) = 2^126.18508 |
| Trails: T theta x 1.001 = 2^45.70 trails | each 500 (start) + 2 (distinguished handling: Frobenius image about 7,700 instructions, scaling, 256-bit hash about 4,000, store insertion and lookup, all inside 40,000 instructions) = 502 units, total 2^54.67 |
| Jump table: s = 2^20 entries of 450 + 12 x 0.5 + 2 units, copied to at most 2^20 processors (5.4 GB each, 2^36 units) | at most 2^50 |
| Pairing e(P0, Q), [x]P0, check g^x = z, final pairing check | at most 2^16 (a pairing is about 1,000 dense multiplications; [x]P0 about 100; generous) |
| Total | W = 2^126.18508 + 2^54.7 + 2^50 + 2^16 = 2^126.18508, below 2^126.19 |

Submitted time_log2 = 126.19 (12619/100). Preprocessing (table, pairing, setup) at most 2^50, charged inside the total. Advice is 2^0 bytes.

Memory. The store holds T theta x 1.001 = 2^45.70 records of 96 bytes: 2^52.28 bytes. Tables are 2^32.4 bytes per processor, so at most 2^52.4 for 2^20 processors. Peak 2^53.34 bytes; claimed 2^53.4 bytes. Raising theta from 2^-40 to 2^-80 lowers the store by 2^40 and costs only the in-flight term (P 2^84 = 2^104 steps, 2^-21.7 of T) and the trail starts, which fall from 2^94.7 to 2^54.7 units.

## 6. Heuristics and limits

1. Random-mapping behaviour of the class walk (Teske-type), covered by a 5% margin; no proof is claimed.
2. Structured fruitless cycles of length above 64 have probability below 2^-450 per step by the index-repeat count of section 3, which treats the jump indices of distinct classes as independent uniform. Random self-collisions of a trail have total probability at most 2^-42.6 over the run. Both are heuristics.
3. The key and jump-index hashes are random-like; ties in keys have negligible probability.
4. Per-step instruction counts are analytic, based on blst's 383 instructions per Fp mult; the 8,000 charge exceeds the estimate by about 55%.
5. A toy run (experiment.py) of the same scheme, including the period-64 detection and the least-key escape, on a group of order 59,929 with a 12-element automorphism group recovers x with no wrong answers, and a merge test (200 pairs of walkers started at one class with different phases, on a 4-jump table so that cycles are frequent) gave 195 pairs ending at the same distinguished point, 5 pairs dropped and none ending apart. It checks the scheme only, not production costs.

References: van Oorschot, Wiener, J. Cryptology 12 (1999); Teske, Math. Comp. 70 (2001); Bos, Kaihara, Kleinjung, Lenstra, Montgomery (2010); Bernstein, Lange, Schwabe, PKC 2011.
