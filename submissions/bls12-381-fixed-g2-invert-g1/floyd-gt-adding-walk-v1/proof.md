# bls12-381-fixed-g2-invert-g1: GT walk bound

This is an analytic accounting argument. It does not evaluate a production pairing inversion, and it does not say that BLS12-381 pairings are broken in practice.

The input is a non-identity point Q in G2 and an element z in GT. Return P in G1 with e(P, Q) = z. The identity z = 1 returns the point at infinity in G1 and does not walk.

The public generator P0 of G1 is fixed by the profile. Compute base = e(P0, Q) with one pairing. Solve z = base^a in GT by the walk below. Return P = [a]P0.

C = 20,000 instructions per reference dense Fp12 multiplication. Walk multiplications are H. They are not unpacked into instructions a second time.

## Walk

Work in the multiplicative group of the r-torsion in GT, of prime order r. Draw one 64-bit seed and build 32 jump elements R_i = base^(c_i) * z^(d_i) with the same linear congruential rule as the elliptic walk (c_i nonzero). Each exponentiation is charged at 510 dense Fp12 multiplications, above the 254 squarings and 254 multiplications of a 255-bit binary schedule. Two exponentiations per jump give 64 * 510 = 32,640 multiplications. One further 510-multiplication exponentiation verifies base^a = z.

The index uses one 64-bit limb of the Fp12 representative: h = ((word * K) >> 27) & 31, with K = 2*seed+1. It does not hash every coefficient. Charging 2,048 instructions per step covers reading several limbs.

Floyd's cycle rule is the same as on the elliptic tracks. The iteration budget is

I = 3 * ceil(sqrt(r)) = 686,966,430,457,948,734,194,560,730,881,182,400,512.

Three group multiplications per iteration give

9 * ceil(sqrt(r)) = 2,060,899,291,373,846,202,583,682,192,643,547,201,536

dense Fp12 multiplications. The budget is charged in full even if a collision appears earlier. A collision with equal exponent-of-z coefficients is a failure of this single seed. There is no restart.

## Pairings

Each of the two pairings (the base pairing, and the verification pairing e(output, input) = z) is charged at 9,343 dense Fp12 multiplications:

- Miller loop on the 64-bit parameter X = 0xd201000000010000, which has Hamming weight 6. Charge 3 products per bit after the leading bit and 2 further products per extra addition bit: 3*(64-1) + 2*(6-1) = 199.
- Final exponent, charged as 2 * 12 * bit_length(p) = 9144.

The sum is 9343. Two pairings are 18,686. This is several times a tight ate schedule. The extra products make the score worse, which is the conservative direction. H includes them.

## Instruction side and the output point

The output scalar multiplication is in G1. It is charged at the same 11,264 instructions per complete Jacobian addition used on the G1 tracks, times 520 additions, then divided by C = 20,000. That contribution is under 300 reference units.

Equality of two Fp12 elements is charged at 2,048 instructions per iteration. Coefficient updates use the same 1,536 instructions per iteration as the elliptic walk. The partition contributes 3 * 2,048 instructions per iteration. A 2^40-instruction reserve is added.

| Component | Count |
| --- | ---: |
| Walk multiplications H | 2,060,899,291,373,846,202,583,682,192,643,547,201,536 |
| Setup multiplications H | 32,640 |
| Verification exponentiation H | 510 |
| Two pairings H | 18,686 |
| Total H | 2,060,899,291,373,846,202,583,682,192,643,547,253,372 |
| Equality instructions | 1,406,907,249,577,879,007,630,460,376,844,661,556,248,576 |
| Coefficient instructions | 1,055,180,437,183,409,255,722,845,282,633,496,167,186,432 |
| Partition instructions | 4,220,721,748,733,637,022,891,381,130,533,984,668,745,728 |
| Output scalar-multiplication instructions | 5,857,280 |
| Reserve | 1,099,511,627,776 |
| Total R | 6,682,809,435,494,925,286,244,686,790,013,241,909,665,792 |

W = H + ceil(R / C) = 2,395,039,763,148,592,466,895,916,532,144,209,348,856.

The two fixed-argument tracks differ only in the output group, so their W integers differ in the low bits. Both ceilings of log2(W) at six decimal places are 130.815246. The declared score is the rational 130815246/1000000.

## Success, memory, preprocessing

The success claim is 0.39 over the uniform non-identity input point, the uniform target z, and the fresh seed. The identity target returns infinity and succeeds. The same random-mapping comparison as in the elliptic write-up applies, with the same disclosed gaps: the adding walk is not proved to be a random mapping, a dependent collision fails the seed, and the toy groups are canonical integers rather than Fp12 elements. The model failure at this budget is about 0.011. The claim stays at 0.39 unless the walk is more than about 2.4 times slower than that model.

The experiment checks the walk on the same three prime-order additive groups and checks that this W recomputes. It does not run a pairing.

Peak memory is 2^20 bytes. Preprocessing that does not depend on the challenge is ceil(2^30 / C) = 53,688 reference units, log2 rounded up to 15.712313. The pairings and the jump elements depend on the challenge and are inside W. Nonuniform advice is 0 bytes.

No Miller inversion, no number-field sieve, and no claim of a practical break is included.
