# Deferred imaginary coordinates and a hard native-product budget

The canonical pair and k are exactly the same as before. This changes when coordinates are computed, and replaces the old class counter by a new hard product counter. It introduces no hash heuristic or unchecked exceptional case. All evidence is static and never executed on the host.

## First stage reads only initialized real coordinates

Let x=a+bu be the selected nonzero b0/b2 pair, y=T(x), and beta the primitive sixth root in Fp. For d=1/5 the nonzero Frobenius constant is (s,-s), so

 real(y)=s(a-b), imaginary(y)=(-s)(a+b).

Initially copy x to candidate k0 (96 bytes). Compute only real(y), beta*a and beta*real(y): three native Fp products and one subtraction. Write these to the real fields of candidates k1,k2,k3. Generate the other eight real fields with eight scalar add/sub/neg operations, using the same beta²=beta-1 identities. Only the real halves of all twelve candidates are initialized at this stage. Their unwritten imaginary halves are not read.

fapi_family_real_stage in canonical-min.S tests real(x) and real(y). If both are nonzero, it finds the even and odd family winners by their distinct real coordinates, and writes their indices to WMETA. If one is zero, it writes the corresponding zero-family mode and initialized zero index fields. WMETA is 24 bytes: evenk, oddk, mode (0 normal,1 even-zero,2 odd-zero). Both real components cannot vanish for nonzero x, by the nonzero diagonal formula above. Stage one reads no imaginary coordinate.

## Selected imaginary programs

In normal mode, ALWAYS compute the imaginary values of BOTH family winners into fixed initialized WIN_E and WIN_O buffers, even if their real values differ. If the even winner k=2m has m=0/3, use b or -b with no product; otherwise compute beta*b once and use beta*b, beta*b-b, -beta*b or b-beta*b as indicated by m. The even cost is I[k_even not in {0,6}] native products.

For the odd winner, first compute imaginary(y)=(-s)(a+b), one product and one addition. If its k is 1/7, use that value or its negative; otherwise compute beta*imaginary(y) once and apply the corresponding one-coordinate rotation. The odd cost is 1+I[k_odd not in {1,7}]. Thus normal imaginary cost is

 g1(X)=1+I[k_even not in {0,6}]+I[k_odd not in {1,7}] <=3.

A direct copy for m0 is an ordinary paid field-copy call when needed; it is not a free read or assignment.

In even-zero mode generate all six imaginary coordinates of that family: compute beta*b (one product), then four scalar rotations. In odd-zero mode compute imaginary(y), beta*imaginary(y) (two products and one addition), then four rotations. These give g1=1 and g1=2 respectively. The chosen zero-family's six imaginary fields are all initialized before the second stage; the other family's unwritten imaginary fields are ignored.

fapi_family_complete then selects the exact winner. Normal mode compares the two real winners, using WIN_E/WIN_O only on a real tie. Zero mode compares the six initialized imaginary values in the unique zero-real family. It copies the winning real and imaginary halves into CANONB and stores k. Its output is exactly the previous real-then-imaginary minimum. No probability claim about ties or zero components is used.

After k is known, transform the other two b pairs by the existing exact signed-constant programs. They require g2(X)=2c(k) products, where c=0 for k0/6,1 for other even k,2 for odd k. All copies/negations and structured add/sub operations are paid. Raw GT points remain complete for the reference multiplication; CANONB and saved native keys remain unchanged.

## Exact orbit averages, including tied and zero cases

The selected b0/b2 mode is constant on a Frobenius orbit, since each coefficient is multiplied/conjugated by a nonzero constant. On every nonidentity orbit, k(sigma X)=k(X)-1. The two sets of six candidate pairs are fixed sets, interchanged and cyclically relabeled as the raw orientation changes.

For a normal orbit, each family has a fixed unique minimum real value. As the input runs through all twelve orientations, the even winner's applied index takes each even value twice, and the odd winner's applied index takes each odd value twice. This holds whether the two family real minima tie or differ. Consequently sum g1 over the orbit is 12+8+8=28, mean 7/3. The algorithm deliberately computes both winner imaginaries so this bound does not assume real ties are rare.

A zero-family orbit alternates six even-zero and six odd-zero orientations. Its mean g1 is (6*1+6*2)/12=3/2. Also mean g2=8/3 over every nonidentity orbit, by the twelve c(k) values. Therefore g=g1+g2 satisfies 0<=g<=7 and has orbit mean at most5. Define g(1)=7. Uniform raw GT gives E g<=5+2/r.

This is exact finite orbit counting, not an assumption about a sampled orientation or a generic projection of GT. It works for every good generator table.

## New hard counter and failure proof

After the initial cold key, initialize B=floor(321K/64) remaining extra products. It has129 bits at this profile. Store it in THREE ordinary64-bit limbs (192 bits); do not assume it fits128 bits or that a192-bit instruction exists.

After real-stage selection and before the imaginary program, subtract its known g1 from B, failing on underflow. After final k and before the other-pair program, subtract g2, likewise failing before products are performed. product-budget.S gives the complete three-limb subtraction/status leaf: ten instructions on success, eight on failure. Argument formation, constant-pointer addressing, BL and returned-status branch fit six more, so pay16 per check. Two checks replace the former one, increasing ordinary controller849 to865. Cold program headers include their exact small product counts; all header loads/selection are paid. No expected-only work is claimed.

The existing stationary joint-state mixing bound has relative error epsilon=2^-127 for lags at least T=32. Since g is bounded by7, long-lag covariance has absolute value <=49epsilon; short-lag covariance <=49. For S=sum_(t=1)^K g(X_t),

 E S <=K(5+2/r),
 Var S <=49(2T+1)K+49epsilon K²=3185K+49epsilon K².

The exact protected-parameter Chebyshev certificate is

 Pr[S>B] <=(3185K+49epsilon K²)/(B-K(5+2/r))² <2^-100.

It includes floor rounding and the identity correction. The bound is uniform on every good table. The actual execution follows a prefix of the ideal stationary chain until first repeat or an earlier failure. Because g>=0, counter exhaustion on a prefix implies S>B on the full ideal extension. This is the same coupling discipline as for the earlier cost cap. Counter decisions inspect group-coordinate geometry only, not the coefficient histories; denominator uniformity is preserved.

On EVERY execution, first-stage products total at most3K and paid extra products at most B. Native-product work is at most

 384(3K+B)<=3078K.

A partially failed helper is conservatively also covered by a fixed cold allowance. The new counter-cap failure<2^-100 REPLACES the old canonical class-cap loss, and is subtracted separately from identity, fingerprint, table, phase, retry and dictionary losses. It does not add a second heuristic.

## Complete staged instruction certificate

Stage one has188 internal instructions on the normal path; zero paths are smaller. Stage two has60 on the normal tie path and at most109 on a zero path, whose first stage is at most26. Thus combined maximum248. canonical-min.S supplies both complete scalar leaves, including all zero/tie branches, initialized metadata, result copies and returns. Two16-instruction interfaces plus initial copy14 give248+32+14=294; retain304. They use caller registers only, no x18, internal calls, frame or hidden ordinary spills. All paired offsets are aligned/in range. No final Z flag is treated as a multiword equality.

Real rotations use8 field operations, selected imaginary rotations at most4, total12. Structured operations are real(y)'s subtraction, imaginary(y)'s addition when used, and at most four for the two other b pairs, total6. Therefore at most18 ordinary field calls at42 each are charged, including copies and negations. Worst field-program call count is12 real-stage calls+7 imaginary calls+8 remaining calls=27. Three interpreter entries and caller interfaces cost5*27+3*7+3*3=165, within the retained180 bound; parent/top preservation still fit256 wrappers.

Existing remaining canonical control133 plus imaginary-program selection16, mode guard4 and metadata context8 totals161, within192. Both minimum interfaces are separately in304; both product-cap checks are separately in controller865. Descriptors encode actual helper ABIs and fixed shared addresses. No adapter, argument, dispatch, program setup, read/write, allocation or incomplete-coordinate access is omitted. For EACH selected b-pair mode there are36normal and2zero imaginary programs,whose constants/source offsets are cold bound to that mode. All programs/code/constants and fixed WIN/WMETA scratch are constructed cold within the unchanged complete2^120 term and shared2^45 allowance.

Nonproduct canonical work<=18*42+304+192+256=1508K. Together with the hard product bound3078K, total canonical<=4586K. Complete ordinary ledger is

 4586+807 tabulation+797 joint-label hash+46 scalar+865 controller
 +192 head+224 payload+32 capped traversals+1 retry=7550.

Retain7568K,18 slack. Complete coldD=540*2^119 includes EVERY head zero, all generator/table/scalar/variant construction, attempts, validation, pairings, code/constants/descriptors, initial state/key/edge/payload, partial failures, inversion, canonical G2 output and exact verification. Peak physical capacity<2^136 bytes, zero advice. At protected C20000 and one reference per accepted transition,

 W<=(27568K+540*2^119)/20000.

Exact integer certificate (27568K+540*2^119)^125<20000^125*2^15768 establishes proposed15768/125=126.144, above diagnostic126.14337795082898. H and explicit finite banked-memory admissibility remain material. No official score, qualification, production execution, fixed-hash theorem or practical hardware availability is claimed.
