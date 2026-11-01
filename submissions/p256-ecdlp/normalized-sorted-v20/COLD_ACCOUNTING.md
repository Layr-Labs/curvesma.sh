# Conservative complete cold allowances

This analytic construction refines unused reserves, not measured production throughput. It uses the exact P-256 field primitives already supplied and checked: multiply116, square92, add23, subtract19, inverse39479, including each primitive's return. Scalars and field values are four64-bit words. No arbitrary-width operation is a unit-cost primitive.

## Affine add/double, including exceptional points

Represent a point by two canonical Montgomery field coordinates and one infinity flag. Check the flags, then coordinate equality. Copy the other point for an identity input; return infinity for opposite points or doubling with y=0. Otherwise use the ordinary slope formulas. Addition uses two differences, one inverse, two multiplies, one square and the remaining coordinate differences. Doubling uses slope=(3*x^2-3)/(2*y): two squares, two multiplies, one inverse and at most eight additions/subtractions. The constants1 and3 are stored in Montgomery representation. Thus either nontrivial path uses at most16 field calls, of which at most one is an inverse. The other fifteen calls may all be charged at116, deliberately exceeding the actual mix.

A direct wrapper places at most16 temporary fields in a fixed frame, keeps live pointers/counters in frame slots across calls, and uses no dynamic allocation. Charge16 non-field instructions per call (pointer loads/address formation, BL and live-state management),128 for flag/coordinate comparisons, branches and copies, and128 for frame/constants/result handling. Four-word equality/zero tests use bounded loads, XOR/OR and branches; coordinates and frame offsets are below4096. These512 wrapper instructions are a conservative explicit construction allowance. Required field loads/stores and RET are already in each field primitive count and are not omitted.

One complete affine operation therefore costs at most39479+15*116+512=41731, below45056. The larger45056 allowance also applies to exceptional branches, which skip arithmetic. Storage is bounded by the fixed frame and does not alter the published node-memory bound.

## Scalar multiplication and generation

Left-to-right binary multiplication over at most256 scalar bits performs at most256 doublings and256 additions, starting from infinity. Allow65536 further instructions for bit extraction, loop control, initialization, conversion/copies and output, separately from the point wrappers. Its work is at most512*45056+65536=23134208.

A generated X=aP+bQ requires at most two such scalar multiplications, one affine addition and262144 further instructions for its stored scalars, header, conversions, bookkeeping and messages not already assigned to transport. Hence it costs at most46575616, below2^26=67108864. Y=cP is cheaper and gets the same uniform charge. The separately charged scalar-sampling row is retained. All generator points, including failed trials, are charged; no preprocessing is amortized over future targets.

## Scalar-order arithmetic

For modulus r, canonical four-word modular addition can add into five words, subtract r extended by a zero word, and select the canonical result. Even charging all loads/stores and control,64 scalar instructions suffice: twelve input/modulus word loads, five additions, five subtractions, four conditional selects, four stores and at most16 pointer/control instructions total46; the64 allowance is larger. The same construction doubles a canonical value. The fifth word handles overflow; inputs below r give a result below2r and one subtraction suffices.

A256-step shift-and-add modular product doubles the accumulator each step and conditionally adds the multiplicand. Charge16 additional instructions per bit for extraction/branch/loop state and64 fixed instructions. Worst work is64+256*(2*64+16)=36928, strictly below2^16. Each addition is canonical, so induction gives the exact product modulo r. Binary exponentiation to r-2 uses at most512 such products. Including32 controller instructions per product and512 fixed instructions costs at most512*(36928+32)+512=18924032, below2^25. Fermat applies because the protected scalar group has prime order r. A zero solve denominator is rejected before inversion.

## Candidate resolution

The source-service and four count-exchange rows remain separately charged in full. The resolver receives the stored generator points and scalar labels, reconstructs its two leaf points with at most two affine operations, tests the negation relation and nonzero scalar denominator, solves for x modulo r and verifies [x]P=Q. Allocate:
- scalar inverse:2^25=33554432;
- verification scalar multiplication:23134208;
- two affine reconstructions:90112;
- all remaining scalar operations, comparisons, canonical output, local message preparation and bookkeeping:2^18=262144.

The last row can cover four full modular products at36928 each, sixteen modular additions at64 each, and more than100000 instructions for the remaining fixed-width control/copy/output operations. The solve itself needs only one product after the inverse and a few modular differences/sign changes. Transport/service counts are not reduced or charged implicitly to this row.

The sum is57040896, below2^26. Charge57040896 for every allowed candidate, including false partial matches and failed verifications. Early exits cannot increase this bound. No claim depends on a prior submission's acceptance, or on the sharper54413312 resolution allowance in public submission1d6d2b74.

cold_accounting_check.py verifies these exact integer inequalities and compares the shift-and-add arithmetic against independent modular products on boundary and deterministic randomized inputs. Those finite tests supplement, not replace, the word-level and loop-invariant arguments above. A full protected package replay and external exact-version review are still required before publication.
