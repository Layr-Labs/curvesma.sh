# v4 design decision — verified release candidate

Verified outcome: claim 3512/25=140.48; exact base bracket 140.2126..140.2127; stress bracket 140.3594..140.3595, leaving >0.1205 bits. Native memory 405504 bytes under 2^20. Two unit tests pass, 1503/2400 faithful toy recoveries, 30063 action comparisons, 72 old/new pipeline comparisons, 106/106 certificate checks. Full toy recapture and non-overwriting replay agree. The independent coordination gate passes. This is 4.99 bits below our accepted 145.47, not the whole-board best. Remaining uncertainty is the unchanged explicit H1, not production-size execution evidence.

Plain input is exactly T=[x]P0, including identity. No A_2, A_3, A_6, SRS secret, or n is supplied. Endomorphism only supplies [lambda*x]P0, not [x^d]P0; it does not create Cheon auxiliary powers. No unsupported powers lane is built.

Preserve accepted prior work 4dd31b95's exact keyed-SHAKE scalar-label map, canonicalizer, sign proof, initial coins, B, Floyd, entry replay and single-lap recovery. Replace scalar action's 511 affine operations (each with a Fermat inversion) by at most 255 Jacobian doublings and 256 complete mixed additions, followed by one normalization. This is an arithmetic implementation change, not a new walk or partition model. No DP store or walk table is introduced.

Kill condition: action differs from affine reference, any cap/evidence check fails, stress ledger plus 0.10-bit reserve reaches 145.47, or sufficient support cannot be packaged. Derived pre-build rows: 5375*383=2058625 multiplication instructions; 511*4096+2*8192+2*65536+131072=2371584 other instructions; total I=4430209. Fixed work 2^40 instructions; W=(6B*I+2^40)/2000. Approximate display only: base 140.2127, 1.2x other stress 140.3595; proposed claim 3512/25=140.48. Exact integer certificate is still required.

Conflicting formatting gates: Wstress>=W makes claim-log2(W)<=0.06 and claim-log2(Wstress)>=0.10 mathematically incompatible. Follow v3's transparent stronger stress reserve instead of manufacturing a certificate. This discrepancy is disclosed in the final note.

Only this new candidate directory is owned. All predecessor and sibling files/selectors remain untouched. Source and staged inventories must agree before submit.

Live lookup correction: 49d0f532-252c-4b9d-8f19-1a61e0dbb606 currently has status rejected with valid=false, not PASS. Its public note discusses the proposed repaired dense grid but its stored rejection discusses the earlier PH defect. It is not a dependency and no success or cost premise from it is imported here. Prior 4dd31b95-073f-4056-9460-c53d2f0ecaa5 was read back as review, valid=true, officialScore 145.47. Other authors' lower scored review tickets exist; 145.47 is our own accepted frontier, not the whole-board best.
