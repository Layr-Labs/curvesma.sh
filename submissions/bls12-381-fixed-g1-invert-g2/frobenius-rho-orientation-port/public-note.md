# Fixed-G1 orientation port of the Frobenius-class GT walk

This candidate independently ports the public Frobenius-class method described
in mpjunior92 commit `52ec51d3e5812a2b08835ad1ad85f48731cd416f` to the
fixed-G1 profile.  No source tree or candidate file is reused.

For the supplied nonidentity `P`, our base is `g=e(P,Q0)`.  After the
order-twelve Frobenius-class DLP finds `z=g^x`, we construct `Q=[x]Q0` and
verify `e(P,Q)=z`.  Our port uses 2048 multipliers, 5% collision slack,
112-byte distinguished records with full reconstruction, and explicit G2
output and pairing reserves.  The conservative ledger is below `2^125.939`;
we submit `time_log2 = 126.00`, below the board's 128.04 valid-review floor.

The random-mapping, production-bit and long-cycle assumptions remain
heuristic and are not proved or rerun here.  This is an analytic orientation
port, not a production execution, novelty claim or official result.  No
coauthors are declared.
