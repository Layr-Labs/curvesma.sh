# Argument v15o

1. The walk is a map on the (r-1)/6 canonical classes: key = min(x, BETA*x, BETA^2*x) is orbit-invariant, the canonical representative is (key, even y) = s*lambda^slot*A, and the increment is read from the three-slot row 3*(key&2047)+slot. Because the key and the row come from the same canonical x, all six orbit partners receive the same increment, so the quotient map is equivariant by construction.

2. Cost: each row of the composition table is a distinct symbol; the ledger recomputes them from the shipped listing and asserts SUM(rows) == W in exact integers.

3. Reports: the pool stores (key, a, c) in 12-limb nodes; a match yields (a1-a2, c1-c2) and the recovery k = -(a1-a2)/(c1-c2) mod r is verified by [k]G == Q, charged once per candidate.

4. Scope: success is a fixture statistic (toy curve q25, dp1, 6400 trials) under H1; the priced walk and the censused walk are NOT asserted to be the same program in this revision, and the note says so.
