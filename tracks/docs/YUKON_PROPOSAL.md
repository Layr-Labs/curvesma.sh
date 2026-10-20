# Yukon challenge application

Repository: [Layr-Labs/curvesma.sh](https://github.com/Layr-Labs/curvesma.sh).
This document describes the intended import, not a completed registration.

Import the four tracks from the root manifest: secp256k1 discrete logarithm,
BLS12-381 G1 discrete logarithm with auxiliary powers, general pairing inversion,
and exact raw Miller inversion. Each uses Pareto mode with a single minimizing
`time_log2` objective and manual promotion. Success must be at least 39%; all
preprocessing, attempts and failures count. BLS scores evaluate W(n) at n=2^20.

One Astra reviewer applies `exploratory-v1`, optionally requesting one follow-up
experiment. Supported plausible heuristics may qualify; acceptance is not proof
or a practical break. There is no separate defense/adjudication service.

Start with closed, fresh registrations. Historical scalar benchmark identities,
fixed-Q0 results and old price-based scores must not carry over. The original
secp256k1 and pairing candidates remain selected for fresh evaluation; G1 and
Miller select unscored templates. Historical qualification is not current qualification.
Complete the [import and launch checks](EXECUTION.md) before opening submissions.
The UI must bind the fresh IDs and read Pareto evaluations, not old scalar scores.
