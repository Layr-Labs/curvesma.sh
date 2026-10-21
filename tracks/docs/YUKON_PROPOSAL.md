# Yukon challenge application

Repository: [Layr-Labs/curvesma.sh](https://github.com/Layr-Labs/curvesma.sh).
This document describes the intended import, not a completed registration.

Import the six tracks from the root manifest: secp256k1 and P-256 discrete
logarithms, BLS12-381 G1 discrete logarithms with and without auxiliary powers,
and both fixed-argument inversion orientations. Each uses Pareto mode with a single minimizing
`time_log2` objective and manual promotion. Success must be at least 39%; all
preprocessing, attempts and failures count. BLS scores with auxiliary powers evaluate W(n) at n=2^20.

One Astra reviewer applies `exploratory-v1`, optionally requesting one follow-up
experiment. Supported plausible heuristics may qualify; acceptance is not proof
or a practical break. There is no separate defense/adjudication service.

Start with closed, fresh registrations. Historical scalar benchmark identities,
fixed-Q0 results and old price-based scores must not carry over. The original
secp256k1 candidate remains selected for fresh evaluation; the remaining tracks
select unscored templates. Historical qualification is not current qualification.
Complete the [import and launch checks](EXECUTION.md) before opening submissions.
The UI must bind the fresh IDs and read Pareto evaluations, not old scalar scores.

The additional `bls12-381-fixed-g1-invert-g2` track supplies P and z and returns Q. It has no auxiliary powers or n parameter and selects an unscored template. Register it separately and qualify a fresh candidate; no existing scores transfer.

The manifest also includes `bls12-381-g1-ecdlp-plain`,
`bls12-381-fixed-g2-invert-g1`, and `p256-ecdlp`. These have no auxiliary powers
or n parameter and select unscored templates. Each needs a separate fresh
registration and qualification before opening; existing results do not transfer.
