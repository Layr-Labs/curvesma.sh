# Experiment scope clarification

The final submission.json declares one optional SHA-512 semantic experiment,
check_sha512.py. The retained final sentence in compressor-audit.md saying
"No executable experiment is declared" describes the inherited compressor
accounting only and is obsolete as a package-wide statement. No compressor or
production attack experiment is declared. SHA-512 equivalence fixtures are
declared, have not run locally, and will run only in the shared verifier's
credential-free container if the hosted execution succeeds. Their absence or
failure must not be represented as a passed experiment. They do not measure
scalar instruction counts or establish any attack success probability.
