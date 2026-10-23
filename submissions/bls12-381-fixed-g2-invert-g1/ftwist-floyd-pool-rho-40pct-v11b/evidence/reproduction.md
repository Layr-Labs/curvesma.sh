# Captured reproduction

Run from this package directory with PYTHONDONTWRITEBYTECODE=1.

1. `python3 scripts/test_mechanism.py`: six mechanism regressions pass.
2. `python3 scripts/faithful_audit.py --trials 2500 --raw-dir <scratch>`: seven configurations,17500 independent sets, fixed seeds2026102310+i;91742554292 nanoseconds captured. Three paired DP densities, plus old Floyd/reset, DP-only and larger-cap diagnostics.
3. `python3 scripts/verify_v11.py --write --raw-dir <scratch>` writes the exact certificate/details/log; exit2 denotes deliberate HOLD.
4. `python3 scripts/verify_v11.py --raw-dir <scratch>` must reproduce the saved certificate exactly, with the same HOLD.

Each ordered config part concatenates losslessly into one JSON record. Scratch per-trial rows are zlib-compressed JSON; SHA256 is over the uncompressed canonical bytes. They are optional reproduction data, not a hidden source dependency. Every actual-tower witness, label equation, count, branch bin and control sum is shipped. Finite evidence supports this mechanism, not production mixing, report first moments, sparse-DP transport or assembly pricing. Prior revisions/selectors were hash-frozen before work and compared after completion. No submission was attempted while the qualification gate failed.
