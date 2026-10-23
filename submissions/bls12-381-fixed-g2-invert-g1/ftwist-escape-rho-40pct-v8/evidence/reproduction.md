## Reproduction and limits of verification

From this package directory, run:

    python3 -B scripts/faithful_audit.py --write --raw-dir ../v8-raw
    python3 -B scripts/algebra_audit.py --write
    python3 -B scripts/verify_v8.py --write --raw-dir ../v8-raw

The final audit took approximately571 seconds on the captured host. Raw rows can be regenerated from fixed seeds; optional external zlib captures are bound by exact in-package digests. The arithmetic verifier rechecks every captured summary and label count against those raw rows when `--raw-dir` is supplied. Without that option it checks in-package summaries and arithmetic only. Its bootstrap mode explicitly marks documents pending; final release requires the ordinary non-bootstrap run. No organizer experiment or production attack is asserted; manifest experiments are empty. The certificate binds all evidence and source hashes, with its details split separately to keep files small. Exact algebra and budget checks do not certify Hmix-guard or Hkey.
