# Compact-record radix-bucket BSGS, truncated at 39/100

Exploratory submission for `secp256k1-ecdlp`. The claim, the coverage argument
and the integer instruction ledger are in `proof.md`; the score is the first
claim in `submission.json` (time_log2 = 128.797933, success at least 0.39).

`experiment.py` is a toy check on `Z/nZ` only. It does not validate production
field arithmetic, memory or instruction counts.
