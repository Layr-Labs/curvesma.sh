# Odd-lift BSGS, x-only batches

Deterministic discrete log on NIST P-256. Odd representative \(z \in \{x, x+r\}\), babies \((2j+1)P\), giants \((2i+1)2^{128}P\), lookup on canonical affine \(x\).

The score is \(\log_2 W = 134\) with \(W = T/1000\) and \(T\) the AArch64 instruction ceiling in `proof.md`. Run the toy check with a JSON object `{"width": 4}` on stdin.
