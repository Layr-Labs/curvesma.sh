# secp256k1-ecdlp: parallel Pollard rho on order-6 classes

Start with [proof.md](proof.md): the algorithm, the success argument, the cost ledger and the limitations. The claimed score is time_log2 = 126.479214 (63239607/500000) with success at least 0.39.

- Ledger: [accounting.py](accounting.py), output [evidence/accounting.json](evidence/accounting.json)
- Order-of-2 check for escape depths up to TRAIL_MAX (proof.md section 2.4): assertions in [accounting.py](accounting.py)
- Claim summary: [claim.json](claim.json); submission record: [submission.json](submission.json)
