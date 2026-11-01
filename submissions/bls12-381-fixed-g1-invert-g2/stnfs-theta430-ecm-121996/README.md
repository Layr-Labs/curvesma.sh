# stnfs-theta430-ecm-121996

This candidate models STNFS for fixed-G1/recover-G2 pairing inversion on bls12-381. It computes g = e(P, Q0), obtains a = log_g(z), and returns Q = [a]Q0. The conditional modeled bound is W ≤ 2^121.996 reference units, with success 0.3922861149 under H1–H5. The complete ledger is log2(W) = 121.995884.

| Input or bound | Frozen value |
|---|---:|
| Relation-collection subregion | θ = 43/100 within a θ = 1/2 sampled pool |
| Relation alpha charged | (5.03, 0.60) |
| Canonical-work allowance | 1.05575; incidence-inclusive SE ≤ 0.02788 |
| ECM schedule search | B1 = 2^k, k = 10…36; wheels 210, 2310, 30030, 510510 |
| Relation count passed to descent | LOG2_RELS = 64.130 |
| Descent alpha and mean norms | (5.03, 1); 748.2/498.9 bits |
| Full-ledger score and modeled success | 121.995884; 0.3922861149 |

This remains a conditional exploratory estimate, not a demonstrated production attack or a claim of platform acceptance. The changed Murphy-alpha values are assumptions: the exact sums through p ≤ 2000 are (4.868096, 0.545620), while a timed-out extension reached only the partial sums (5.005442303, 0.461704663) through p ≤ 5000. Neither finite sum bounds the omitted prime tail. The modeled yield also has a thin 0.3002-bit feasibility margin against a 0.30-bit floor, with a finite-pool SE of 0.1325 bits. Read [proof.md](proof.md) for the complete argument and limitations, [DESCENT.md](DESCENT.md) for descent accounting, and [STORAGE.md](STORAGE.md) and [KERNELS.md](KERNELS.md) for resource and kernel details.

## Reproduction

Use the verifier revision pinned by the challenge README and its credential-free Docker runtime. Do not run candidate code on the host or with credentials. The declared experiment in submission.json invokes evidence16.py with samples = 16 and seed = 1. The pinned workflow is:

    llm-verifier check --root .
    llm-verifier intake --root . --out /work/intake.json
    llm-verifier prepare-runtime --intake /work/intake.json
    llm-verifier experiments --intake /work/intake.json --out /work/initial.json

An independent extended credential-free run of the frozen candidate took 106.8 seconds and returned checks_ok = true, including the declared ledger, relation, descent, ECM, and success checks. A separate run under the protected 20-second per-run timeout ended with status timeout and exit -9. The extended result does not convert that timeout into a pass under the protected runtime. Neither run validates the alpha tail or establishes production success.

The baseline source is the tree at commit 852f57de5e67dd4a066e37e3f3d926c1e4bec592, associated with submission 3a2e69ca-18c2-4f1d-b3dc-5f54c81eece0. This derivative freezes the new configuration, reruns the declared evidence, and corrects its documentation. The predecessor experiment's authorship is not established by this local record. No platform acceptance is claimed.
