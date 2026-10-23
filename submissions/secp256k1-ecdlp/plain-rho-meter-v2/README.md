# plain-rho-meter-v2 (secp256k1-ecdlp)

Claim: time_log2 = 12900/100 (129.00). This is plain R=2048-adding parallel
Pollard rho on secp256k1, with distinguished points and replay recovery. It uses
no negation map. A round meter and event meter make the work bound hold on
every run. Success ≥ 0.39 rests on the random-mapping heuristic H1 with a
factor-2 margin. See `proof.md`.

Files:
- `proof.md`: algorithm, meter, success argument, instruction allowances,
  memory, evidence.
- `experiment.py`: standard-library Python. Reads one JSON object on stdin.
  - `{"mode":"production","track":"secp256k1"}`: production checks, exact ledger,
    and the score inequality.
  - `{"mode":"toy",...}`: end-to-end runs of the algorithm on prime-order toy
    curves of the same shape.

The offline toy rows in proof.md section 7 come from this same script, with
these inputs (one process each):

```
{"mode":"toy","track":"secp256k1","bits":B,"trials":N,"theta_log2_neg":T,"walkers":W,"R":R,"seed":S}
(B,N,T,W,R,S) in (20,400,4,16,2048,1) (24,300,5,32,2048,2) (24,300,5,32,16,6)
                 (28,200,6,128,2048,3) (28,200,6,128,2048,4)
                 (32,60,6,256,2048,5) (32,60,6,256,2048,7)
```

The toy runs use Python's `random` with fixed seeds, so they are reproducible
experimental evidence only.
