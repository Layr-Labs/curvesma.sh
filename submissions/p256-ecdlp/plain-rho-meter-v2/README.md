# plain-rho-meter-v2 (p256-ecdlp)

Claim: time_log2 = 13047/100 (130.47). This is plain R=2048-adding parallel
Pollard rho on P-256, with distinguished points and replay recovery. It uses
no negation map. A round meter and event meter make the work bound hold on
every run. Success ≥ 0.39 rests on the random-mapping heuristic H1 with a
factor-2 margin. See `proof.md`.

Files:
- `proof.md`: algorithm, meter, success argument, instruction allowances,
  memory, evidence.
- `experiment.py`: standard-library Python. Reads one JSON object on stdin.
  - `{"mode":"production","track":"p256"}`: production checks, exact ledger,
    and the score inequality.
  - `{"mode":"toy",...}`: end-to-end runs of the algorithm on prime-order toy
    curves of the same shape.

The offline toy rows in proof.md section 7 come from this same script, with
these inputs (one process each):

```
{"mode":"toy","track":"p256","bits":B,"trials":N,"theta_log2_neg":T,"walkers":W,"R":R,"seed":S}
(B,N,T,W,R,S) in (20,400,4,16,2048,1) (24,300,5,32,2048,2) (24,300,5,32,16,6)
                 (28,200,6,128,2048,3) (28,200,6,128,2048,4)
                 (32,60,6,256,2048,5) (32,60,6,256,2048,7)
```

The toy runs use Python's `random` with fixed seeds, so they are reproducible
experimental evidence only.
