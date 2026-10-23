# gt-plain-rho-meter-v1 (bls12-381-fixed-g1-invert-g2)

Claim: time_log2 = 12804/100 (128.04). The attack computes g = e(P,Q0), runs
plain R=2048-multiplying parallel Pollard rho in GT to find a with z = g^a,
and returns Q = [a]Q0. Distinguished points go to 2^40 storage nodes, and a
round and event meter makes the work bound hold on every execution. See
`proof.md`.

`experiment.py` uses the standard library only and reads one JSON object on
stdin:
- `{"mode":"production"}`: parameter checks, exact ledger, score inequality.
- `{"mode":"toy",...}`: end-to-end runs in prime-order subgroups of F_q^*.

The offline rows in proof.md section 7 come from this script with these inputs:

```
{"mode":"toy","bits":B,"trials":N,"theta_log2_neg":T,"walkers":W,"R":R,"seed":S}
(B,N,T,W,R,S) in (20,400,4,16,2048,1) (24,300,5,32,2048,2) (24,300,5,32,16,6)
                 (28,200,6,128,2048,3) (28,200,6,128,2048,4)
                 (32,80,6,256,2048,5) (32,80,6,256,2048,7)
```

The seeds are fixed, so these runs are reproducible experimental evidence only.
