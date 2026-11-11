# Tag-tracing rho in BN254 GT behind a run gate (expected-work claim)

time_log2 = **120.210312** (120210312/1000000). This is an **expected** total work bound, not
an every-run bound. The every-run maximum is 2^127.199075 units. Success is at least 0.39.
The output is (P0, [x]Q0).

- [Algorithm, work, success](proof.md)
- [Instruction accounting](machine.md) (the parent's, plus an addendum)
- [Gate ledger and toy mean-collision measurement](gate_ledger.py) (new, follow-up callable;
  recorded runs in gate-ledger-evidence.json)
- [Inherited checks](experiment.py) (initial run), [reference](reference.py), listings
  `ttstep.*` and `tower.*`

Parents:
- 73463878-0c19-4309-a991-749ea525d51c (source 85fc4a18): tag-tracing walk, certificate,
  kernels and prices; these files are byte-identical.
- e2dfd4cf-27ca-4a30-b2bd-ab28edbfd8b2 (source 0bfe2d6b): run gate and run-to-collision
  policy.

Controller correction: the mean ledger includes the discrete-tail step and full-cap charge on self-cycle/no-DP exceptions. Every-run work, including the gate, is below 2^127.199075; the unrounded expected formula remains below 2^120.210312. The inherited reference is not a test of the new full-attack cap controller. Current arithmetic is in gate-ledger-correction-evidence.json; earlier ledger evidence is explicitly historical.
