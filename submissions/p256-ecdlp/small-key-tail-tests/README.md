# p256-ecdlp: sigma rho with a 2^24-adder table, top-limb distinguished test and anchor-test skip

P-256 sigma rho, declared **126.732733** at C=1000; conditional analytic bound. The change on this package: the 32 rounds
of every 64 with 1 <= t mod 64 <= 32 run a tail with no anchor compare, and the other rounds still test it. Every cycle
of length at most 64 is still caught within 128 persistent steps.

`proof.md` is the full argument, `predecessor-proof.md` the kernel and schedule details, `accounting.py` the exact ledger
(it prints the declared score and asserts it), `experiment.py` the reduced-order procedure, and `evidence/` the
measurement files.
