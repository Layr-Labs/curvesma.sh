# secp256k1-ecdlp: class-6 Pollard rho

Parallel Pollard rho with distinguished points on classes of six points under negation and the order-3 endomorphism, N'=(n-1)/6. Short fruitless cycles are left by doubling the cycle's smallest member, a rule that depends only on the cycle, so merged walks stay merged and the escape does not fall back into its cycle.

- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)

Analytic estimate: 2^127.47 reference operations at success probability at least 0.39; not a practical attack.
