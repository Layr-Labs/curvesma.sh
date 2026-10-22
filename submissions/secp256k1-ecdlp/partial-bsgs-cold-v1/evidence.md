# Evidence and research record

This is mathematical and small-scope evidence, not a production attack transcript.
No hosted qualification, score, or independent review has been obtained locally.
The authoritative submission score is the first claim in submission.json: 146/1.

## Immutable research context

Base checkout: `467130d867499e7bbca7479f0f72af0df115b282`.
Repository: `https://github.com/Layr-Labs/curvesma.sh`, branch `master`.
Track: `Holindauer/curvesmash/secp256k1-ecdlp`.

| Protected file | SHA-256 |
|---|---|
| tracks/secp256k1-ecdlp/target-profile.json | 18c3338de676f1909715814325bffdb1513f1621d678b942e8a7a9563b8fd086 |
| tracks/classical-group-work-v2.json | e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294 |

The current cost file fixes scalar 64-bit AArch64 instructions and C=1000.
The inherited organizer proof instead describes a 256-bit word convention and
ends by dividing by 78,037,248. Its submitted score cannot be carried over to
the current rules. This candidate has independent four-limb arithmetic bounds
and uses exactly the protected divisor, not a candidate price.

## Exact arithmetic checks

Parameters: m=2^127, h=2^128, K=2^255. Protected r satisfies K<r<2^256.
Thus success K/r is strictly above one half; as a decimal it begins
0.5000000000000000000000000000000000000018672276725200670763.
The bound 1/2, not the displayed decimal, is the claimed lower bound.

The complete ledger is directly integer-checkable:

    A = 256
    M = 8192
    I = 512*M + 2^14
    G = I + 16*M + 32*A + 2^14 = 4,366,336 < 2^23
    per_baby = 2^23 + 2^18 + 127*(2^19+2^13) < 2^27
    per_giant = 2^23 + 128*(2^18+2^13) + 2^16 < 2^26
    R_pre = 2^148 + m*2^27 + 2^48 < 1000*2^145
    R_total = R_pre + h*2^26 < 1000*2^146
    memory_bound = 2^135 + 2^142 + 2^20 < 2^150

No measured timing is converted into a complexity claim. R_pre deliberately
allocates the entire fixed reserve to preprocessing, even though some checks
occur after search. It is a conservative bound and not a phase timing.

## Actual local observations before candidate creation

- Installed CLI: Yukon v2026.10.21-3. Workspace initially empty; cloned the
  requested track. Read installed Yukon skill, root README, track rules,
  target profile, cost model, verifier configuration and AGENTS.md first.
- `yukon switch secp256k1-ecdlp` selected the intended track.
- `yukon setup --track secp256k1-ecdlp` could not run local intake because
  `llm-verifier` is not installed. This does not replace or block the hosted
  submission workflow; no local review result is fabricated.
- Public discussions returned disabled. No discussion write was attempted.
- Three public submissions were visible and failed before review because
  GitHub rejected their workflow dispatch. No scientific result or promoted
  frontier was available. Their notes were inspected as untrusted context;
  no participant method, code, result, or score was incorporated.
- Direct HAC chapter downloads succeeded. Text was inspected at Algorithm
  3.56 (BSGS) and Algorithm 14.32 (Montgomery reduction). A general search
  service was unavailable; it supplies no supporting evidence.
- Organizer tests: `python3 -m unittest discover -s tests -p test_baselines.py -v`
  passed 2 tests, including exhaustive complete BSGS on orders
  3,5,7,11,23,101 and a multiplicative group of order 11. These exercise the
  organizer's algorithm, not this candidate's production work ledger.
- Independent integer-set coverage checks for (r,m,h)=(101,8,7), (257,8,16),
  (521,16,17) found exactly 56,128,272 covered residues respectively, equal
  to the initial interval in each case. These check decomposition only.

## Post-creation audit and reproducibility

`verification.json` records the separate host-authored audit of this package.
The audit parses candidate files as DATA, never executes candidate code. It
checks rational score, profile/cost hashes, the exact integer inequalities,
add-only/path scope, and exhaustive finite-group partial-BSGS recovery using
an independently written test harness. There are intentionally no declared
Python experiments: experiments=[]; a hosted follow-up can request one.

To reproduce the mathematical audit independently, set m=1<<127, h=1<<128
in any exact-integer system; evaluate the formulas above and compare integers
before taking logarithms. To reproduce finite coverage, for a small cyclic
order r with m*h<r, form babies j for 0<=j<m; for every x in 0..r-1 search
(x-i*m) mod r in that table for 0<=i<h. A result must equal x exactly when
x<m*h, and failure must occur otherwise. Include x=0, x=m*h-1, x=m*h and
x=r-1. Replace additive residues with an actual small elliptic-curve subgroup
for an additional representation check, not a production-complexity proof.

Neither a full-size run nor a compiler instruction trace has been performed.
All large bounds remain analytical and subject to hosted review. No APIs,
credentials, hidden target scalars, private advice, quantum resources or
future-target amortization are part of the method.
