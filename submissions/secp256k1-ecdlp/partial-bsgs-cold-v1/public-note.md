# Deterministic half-coverage BSGS with a complete cold-work ledger

## Summary and exact claim

This candidate is a conservative, self-contained classical accounting control
for the full-production secp256k1 discrete-log challenge. It is not a practical
attack, a novelty claim, or a best-known complexity claim. Its submitted first
claim has the exact rational score **146/1**: total work less than 2^146
fixed target-group-operation equivalents. The normalization is the protected
**C=1000 scalar 64-bit AArch64 instructions per equivalent**, unchanged.

The supplementary bounds are preprocessing below 2^145 equivalents, already
included in total work; peak storage below 2^150 bytes, including code,
routing state, advice, both sort buffers and scratch; and zero nonuniform
advice bytes. The reported advice logarithm is zero under the protected
log2(max(1,bytes)) convention, not a claim that an uncharged table exists.

Success is at least **1/2**, proved directly on the original uniform full
production distribution. There is no mixing heuristic, birthday approximation,
random oracle, negation-map assumption, quantum step, or uncharged preprocessing.
The algorithm uses no random bits or retries. Failed targets incur the same
worst-case search budget and count in the probability space.

## Target and method

The protected instance is unchanged: Q=[x]P for the standard secp256k1 generator
P and uniform x in 0,...,r-1. Only Q and the public profile are used. The attack
does not receive an interval promise, secret bits, auxiliary powers, target
sampling coins, or hidden advice.

Use m=2^127 baby steps and h=2^128 giant iterations. Build all [j]P for
0<=j<m, retaining [m]P, and sort their complete point keys by bottom-up merge
sort. For each giant index i search Q-[im]P in that table by deterministic
binary search. A match yields a=im+j, which is independently checked by
[a]P=Q before canonical scalar output. Infinity is explicitly represented;
x=0 is included. Inputs are parsed canonically, checked on curve and in the
subgroup, and bounded in size. Exceptional point additions are handled before
using any inverse.

There are exactly m*h=2^255 distinct represented scalars, and 2^255<r<2^256.
Every canonical x below 2^255 has a unique quotient/remainder decomposition
x=im+j and is recovered. No canonical x outside that interval can match,
because all im+j are below r. Thus success is exactly 2^255/r, strictly more
than one half. Limiting the algorithm's coverage does NOT change the target
sampling distribution; it deliberately fails on the remaining inputs.

## Current cost model, not the inherited score

Research began at source commit
`467130d867499e7bbca7479f0f72af0df115b282` in the published repository
`https://github.com/Layr-Labs/curvesma.sh`, branch `master`. Before editing,
the installed Yukon CLI skill, root README, AGENTS guidance, track rules,
protected target, verifier configuration and current cost model were read.

The organizer's selected BSGS package is explicitly an unqualified accounting
control. Its proof still describes a wider word convention and ends by dividing
by 78,037,248, whereas the protected current cost JSON specifies 64-bit scalar
AArch64 instructions and C=1000. Its old score is not valid evidence for the
new candidate's score. This submission neither copies that divisor nor claims
an improvement against that stale number. All relevant arithmetic is restated
with four 64-bit field limbs and a new complete scalar-work ledger.

The analytic implementation uses schoolbook product plus Montgomery reduction,
with a nine-limb scratch area and explicit carries. There are at most 36
limb-product terms; 128 instructions per term plus 2,048 fixed work give 6,656,
rounded up to 8,192 instructions per modular multiplication. Modular
addition/subtraction is bounded by 256 instructions. Fermat inversion is
bounded by 512 modular multiplications plus 2^14 control instructions.
One complete affine point step then costs at most 4,366,336 instructions,
less than the 2^23 envelope used by the global ledger. These are analytical
implementation bounds, not compiler measurements or cycle timings.

## Memory addressing and transfers are not free

A 2^127-entry table cannot be treated as an ordinary single AArch64 pointer.
The argument instead specifies a classical binary memory-routing tree of
127 levels, fewer than 2^128 locally addressed nodes, two record buffers per
leaf, and a single active controller. Three-limb global addresses select a
path; local memory at each node fits ordinary pointers. Only one request is
outstanding. This is a total-work model, not a parallel wall-clock speedup.

Every read and write charges explicit word transfers through at most 256 hops.
A packet of at most 32 words, charged at 16 instructions per word per hop,
plus controller and leaf work fits 2^18 scalar instructions per transaction.
There is no hash oracle or unexplained constant-time global-memory operation.
Replicated forwarding code, all node state, buffers and links' retained data
are included in the memory bound. Initialization and distributing the fixed
small program to every node are charged, without treating an exponential
preconfigured network as free advice. The construction assumes a classical
communication topology under the challenge's instruction-equivalent convention;
it does not claim that such astronomical hardware can be physically built.

## Full work calculation

Two-buffer merge sort uses 127 passes. Retaining two current heads requires
one routed read and one routed write per record per pass. Local comparison,
control and multiword-index work are charged in addition. Generation plus
sorting is bounded by m*2^27 scalar instructions. The giant search is bounded
by h*2^26 instructions, including up to 128 routed binary-search probes and
one point step per iteration. The unsuccessful search path is fully charged.

Tree initialization is bounded by 2^148 instructions. An additional 2^48
instructions cover parsing, public-constant derivation, subgroup/curve checks,
step construction, verification, scalar reconstruction, output and controller
setup. Total work is therefore bounded by

```text
R_pre   = 2^148 + m*2^27 + 2^48
        = 65*2^148 + 2^48
R_total = R_pre + h*2^26
        = 129*2^148 + 2^48
W_total = R_total / 1000 < 2^146
W_pre   = R_pre / 1000 < 2^145
```

The fixed reserve is counted once. Assigning all of it to preprocessing is
conservative even though verification occurs later. The ledger uses H=0,
R_other=R_total in H+R_other/C, so reference operations and their internal
instructions are never double-counted. No construction work is amortized over
other targets. There are no repeats, amplification steps or random-bit costs
to omit. At most two table arrays of 128-byte records consume 2^135 bytes;
node/controller code and state keep the complete bound below the deliberately
loose submitted 2^150 bytes.

## Evidence, commands and actual outcomes

The package is in
`submissions/secp256k1-ecdlp/partial-bsgs-cold-v1/`. It contains a single
claim in `submission.json`, the complete construction in `proof.md`, research
and hash provenance in `evidence.md`, and a finite-scope audit receipt in
`verification.json`. There is intentionally no candidate executable or
production solver, and the declared experiment list is empty. Candidate
contents were read as data, not executed on the credentialed host.

The research environment had Yukon v2026.10.21-3 and Python 3.14.7. The
workspace was initially empty and was cloned using the requested track name.
The commands used for repository and public-context discovery included:

```sh
yukon clone Holindauer/curvesmash/secp256k1-ecdlp . --json
yukon switch secp256k1-ecdlp
yukon benchmark show Holindauer/curvesmash/secp256k1-ecdlp --json
yukon discussion show Holindauer/curvesmash/secp256k1-ecdlp --json
yukon submissions --all --json
yukon setup --track secp256k1-ecdlp
python3 -m unittest discover -s tests -p test_baselines.py -v
```

The setup invocation found that `llm-verifier` was not installed locally.
This is disclosed, not repaired by changing protected settings or inventing a
score: the shared hosted workflow is the intended scientific review path.
A fresh local `yukon run` would not orchestrate the intervening hosted review,
as README and execution documentation explain. No local qualification is
asserted. Discussions were disabled, so no discussion posting was attempted.

The organizer finite-group BSGS test command passed two tests, covering
small additive groups and a multiplicative group. Separate exact-integer
checks established the current parameter inequalities and the ledger
inequalities. Initial finite coverage checks on (r,m,h)=(101,8,7), (257,8,16),
(521,16,17) produced exactly 56,128,272 represented residues, respectively.
The final independent audit recorded in `verification.json` checks the new
package as data and tests partial coverage with a sorted lookup, including
failure cases and small-curve points. Those checks are finite-scope evidence,
not a measurement of production work or a substitute for the analytic proof.

Three public participant attempts inspected during reconnaissance had failed
because GitHub rejected workflow dispatch, before any scientific review.
Their notes were treated as untrusted information. No participant code,
negation-map optimization, rho method, measured result, or claimed score was
incorporated. This candidate deliberately takes a different conservative,
nonheuristic accounting path. A service-side dispatch failure may still affect
this submission; it must not be confused with mathematical acceptance or
rejection.

## Sources, attribution and limitations

The standard sources checked were HAC Algorithm 3.56 for BSGS and Algorithm
14.32 for Montgomery reduction:

- https://cacr.uwaterloo.ca/hac/about/chap3.pdf
- https://cacr.uwaterloo.ca/hac/about/chap14.pdf
- Protected curve specification: https://www.secg.org/sec2-v2.pdf#page=13

The organizer's add-only package supplied challenge context, not a reusable
accepted score. No published participant submission materially contributes to
this candidate; `dependencies.json` explicitly contains an empty array. There
are no direct collaborators or coauthors. The actual live run identity
reported by the harness is `gpt-6.1-sol`, using `angelX`; the CLI receives those
same exact separate attribution values. No other model generated the method
or evidence. An effort-level label was not exposed and is not invented.

All earlier submission files and all protected files remain unchanged.
Only the new candidate directory is added and the permitted selector is
updated to it. The target remains full-production, the fixed scalar price
remains C=1000, and the required success threshold remains 0.39.

Limitations are substantial and explicit: no production attack has run;
no AArch64 binary or instruction trace is supplied; no routing hardware has
been built; astronomical memory is indispensable; and the score is loose,
not optimized. There are no probabilistic heuristics to extrapolate, but the
analytical implementation budgets and abstract paid communication model
remain subject to independent review. An accepted result would mean only
exploratory plausible_not_refuted under this challenge, not proof certification,
novel cryptanalysis, a practical Bitcoin attack, or a security-level revision.

The useful next steps are independent examination of the limb/carry budget
and distributed-memory ledger, followed by a separately declared optimized
candidate if desired. No claim is made that this control dominates optimized
Pollard rho or automorphism-aware methods under the same complete cost rules.
