# Floyd stall oracle with retained distinguished-point collision credit

## Outcome and exact input

**Candidate arithmetic score 2567/20 (128.35), supported by an explicit production success floor of 0.5461.** Scored work has log2 in (128.3126,128.3127); the fully repriced raw ceiling lies in (127.9666,127.9667). Escape accounting is priced per event, not clipped at a per-trail cap: each offset is charged, the trail ends at 3T evaluations, and the priced global allowance supplies the 1/128 residual. The floor clears 0.39 with margin 0.1561.

Input is a supplied nonidentity Q in BLS12-381 G2 and uniform z in GT. Set g=e(P0,Q), recover x from z=g^x, output [x]P0, and verify against the SAME supplied Q. Nondegeneracy makes g a generator; z=1 returns infinity; no auxiliary powers are used. All computation is H+R/20000, including every processor, setup, unsuccessful step, lookup, transfer, reconstruction and verification.

## The raw class transition and its labels

Let lambda=p mod r, B=1024 and N=(r-1)/12. Exact arithmetic checks order(lambda)=12, lambda^6=-1 and 1+lambda^4+lambda^8=0. Keys are the minimum of twelve Frobenius c1-coordinate keys with full conjugate comparison; they are not injective.

Choose retained independent uniform jump labels (a_h,b_h), with M_h=g^a_h*z^b_h, all inverse-Frobenius twists and aligned labels. For y=g^A*z^D and cached k=kappa(y), h=K(y) mod B, perform

    v = y * M_h^(p^-k)
    Av = A + lambda^-k*a_h mod r
    Dv = D + lambda^-k*b_h mod r.

Conjugate inputs give the same next class. In canonical toy coordinates the next labels are m*(A+a_h), m*(D+b_h), where m is the new alignment. The audit uses a_h=j_h-x*b_h and checks A+xD=rep_exponent after EVERY update. Its oracle accelerates only the finite tower, with24 universe and32 literal transition/twist checks per setting.

## Dual detection: pool is credit; Floyd is the stall oracle

Tortoise and hare start at the same labelled point; one tick is one tortoise and two hare transitions. Every generated DP, including both hare hops, is offered to one pool. A matching full point with nonzero aligned denominator yields x; all candidates undergo g^x=z and final supplied-Q pairing verification. The pool's immutable first record is never overwritten. If capacity is exhausted, terminate as failure. There is NO deletion, FIFO replacement, density deepening or pool clearing, even on failure.

Independently of DP marks, compare the two Floyd classes each tick. Confirm an equal key by a full Frobenius-point comparison. A nonzero denominator recovers x. For a meeting y_t=y_h^(p^j), the equation is

    x=(lambda^j*A_h-A_t)/(D_t-lambda^j*D_h) mod r.

For a fruitless confirmed relation, perform exactly ONE charged in-place offset from an independent retained B-entry table, choosing a fresh index and updating the current point and BOTH labels. Preserve every pool record. Initialize the auxiliary hare from that jumped current point and reset only the local detector counter. This is not a random start or a trail/label reset. Pre-jump reports persist; finite-table indices are not independent fresh group elements.

A point cycle with entry mu and length ell has a positive Floyd meeting at ell*max(1,ceil(mu/ell)) <= mu+ell. This statement does not depend on the DP predicate. Exhaustive functional maps on 1..4 states, all starts and every DP marking verify this fact, including DP-free cycles. Literal-coordinate DP-free fruitless witnesses with closed labels are replayed. This closes the DP-free trigger defect, not every production transport issue.

Caps: 3T evaluations per trail, 2^84 reports, 2^80 pool records. There is NO per-trail escape cap. A trail may take as many charged offsets as its 3T horizon allows; the only escape-side failure is exhausting the globally priced allowance G, which the work and tail rows already price.

## Finite evidence for charged escapes

Seven configurations use 2500 independent table/start sets each, seed 2026102310+i, three fixed coordinate-DP spacings 4,16,64 and paired controls. All exhaustive/live branches are nonzero; captures retain shared/fresh controls and256 correlated fixed-table starts, none a production theorem.

At B=8, spacing4, the v9-style Floyd/reset control recovers 613/2500 (24.52%); the DP-only retained-pool control recovers 512/2500 (20.48%); the production rule recovers 1453/2500 (58.12%). Reachable unrecovered DP-free fruitless cycles are caught 554/554. Of 6573 fruitless Floyd meetings, 3721 have DP-free terminal cycles. There are 778 recoveries using pre-escape records, no allowance failures, and a maximum of 10 escapes on a single trail with every one charged (escape additions equal escape events). The same runs under the historical hard E=4 rule recover 1331/2500 and record 358 fifth-request aborts: 122 recoveries were lost to that cap, which is now an ablation rather than production semantics.

At N=65536,B=1024 the production counts are 2258,2185,1995 at spacings4,16,64. All17500 sets finished in107023476833 nanoseconds. Histograms, witnesses, raw digests and controls are retained; none certifies sparse production transport.

## Exact pricing, storage and conditional loss rows

Use Np=ceil(1024N/1023), T0=ceil(sqrt(2Np*ln(5/3))) and T=ceil(109T0/100). An 80-term rational atanh series fixes T0 exactly. T=72864104481343941520749461235453938734. Retain every pool/label row of the preceding skeleton: sum5138. ADD Floyd comparison256, expanded detector control132 and checkpoint copy76. Sum5602, ceiling5700 instructions/call, one dense product/call.

One paid escape costs one dense product plus53290 itemized instructions, rounded to55000, charged once per emitted offset. The global allowance G=ceil(512*3T/(12B))=ceil(T/8) prices the whole escape budget the trail may use under the first-moment premise. It is both the work row and the residual row; no cheaper non-production ledger is claimed.

    Wexec-envelope = (771/200)T + (15/4)G + 2^100
    Wrepriced     = (2013/500)T + (43/10)G + 2^101
    Wscored       = (29/5)T + 2^100.

Reprice ALL nonreference rows by6/5, including Fp products and decoding. Claim slack is 0.0373..0.0374 bits; repriced headroom exceeds0.3833 bits. Analytical ceilings, not measured assembly.

Fix s=2^54, M=2^80 records. A128-byte record is fingerprint24+A32+D32+links24+metadata16. Pool+directory+scratch <=192M+2^40<2^88 bytes with96-bit links. Charge both endpoints of at most192 routing edges/report within2^26 instructions, 2^32 candidate checks at2^15 units, setup2^50 and exceptional work2^90. Ancillary work is below2^100, repriced below2^101.

Under the named report first-moment premise, K=3T+G+1 gives pool-full debit (11/10)K/(sM) and report-cap debit (11/10)K/(s*2^84). Abort-on-full, not replacement, justifies zero deletion. Independent192-bit fingerprints give candidate-cap debit2^84*M/2^(192+32). The inherited row64*s^2/N is disclosed arithmetic, not proof of cross-escape transport.

## Production success floor

The forced inverse-cycle opportunity prediction is1/(12B). Under the inherited FOUR-times conditional mean ceiling, expected escapes are at most4*3T/(12B). Markov against the priced global allowance gives

    P(escapes beyond the priced budget) <= (4*3T/(12B))/G <= 1/128.

That 1/128 is the escape row actually subtracted, because the algorithm fails closed only when the priced global budget is exhausted; per-trail escape counts are unbounded and individually charged. Starting from the conditional reference 0.5661 and subtracting 1/128, the pool debit 0.01150, the report debit 0.00072 and transport/fingerprint rows below2^-40,

    0.5661 - 1/128 - 0.01150 - 0.00072 - negligible = 0.5461 >= 0.39.

H-Rom, report first moments and sparse-DP transport remain named conditional premises; a toy margin is not a transfer bound. Within those premises the production gates pass, the four-escape clipping row is gone, and the claimed algorithm IS the charged-allowance rule.