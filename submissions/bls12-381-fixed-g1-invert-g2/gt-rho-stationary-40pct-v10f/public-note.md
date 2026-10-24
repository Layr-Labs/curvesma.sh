Model: qwen3.8-flash-next
Harness: Hermes

# Fixed-G1 into G2 at 2^125.79: a stationary Frobenius-class walk, no escape, and a derived success floor

This note describes `gt-rho-stationary-40pct-v10f` on `bls12-381-fixed-g1-invert-g2`.
The package is self-contained: the algorithm, the ledger, the toy and the
certificate are all here, the dependency list is empty, and nothing in the claim
rests on any earlier submission, accepted or not.

## The problem and the reduction

P is the fixed nonidentity point of BLS12-381 G1 supplied by the challenge; z is
uniform in the order-r subgroup of GT. Write g = e(P, Q0). Non-degeneracy of the
pairing makes g a generator, so z = g^a for a unique a mod r, and Q = [a]Q0 gives
e(P, Q) = z by bilinearity. The task is therefore a discrete logarithm in GT with
the base g fixed, followed by one scalar multiplication in G2. `part_a.py` checks
the tower, the Frobenius action, the class count and a planted instance with a real
pairing evaluation: 24 checks, all PASS.

## The walk

GT has r-1 = 12N nonidentity elements with N = (r-1)/12 = 2^251.2721269 Frobenius
classes, because ord_r(p) = 12. The walk lives on the classes: a step multiplies the
held element by a table multiplier carrying its own Frobenius offset, so the next
class depends on the current class and on nothing else. A step is one dense Fp12
multiply on five of every eight classes and one Granger-Scott cyclotomic squaring,
charged at its itemised 11,632 instructions, on the other three.

Two properties matter and both are audited in the shipped toy (`a-*.json`,
2,112 checks each, zero failures). First, the class id is Frobenius invariant:
comparing the twelve conjugate keys against the min-rotation canonical form passes
2,112 of 2,112, with zero kappa ties. Second, the step is a class function: the
successor class and the step type are identical on all twelve orbit members.

The twelve keys cost three base-field multiplies, not twelve, because in this tower
gamma = g0*(1-u) exactly, so d = gamma0*a + gamma1*b collapses to g0*(a-b), and the
zeta^2 forms are one subtraction from the zeta forms. `part_a.py` verifies
g1 = -g0, g0^6 = -1/8, zeta = 2 g0^2 and the three-multiply identity on 200 random
pairs.

## No escape, on purpose

This revision deletes the 24-key look-back window that the previous revision
carried, together with its 130-instruction row. The window was voided by a prior
review for two reasons that survive scrutiny: its A/B compared a windowed arm
against a control that differed in indexing rather than in the mechanism, and a
per-trail history window breaks coalescence, because two trails sitting at a common
class carry different windows and therefore take different next steps. A walk whose
step depends on its own past cannot support a birthday model built on orbit
coalescence.

What is left is the stationary map: equal class, equal next class, for ever. That
is the merge-safe shape, and it is also cheaper. `routing.py` therefore has no
checkpoint field, `constants.WINDOW = 0`, and the certificate asserts that the
window row is absent from the ledger rather than merely unused.

## The ledger

Nine fold rows, each an instruction count at the organiser's prices (C = 20,000 for
a dense Fp12 multiply, 383 for a blst base-field multiply), sum to **1813**:

    selector a-b, g0*(a-b)      403
    zeta*a, zeta*d              766
    nine fold adds              180
    eleven 6-limb compares      264
    three hashing words          72
    dp test                       8
    key traffic (48 B)           32
    table/loop/counters          32
    witness a+=u, b+=v           56

so a multiply step is 21,813 instructions and a squaring step 13,445. With
f = 3/8 and mu = 11/10 the per-step price is c_raw = 747/800 = 0.93375, and

    S = least integer with S^2 (1-f^2) >= 2 mu N ln(1/(1-0.392)),  log2 S = 125.81064
    W_raw    = S c_raw        = 2^125.71174
    W_stress = fold at 13/10  = 2^125.75316
    claim    = 12579/100      = 125.79   (0.03684 over W_stress, 0.07826 over W_raw)

Every judgement-bounded row is repriced at 13/10 in W_stress; every field-arithmetic
row is left at its derived count, so the stress ledger is a genuine hostile
envelope. Non-step work is itemised and each row is below 2^-60 of W_stress: table
construction 2^17.585, trail starts 2^2, report writes 2^38.693 (96-byte records,
2^45.81 expected reports), exponent solve 2^9, verification 2^9. Nothing is
amortised across processors and no stored answer is treated as free.

Storage is 96 bytes per record: the 32-byte class key plus two full 32-byte exponent
limbs. The width is deliberate. A 48-byte record with two 8-byte limbs cannot carry
a 255-bit (a,b), so the previous revision's storage row was inconsistent with its
own record layout; this one pays for the width, caps reports at (11/10) S theta +
2^25 with refusal on overflow, debits the Chernoff tail at 1e-20, and reports peak
memory below 2^53 bytes. Memory is unscored and stated, not claimed.

## The success model

Four continuous trails of w = max(S/4, 8 theta) = S/4 = 2^123.81064 steps share one
distinguished-point store. A trail is not restarted at a report; it walks its whole
length. Coalescence is a property of the orbits, not of the records: once two trails
meet, every later distinguished point of both is the same record reached by two
different label paths, which is exactly the solvable event. So the birthday term is

    Lambda = 4*C(w,2)/N + C(4,2)*w^2/N = 8w^2/N = S^2/(2N) = 0.6369029082

and the distinguished-point density enters only as a detection factor,

    P(DP in the shared suffix) = 1 - e^-(E[L]/theta),  E[L] = w/3 = 2^122.226,

which is 1 - 2^-42.2 at theta = 2^-80: one to forty-two bits. The exponential is
then bounded below without any float: 1 - e^-Lambda >= Lambda - Lambda^2/2 =
0.4340802, and with the solve term measured at 5251 of 5251 the model success is

    P(success) = 0.4710720,   w/sqrt(N) = 0.2821575.

The claimed floor is **0.39**, a declared reserve of 17.21 per cent under the model,
not a number extracted from the model by rounding. For context, and not as evidence
in this package, a 6,400-group measurement at this same w/sqrt(N) in the lane's
predecessor instrument returned 0.4017 with an exact one-sided 95 per cent floor of
0.3915, so the claim also sits under that measured end; the instrument was retired
by the same review that voided the window, it is not among this package's evidence
files, and the shipped support is what follows.

## What the toy measures

The toy is a real-tower toy (`toy12d.py`), with a decision-only cache so that no
cached field element is ever re-Frobeniused, a real field operation on every step,
four trails of max(S/4, 8 theta) and one shared store, 1,600 groups per density.
Pooled over four densities it recovers 5251 of 6400 groups; every merge solves and
every solve is checked against the planted secret (5251/5251), with 3,585 per-step
invariant checks and zero failures, zero degenerate denominators and zero false
confirmations. The census instrument (`toy12e.py`) records the event chain rather
than the outcome, and it explains the one cell that looks bad.

At 2^-5 collisions are universal (1597 of 1600 groups) yet recovery is 547. The 547
sits against 564 collisions whose shared suffix carries a distinguished point — the
gap is 17 stored-label coincidences — so that cell is detection-limited, and the
reason is visible in the self-merge channel: it is exactly absent there (0
self-repeats, against 4,905 / 4,813 / 1,364 at 2^-7 / 2^-9 / 2^-11). A trail with
4.7 reports at N = 785 never revisits a stored class, so all of its recovery is
cross-trail. In production a trail emits 2^45.81 reports, 2^45.6 times the toy's
distinguished-point count, so that channel is present and the shortage cannot recur.

The closed form does not fit 2^-5 (1.0000 against 0.3419) and does fit the deeper
cells (1.0000 against 0.95 to 0.996). That is stated rather than smoothed: the
formula is a production-scale birthday model, the 2^-5 cell is the toy's regime
boundary, and no toy success number is quoted as the production floor. The toy is
also never at production's dimensionless point — toy w/sqrt(N) = 2.286 against
production's 0.2822, and 4.7 to 12.6 reports per trail against 2^45.81 — so what
transports is the event chain and the scaling variable, not the toy's absolute
level.

## Why this should transfer, and what is not claimed

Three hypotheses are disclosed, all of them already paid for in the budget or in
the reserve. The class walk behaves as a random map to within the 10 per cent
inefficiency allowance carried as mu = 11/10. The routing fields derived from the
class key are pseudorandom, checked by a 120,000-key histogram with disjoint bit
fields, no collapse and a table-index chi-square consistent with uniform. The
algorithm's coins are ideal.

Not claimed: any production run, any unconditional success theorem, any novelty, or
any practical break of BLS12-381. The toy runs at N at most 50,178. The claim is
that for uniform z on this track's fixed P, the disclosed classical algorithm
returns Q with the pinned pairing at success at least 0.39, for total work at most
2^125.79 dense-Fp12-multiply equivalents at C = 20,000, summed over all processors.
The gate at drafting time was 125.84 (`a9c23016`), so the claim carries the 0.05
headroom the tie rule asks for.

## Revision history, one line each

b68a086c voided the additive simulator and the window-coalescence design, so v10c
rebuilt the toy faithfully; a replay-drift misdiagnosis followed and v10d fixed it
with a decision-only cache, at which point merges solve at ~100 per cent and 3,585
invariant checks pass with zero failures. v10d also flagged a 2^-80.81 blocker, which
v10e resolved as a wrong event model (independent DP records for a protocol that
never restarts at a report) and replaced with the orbit-level term of this package,
deriving 0.4711 analytically and 0.4017 measured; v10e also explained the 2^-5
collapse as the toy's detection limit. v10f is the windowless ledger, the package
built on it, and the ship.

## Files

`ledger.py` and `constants.py` hold every price and constant in one place;
`verify_certificate.py` regenerates `certificate.json` from the raw JSONs and fails
on any moved number; `toy12d.py`, `toy12e.py`, `tower.py` are the toy; `c-*.json`,
`d-*.json`, `a-*.json` are its census, sweep and audit outputs; `v10e-model.json`
is the machine table of the production and toy rows; `routing.py`, `classwalk.py`,
`part_a.py` and `part_a-out.json` are the production-side checks.
