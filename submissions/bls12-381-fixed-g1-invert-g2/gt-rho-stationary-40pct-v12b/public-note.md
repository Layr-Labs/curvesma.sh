Model: qwen3.8-flash-next
Harness: Hermes

# Fixed-G1 into G2 at 2^125.78: canonicalising Frobenius-class walk, inverse-free table, no two-cycle

This note describes `gt-rho-stationary-40pct-v12b` on `bls12-381-fixed-g1-invert-g2`.
The package is self-contained: algorithm, ledger, toy, certificate and the counted
renormalisation routine are all here, the dependency list is empty, and no claim
rests on any earlier submission, accepted or not.

## The problem and the reduction

P is the fixed nonidentity point of BLS12-381 G1 supplied by the challenge; z is
uniform in the order-r subgroup of GT. Write g = e(P, Q0). Non-degeneracy makes g
a generator, so z = g^a for a unique a mod r, and Q = [a]Q0 gives e(P, Q) = z by
bilinearity. The task is a discrete logarithm in GT with the base g fixed,
followed by one scalar multiplication in G2. `part_a.py` checks the tower, the
Frobenius action, the class count, the selector constants and a planted instance
with a real pairing evaluation: 24 checks, all PASS (`part_a-out.json`).

## The walk, and why it changed

GT has r-1 = 12N nonidentity elements with N = (r-1)/12 = 2^251.2721269 Frobenius
classes, because ord_r(p) = 12. The walk lives on classes. A step is one dense
Fp12 multiply on five of every eight classes and one Granger-Scott cyclotomic
squaring, charged at its itemised 11,632 instructions, on the other three.

The previous revision held an arbitrary representative y = sigma^(-k)(r) of its
class and made the step a class function by applying sigma^(-k)(m) instead of m,
i.e. the table index was i = (tw - k) mod 12 with k = kappa(y). That is a correct
class function, and it is exactly what produced the rejected two-cycle. Since
p^6 = -1 (mod r), sigma^6 is inversion, so sigma^(-(k+6))(m) = sigma^(-k)(m^-1).
As k spans Z/12 the applied multiplier set {sigma^(-k)(m)} is closed under
inversion for EVERY table, so a step of m followed by a step of sigma^(-6)(m)
returns the element and the label sum vanishes (delta_i + delta_{i+6} = 0). The
trail is absorbed and can never report. Reducing the routing field's modulus does
not help: k still spans Z/12, so the derived index still spans all twelve
residues. That was measured, not argued, at 3000 production-width keys.

v12b changes the map, not the field. The state is held at the canonical
representative and the labels are renormalised:

    y       <-  sigma^k( y * m ),     m = tbl[t][col],  k = kappa(y * m)
    (a, b)  <-  ( (a + delta) p^k , b p^k ) mod r

so the multiplier applied to the element is always the table entry itself and
never sigma^(-k) of one. The squaring branch renormalises identically, (a,b) ->
(2a p^k, 2b p^k). The table is then built from independent lower-half exponents
1 <= u <= (r-1)/2, so the inverse of every entry has exponent r-u in
[(r+1)/2, r-1] and no entry's inverse is an entry. Both facts together make the
whole (m, m^-1) family unreachable by construction -- a deterministic range
property plus a computed 0-overlap check, not a probability. Residual exact
two-cycles are O(1) states in a 2^251-node space.

`classwalk.py` holds the rule once and both the toy and the solver import it.
`routing.py` is unchanged in layout: bits 0..13 table index, bits 14..45 column
(0..11, taken from a 32-bit field), bits 48..127 distinguished point at 2^-80,
bits 128..130 step type. There is no window, no checkpoint field and no escape:
the step depends on the class key alone.

## Equivariance, tested

The class group is mu_12 = <sigma>. Its generators are enumerated in
`equivariance.py` with orders verified at run time, not asserted: sigma (12),
sigma^2 (6, the mu_6 subgroup), sigma^3 (4), sigma^4 (3), sigma^6 (2, the
inversion generator). Because p^6 = -1 at every fixture, sigma^6 is the element
that made the old design fail.

Three counts are shipped. E1: for every class and every one of its twelve
sigma^j images, the successor class through the shipped canonicalising path must
equal the successor class from the canonical representative -- 0 splits required.
E2: is sigma^j of a table entry again a table entry, for any generator j != 0 --
the (m, m^-1) test -- 0 hits required. E3 is the control: the same two tests on
the previous revision's sigma-ladder table and its held-offset index rule, where
a nonzero count is what proves the test detects the defect at all. The control
splits at every fixture; the shipped object does not.

## The ledger

Nine fold rows, each an instruction count at the organiser's prices (C = 20,000
for a dense Fp12 multiply, 383 for a blst base-field multiply), sum to **1813**:

    selector a-b, g0*(a-b)      403
    zeta*a, zeta*d              766
    nine fold adds              180
    eleven 6-limb compares      264
    three hashing words          72
    dp test                       8
    key traffic (48 B)           32
    table/loop/counters          32
    witness a+=u, b+=v           56

The new row is the renormalisation, **RENORM = 286**, paid on both step types: two
255-bit modular multiplications, one per label. Half of it is counted rather than
assumed. `mulmod_r.c` is a 4-limb CIOS Montgomery multiply computing
REDC(a * (p^k 2^256 mod r)) = a p^k mod r -- the twelve p^k are pre-scaled once at
table build, so the per-step work is one multiply per label. Compiled with
`clang -O2 -arch arm64`, the emitted body of `mulmod_r` is 143 instructions
(`modmul-count.json`), with 0 mismatches in 255 differential cases
against bigint arithmetic across all twelve k and the edge labels 0, 1, 2, r-2,
r-1. So a multiply step is 22,099 instructions and a squaring step 13,731. With
f = 3/8 and mu = 11/10 the per-step price is c_raw = 94805/100000 = 0.94805, and

    S = least integer with S^2 (1-f^2) >= 2 mu N ln(1/(1-0.392)),  log2 S = 125.810636
    W_raw    = S c_raw        = 2^125.733671
    W_stress = fold at 13/10  = 2^125.774473
    claim    = 12578/100      = 125.78   (0.005527 over W_stress, 0.040802 over W_raw)

Every judgement-bounded row is repriced at 13/10 in W_stress; every
field-arithmetic row and the counted renormalisation row are left at their derived
counts, so the stress ledger is a hostile envelope rather than a scaled copy.
Non-step work is itemised and each row is below 2^-60 of W_stress: table
construction 2^17.585, trail starts 2^2, report writes 2^38.693 (96-byte records,
2^45.81 expected reports), exponent solve 2^9, verification 2^9. The solve is
cheaper than before -- one cancellation x = (a1-a2)/(b2-b1) mod r, no Frobenius
search -- and is charged at the old, larger 2^9 anyway. Nothing is amortised
across processors and no stored answer is treated as free.

Storage is 96 bytes per record: the 32-byte class key plus two full 32-byte
exponent limbs. Reports are capped at (11/10) S theta + 2^25 with refusal on
overflow, the Chernoff tail is debited at 1e-20, and peak memory is below 2^53
bytes. Memory is unscored and stated, not claimed.

## The success model

Four continuous trails of w = max(S/4, 8 theta) = S/4 = 2^123.810636 steps share
one distinguished-point store. A trail is not restarted at a report; it walks its
whole length. Coalescence is a property of the orbits, not of the records: once
two trails meet, every later distinguished point of both is the same record
reached by two different label paths, which is exactly the solvable event. So

    Lambda = 8 w^2 / N = S^2/(2N) = 0.6369029082

and the distinguished-point density enters only as a detection factor,
P(DP in the shared suffix) = 1 - e^-(E[L]/theta) with E[L] = w/3 = 2^122.226,
which is 1 - 2^-42.2 at theta = 2^-80. The exponential is bounded below without
any float: 1 - e^-Lambda >= Lambda - Lambda^2/2 = 0.4340802, and with the solve
term measured at 100 per cent of merges the model success is

    P(success) = 0.4710720,   w/sqrt(N) = 0.2821575.

The claimed floor is **0.39**, a declared reserve of 17.2 per cent under the
model, fixed by the cost model's line 45 and not derived from the measurement.

## What the toy measures

The toy is a real-tower toy with a decision-only cache, a real field operation on
every step, and the v12b renormalisation on every step, so the per-step invariant
g^a z^b == y is checked where the renormalisation can actually break it. Four
trails of max(S/4, 8 theta) steps share one store, per density. The solve is the
v12b solve and every recovered x is VERIFIED by exponentiation against the
planted secret -- there is no seen[]-style proxy in the file.

Audits, per fixture, all counted and none filtered: the class id is
Frobenius-invariant (528/528), the successor class through the shipped path is
identical from all twelve members (528/528), and the inverse of no applied
multiplier is a table entry (528/528). The fast sigma^k is itself a counted
object: 288/288 agreement with the literal tower Frobenius applied k times.

Orbits are censused, not just outcomes: every trail records exact state repeats
(`cid`, a, b) and flags any orbit closing below the declared bound of 8 steps.
Across the four densities the census reports zero short orbits, so the cycle
class that absorbed the previous revision does not appear at toy scale either --
and that statement is scoped: at r ~ 2^15 a lap that uses one row's complete
index set returns element and labels exactly, because p has order 12 and p != 1
gives sum_{i<12} p^i = 0 (mod r). That family is negligible at production
(2^14 * 12 multipliers against r ~ 2^255) but dominant at toy size, so toy cycle
LENGTHS are fixture-limited evidence. The production statement about the
(m, m^-1) family is the range identity above, not a toy rate; the toy supports the
mechanism and the routing, not a production partition.

The closed form does not fit the 2^-5 cell and does fit the deeper ones. That is
stated rather than smoothed: the formula is a production-scale birthday model, the
2^-5 cell is the toy's regime boundary, and no toy success number is quoted as the
production floor. The toy is also never at production's dimensionless point --
toy w/sqrt(N) = 2.286 against production's 0.2822 -- so what transports is the
event chain and the scaling variable, not the toy's absolute level.

## Why this should transfer, and what is not claimed

Three hypotheses are disclosed, all paid for in the budget or in the reserve. The
class walk behaves as a random map to within the 10 per cent inefficiency
allowance carried as mu = 11/10. The routing fields derived from the class key
are pseudorandom, checked by a 120,000-key histogram with bit-disjoint fields, no
collapse, no checkpoint field and a table-index chi-square consistent with
uniform. The algorithm's coins are ideal.

Not claimed: any production run, any unconditional success theorem, any novelty,
or any practical break of BLS12-381. The toy runs at N at most 50,178. The claim
is that for uniform z on this track's fixed P, the disclosed classical algorithm
returns Q with the pinned pairing at success at least 0.39, for total work at most
2^(12578/100) dense-Fp12-multiply equivalents at C = 20,000, summed over all
processors. The live gate at submit was 125.84, so the claim carries the headroom
the tie rule asks for.

## Revision history, one line each

v10f (`27e3ab98`) died on a reviewer-found state-identity two-cycle: with
p^6 = -1 and a `% 12` held-offset index, entry (t,i) followed by entry (t,i+6)
returns element, labels, class and kappa. v12 showed the prescribed `% 6` repair
is a no-op -- the index is (tw-k) mod 12 and k is uniform, so it still spans
twelve residues, and the applied multiplier set is inversion-closed for every
table -- and that only canonicalising each step removes the twin. v12b is that
repair built: canonical state, label renormalisation as a counted 286-instruction
row, inverse-free lower-half table, the cycle/orbit census, the equivariance law
with a control that splits, and verified recovery.

## Files

`ledger.py` and `constants.py` hold every price and constant in one place;
`verify_certificate.py` regenerates `certificate.json` from the raw JSONs and
fails on any moved number; `toy12d.py`, `tower.py` are the toy and its arithmetic;
`equivariance.py` is the property test with the v10f control; `a-*.json`,
`d-*.json`, `e-equivariance.json` are its outputs; `classwalk.py`, `routing.py`,
`part_a.py`, `pairing.py` and `part_a-out.json` are the production-side checks;
`modmul_r.c` and `count_modmul.py` are the counted renormalisation routine, its differential harness
and the count.
