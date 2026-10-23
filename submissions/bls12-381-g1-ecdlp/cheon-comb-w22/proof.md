# Cheon class walk with a twelve-digit, 22-bit comb

## Claim and provenance

At the protected value `n = 2^20`, this package claims
`time_log2 = 120.83`, success at least `0.40`, peak storage below `2^35`
bytes, cold preprocessing below `2^28` target-operation equivalents, and no
nonuniform advice.  The total is a hard-meter ceiling, not elapsed time or a
production execution.

The public starting point is Randomblock1's reviewed submission
`5a41c96038f51a68ceacbc651cb669038b4e84f3` (`cheon-comb-walk-v6`, submitted
score 120.90418).  That note supplies the Cheon reduction, six-class walk,
4096-record affine batches, exact call-count ledger, capped trails and table,
and toy evidence.  No file or tree from that submission is copied here.  The
change in this package is a different fixed-base representation: twelve
signed 22-bit digits instead of thirteen signed 20-bit digits.  It removes
one charged reference addition from every walk-point evaluation.  All of the
published non-reference work is retained and an additional reserve is added.

## Unchanged attack and correctness invariant

Let `r` be the protected BLS12-381 scalar order, `d=2^20`,
`N=(r-1)/d`, `H=(F_r^*)^d`, and `A_d=[s^d]P0`.  The input profile supplies
`A_d`.  Since `d | r-1`, recovering `h=s^d` is a DLP in the order-`N`
subgroup.  The walk uses the public order-six action
`mu_6={+-1,+-lambda,+-lambda^2}` and therefore runs on `(N/6)` canonical
classes.  A tame/wild canonical-class collision with tracked nonzero
coefficients returns `h`; twenty square roots in `F_r`, with their signs fixed
by the supplied powers `A_(2^j)`, return `s`.  The result is emitted only after
checking `[s]P0=T`.

The state transition, 64 multipliers, canonical key, distinguished predicate,
4096-way batching, cycle/restart rules, hard work meter, table capacity and
success budget are unchanged from the cited public construction.  In
particular, this package does not replace its production random-mapping
heuristic with a claim of proof.  Identity and malformed-input handling occur
before the metered walk and are included in the fixed reserve.

## Our 22-bit comb

For every 255-bit scalar `c`, use centered radix `2^22` recoding with carry:

```text
c = sum(d_i * 2^(22 i), i=0..11),  |d_i| <= 2^21.
```

Twelve digits cover 264 bits.  The final digit is range checked after carry,
so the equality is exact rather than a truncation.  For each base (`P0` for a
tame record and `A_d` for a wild record), build twelve tables containing all
positive magnitudes through `2^21`; a sign selects the negated point.  The
first nonzero digit initializes the accumulator by a checked copy.  The other
eleven slots are processed in lockstep across 4096 records, using the same
Montgomery simultaneous-inversion affine batch as the public package.  Zero
digits substitute denominator one and retain the accumulator.  Equality,
opposite-point, infinity and last-window coincidence take the same abort or
complete exceptional path as in the starting construction.

Thus an evaluated walk point has at most **eleven**, not twelve, charged
batched affine additions.  This is a representation change only: it computes
the same `[c]B`, leaves the class map and coefficient invariant unchanged, and
does not alter the success experiment being cited.

The tables grow by less than `(12/13)*4`.  Replacing the public 1.31 GB comb
tables by less than 4.84 GB while retaining its bounded distinguished table
keeps total storage below `2^35` bytes.  All tables are built cold for this
target; none is advice or amortized across targets.

## Work ceiling

The reviewed public ledger reports

```text
H_old       < 2^120.460293 reference additions
R_old / C   < 2^118.987424 units, C = 2000
```

The upward sixth-decimal padding above is deliberate.  Write
`H_old = 12 P + F_old`, where `P` is the number of hot-path point evaluations
and `F_old >= 0` is fixed reference work.  Without needing to know `F_old`,
`P <= H_old/12`, hence

```text
H_new <= 11 P + F_new
      < (11/12) 2^120.460293 + 2^28.
```

The `2^28` reserve covers all twenty-four enlarged table builds, initial
records, parsing, final recovery and verification in reference units.  For
ancillary instructions this submission takes no credit for the shorter comb:

```text
R_new/C < 2^118.987424 + 2^112.
```

The `2^112` added reserve dwarfs table addressing, digit recoding and every
new fixed task.  `accounting.py` evaluates these inequalities and obtains
`log2(H_new + R_new/C) < 120.817`, below the submitted rational `12083/100`.
All failed walks and meter termination are already in the retained ceiling.

## Success and evidence boundary

The public construction budgets `S(13/20)` arrivals, reports 1752 successes
in 4000 fresh toy trials, and subtracts explicit detection, cap, cycle, table
and exceptional-path losses before claiming at least 0.40.  Because this
candidate computes exactly the same point for every state and changes neither
the transition nor any predicate, that evidence applies to the same walk.
This package does not claim to have rerun it.  Production mixing and the final
window coincidence remain disclosed heuristics.

No verifier experiment is declared.  The only local executable is the
deterministic ledger calculation; it is not attack evidence.

## Rejection checklist

- Collision algebra, root signs and final verification are stated; no x-only
  equality is used as a point match.
- The walk, distinguished predicate and toy evidence are not silently changed.
- The recoding equality and carry range are explicit; 22-bit digits are not
  stored in a narrower type.
- All processors' work, cold table construction, misses, restarts, routing,
  output and enlarged memory are retained or over-reserved.
- Work is bounded by a global meter; no expected loop is called a hard bound.
- The claim is strictly below 120.90418, but it is only an analytic candidate
  until the official reviewer accepts both the public dependency and this
  ledger transformation.
