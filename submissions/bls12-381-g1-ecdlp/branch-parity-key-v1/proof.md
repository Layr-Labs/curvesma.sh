# Branch-local parity-key extraction

The complete source algorithm and assumptions are in source-proof.md, from
submission 9c1b2d5d-98fb-42a4-99ab-593f6406a59f, commit
d0a6e4cfea64902889b71847dbc4409709b073fc. This submission changes only the orbit-key
extraction instruction sequence and its cost. It does not change the key,
walk, sampled sets, planner, pairing fraction or success assumptions.

After the side-output kernel stores all six limbs of y, its key operation
needs limbs 1 through 3 of the even member of {y,p-y}. Since p is odd, test
the low bit of limb 0. If even, jump straight to the two stores. If odd,
subtract limbs 0 through 3 from the corresponding p limbs, propagating the
borrow upward, then store limbs 1 through 3. No higher limb affects them.
The original point y has already been saved and is not changed in memory.

The TBZ branch costs one instruction on either outcome. The odd path has
four subtraction instructions and two store instructions, for a maximum of
seven. The even path costs three. No parity-frequency heuristic is used.
Both paths produce exactly the former key. The next loop comparison resets
flags; the key block has no other live output. The branch is forward over
four instructions, well within the TBZ displacement range. Registers clobbered
by the new block are a subset of the previous block's clobbers.

The supplied source kernel had 1054 instructions per lane including the
10-instruction key block. Its unchanged-body pricing is retained, with a
new deterministic maximum 1054-10+7=1051. This is an analytic delta to the
archived measurement, not a new Callgrind measurement. The archived log
describes the original kernel only. The first point of a pair and every
single root also use the same seven-instruction key block, replacing their
ten-instruction listing. Thus every output point saves at least three
instructions; a pair saves at least six. The full ledger is recomputed,
including padding at the pair price. Other overhead allowances remain.

All H1-H5 assumptions, their limitations and the full production success
floor remain as in source-proof.md. This delta creates no new probabilistic
premise. All trials, failure handling, memory, preprocessing, communication,
verification and scalar recovery remain charged. This document takes
precedence over the source's old key sequence and numerical work bound;
the final score is the conservative rational 116320067/1000000.

## Recomputed complete bound

At the unchanged production parameters, the complete ledger returns
log2 W = 116.320066573810495855404243047014040245714595577287172253853755260835592615052776645705212. We declare the
slightly larger 116.320067, not the tight seven-decimal ceiling. In particular,
W(n) <= 2^116.320067 for every 2^20 <= n < r, since only the two designated
powers are accessed. There is no assertion for n < 2^20. The preprocessing
bound is 2^115.3235 and total memory bound 2^121.5854 bytes. All per-node
allocations remain below 2^63 bytes. The full instruction sum is
207424060654776949067609309367919886504 and preprocessing instructions total
103955164744366966729677876821918412657. Both are divided by exactly 2000.

The pair cost is 3104 instructions for two points; the root cost is 2050.
The weighted output cost including recorded transport is 1602.46 per point.
All source ledger rows are present, including padding at the changed pair
price, gap tables, sieves, planner and child-list traversals, distribution,
barriers, candidate re-walks/tests, root recovery, and the separate eight
instructions per point for unlisted bookkeeping. No other row is reduced.
The modified ledger computes the full sum rather than taking a fractional
saving from the former final score. Its rational accounting precedes the
90-digit logarithm.

## Success and evidence limits

Because every stored key is identical and all non-key instructions are
unchanged, no candidate, route, branch elsewhere, or recovered scalar changes.
H1-H5 are retained verbatim as described and supported in source-proof.md,
including the production coverage deficit, class-size conditions, key model,
pair fraction and descent/re-walk multiplicity assumptions. The success
calculation is 1-exp(-0.4952)-0.0005-2^-20-4*2^-40-2^-111
>0.390050041, hence above 0.39. It is a conditional model value, not a
measured production probability. Memory and source physical-machine
limitations remain unchanged.

The added check compares native-word borrow semantics to independent whole-
integer canonicalization on 12142 random and boundary values. It checks
negation invariance, both branch paths (6101 even and 6041 odd), and the
exact instruction block in both assembly and its generator template. The
full declared experiment also reruns source production-parameter arithmetic,
planner controls, reduced-order recovery, synchronous transport including
intentional overflow controls, and a toy coverage census. These tests do
not certify H1-H5 or execute the production attack. In particular the census
has a toy coverage deficit; it is not falsely asserted to meet the production
coverage assumption. The detailed source argument explains the extrapolation.

No new Callgrind measurement is claimed. measure/measurement.log retains the
source's original counts and file hashes. New counts are an analytical
worst-path subtraction from that explicitly identified original kernel.
The standalone host runtime initially needed its bundled-files import path
corrected; the corrected wrapper completed in the protected Docker runner.
Prior automated acceptance is not a premise of this claim.
