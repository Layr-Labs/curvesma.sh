#!/usr/bin/env python3
"""The class-step rule, written once and imported by both the toy and the simulation.

The walk lives on Frobenius classes: the class of y is {sigma^j(y)}, j = 0..11, where
sigma is the p-power Frobenius.  A step must be a *class function* -- two elements of the
same class must step to the same class -- and the only way to get that is to make the
multiplier index depend on the held element's own Frobenius offset.

Let kappa(y) = k be the position of the minimum of the twelve conjugate keys, i.e. the
power that carries y to the canonical representative r:

    sigma^k(y) = r,        y = sigma^(-k)(r).

A table step multiplies by a multiplier indexed by the routing twist tw of the class key.
We need the multiplier to be sigma^(-k)(m) with m = sigma^tw(m0), because then

    y' = sigma^(-k)(r) * sigma^(-k)(m) = sigma^(-k)( r * sigma^tw(m0) ),

so the *class* of y' is the class of r * sigma^tw(m0): it depends only on the class of y
(through r and tw) and not on which of the twelve images was held.  Two trails that meet
at a common class therefore take the identical next class -- the merging lemma -- and the
sqrt(12) class quotient is actually realised by the implementation.

Hence the index, with the table holding sigma^i(m0) at position i:

    i = (tw - k) mod 12.

The same rule in exponent terms: if y = g^a z^b then sigma^c(y) = g^(a p^c) z^(b p^c), so
the entry at position i supplies (u p^i, v p^i) and (a, b) -> (a + u p^i, b + v p^i).

The v5 package instead froze a per-trail twist j = mix(key, seg) % 12 at trail start and
used index (tw + j).  j is not the held offset, so the step was not a class function: 59
per cent of same-class Frobenius-image pairs stepped to different classes, the walk
effectively lived on r = 12N elements and the measured merge rate was 0.037, i.e. the
model divided by twelve.  That was our defect, found by measurement; this module is the
repair.  The field arithmetic is unchanged (the kappa minimum over the twelve keys was
already paid for in the fold row), so no ledger row moves.  sim.py keeps the old rule as
a `legacy=True` control on the same fields, tables and seeds, so the cost of the defect is
measured rather than asserted.
"""


def twist_index(routing_twist, held_twist):
    """Table twist index for a step taken from a state held at `held_twist` = kappa(y).

    Both arguments are integers mod 12.  The index is (tw - k) mod 12, the unique choice
    making the step a class function (derivation in the module docstring).
    """
    return (routing_twist - held_twist) % 12
