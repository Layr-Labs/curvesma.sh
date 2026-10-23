#!/usr/bin/env python3
"""The adopted escape: a look-back window of 24 class keys with 8 deterministic probes.

Provenance.  The mechanism is ported from the published description of
submission ce215151 (jungjipdo) and its sibling 793a9dfa, and from 1b6c7d78
(kaankolcu), which state the same rule: a step whose resulting class key is one
of the last 24 keys is rejected, and up to 8 table indices are tried in
deterministic order; the first probe landing outside the window is taken.

Why this and not the v9 class predicate.  v9's escape was P_ck(K), a predicate on
the class alone.  That makes the walk an AUTONOMOUS map on a finite set, hence
eventually periodic, and the v9 sweep measured the consequence: 0.095 distinct
classes per step at 2^-11.  Removing history removes re-randomisation.  A window
of 24 keys forbids every closed pattern of length <= 24 -- in particular the
2-cycles a multiplicative table produces at rate about 1/(12R) -- and its
influence decays: two trails at a common class have windows that agree after 24
steps, so their paths differ for at most 24 steps.  Crucially the walk is still
deterministic per trail GIVEN ITS OWN PAST, which is all replay needs: a recovery
replays both trails from their own starts with their own windows, so each path is
reproduced exactly and the (a,b) pairs are exact.  What the v7 objection
("history-dependent escape breaks coalescence") actually needs is per-trail
determinism under replay, not shared futures -- v9 bought the latter and lost the
former, and lost every recovery with it.

Cost.  The ring holds 24 32-bit tags (a tag comparison is 3 instructions: load,
compare, conditional move), the ring pointer update is 2, and the probe loop
control on a fire is ~20.  Priced at 130 instructions per step in ledger.py, the
same figure the two published notes carry for the 24-key window.

This module is an independent reimplementation from the published prose.  No code
was copied from either package.
"""

WINDOW = 24          # keys remembered
PROBES = 8           # deterministic probes when the first target is inside the window
TAG_BITS = 32


def tag(key):
    """32-bit tag of a class key: what the ring stores and compares."""
    k = key & ((1 << TAG_BITS) - 1)
    k ^= k >> 16
    k = (k * 0x7FEB352D) & ((1 << TAG_BITS) - 1)
    k ^= k >> 15
    return k


class WindowRing:
    """Last WINDOW class tags, plus the deterministic probe rule.

    `admit(candidates)` takes the probed target classes in order and returns
    (index_of_accepted_probe, accepted_class).  If all PROBES are inside the
    window the first is taken: the rule is total, every step is defined, and the
    walk never aborts on a blocked window (an abort would be a lost trail and a
    ledger row this package does not need).
    """

    __slots__ = ("ring", "ptr", "n", "fires", "blocked_all")

    def __init__(self):
        self.ring = [-1] * WINDOW
        self.ptr = 0
        self.n = 0
        self.fires = 0
        self.blocked_all = 0

    def start(self):
        self.ring = [-1] * WINDOW
        self.ptr = 0
        self.n = 0

    def contains(self, tg):
        r = self.ring
        for i in range(WINDOW):
            if r[i] == tg:
                return True
        return False

    def push(self, tg):
        self.ring[self.ptr] = tg
        self.ptr += 1
        if self.ptr == WINDOW:
            self.ptr = 0
        if self.n < WINDOW:
            self.n += 1

    def admit(self, candidates):
        """candidates: list of (class, tag) in deterministic probe order."""
        for p in range(PROBES):
            c, tg = candidates[p]
            if not self.contains(tg):
                if p:
                    self.fires += 1
                self.push(tg)
                return p, c
        self.blocked_all += 1
        if self.fires == 0:
            self.fires = 1
        self.push(candidates[0][1])
        return 0, candidates[0][0]

    def future_equal(self, other):
        """True iff the two rings hold the same multiset -- the resynchronisation test."""
        return sorted(self.ring) == sorted(other.ring)
