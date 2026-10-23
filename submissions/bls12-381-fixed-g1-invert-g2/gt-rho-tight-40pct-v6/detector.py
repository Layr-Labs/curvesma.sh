#!/usr/bin/env python3
"""The cycle detector: one implementation, imported by toy.py and sim.py.

Nothing here is a model of the rule -- the same object is executed by the toy and by
the simulation, and the per-step cost of every branch is a named row in ledger.py.

Three parts, all class-deterministic (the rule reads only the class key and the step
index, both class functions), so two trails sitting at a common class have identical
history from that class on and take identical detector actions.  That is the merging
lemma the escape argument needs.

  1. two-cycle rule.  fp(K(t)) is compared with fp(K(t-2)); a fingerprint match is
     confirmed against the stored full key at t-2 before anything happens, so a
     fingerprint false match costs one comparison and can never cause a wrong action.

  2. checkpoint rule.  Checkpoints sit at steps t = phi (mod W_CAP), where W_CAP = 64
     and phi is derived from the hash of the state at the current window start (so
     every phase in [0, W_CAP) is reachable, and re-deriving phi after an escape makes
     repeated escapes non-deterministic).  At a boundary the window's minimum class
     key is compared with the previous window's minimum, and the boundary key with
     the previous boundary key.  Correctness of this rule does not depend on the
     phase: a closed cycle of length L <= W_CAP lies entirely inside every window of
     W_CAP consecutive steps, so the window minimum equals the cycle minimum at every
     boundary and consecutive minima are equal -- the rule fires within at most
     2*W_CAP steps of entering any such cycle.  The phase is not what makes it fire;
     the phase is what stops an escape from landing back on the same window position
     forever (the fixed-phase trap: a 4-cycle A->B->C->D->A with the boundary at B
     keeps the same phase after the escape and re-enters the same trap).

  3. escape = double traversal-min, a CLASS MAP.  Index the multiplier table by the
     current class key K -- the traversal minimum of the twelve conjugate keys, i.e.
     the class's own minimum -- and multiply the held state by that multiplier twice,
     then clear the fingerprint ring, the window minimum and the boundary key, and
     re-derive the phase from the new class key's hash.  The action reads only K, so
     two trails sitting at a common class perform the identical action; nothing about
     it is history dependent.  The class of (c * m(K)^2) is determined by the class of
     c alone, and one escape leaves the closed pattern unless m(K)^2 is the identity,
     i.e. unless u(K) + v(K)*s = 0 mod r for the trail's own secret s, which has
     probability <= 1/r over the table draw.

     Trigger (a), the lag-2 rule, is itself merge consistent: merged trails follow the
     same class path, so a lag-2 recurrence seen by one is seen by the other at the
     same class.  Trigger (b), the window rule, is *not* merge consistent -- it reads
     a window the two trails did not share -- so a trail that fires while another sits
     at the common class splits away from it.  That is a repair, not a merge
     detection: the split cost is carried in the loss budget, and the merge detection
     itself is the distinguished-point store, which needs no detector agreement.

  ESC_CAP escapes are allowed per trail; the next trigger aborts the trail, which is
  charged as a lost trail in the loss budget rather than left unbounded.
"""

W_CAP = 64          # checkpoint window length, in steps
ESC_CAP = 4         # escapes allowed per trail before the trail is abandoned
FP_BITS = 62        # fingerprint width (production key is 381 bits, so fp(K) does not determine K)


def fp(key, bits=FP_BITS):
    """Fingerprint of a class key: the low `bits` bits, independently mixed upstream.

    The machining cost of forming it is a hashing row in the ledger.  It is *not*
    claimed injective at production width: a match is always confirmed on the full key.
    """
    return key & ((1 << bits) - 1)


class Detector:
    """Per-trail detector state.  `observe` is called once per step, after the key."""

    def __init__(self, w_cap=W_CAP, esc_cap=ESC_CAP, bits=FP_BITS):
        self.w = w_cap
        self.esc_cap = esc_cap
        self.bits = bits
        self.phase = 0
        self.escapes = 0
        self.aborted = False
        self.escape_key = None
        self.fires_two = 0
        self.fires_ckpt = 0
        self.confirm_rejects = 0
        self.start(0, None, 0)

    # -- history -------------------------------------------------------------
    def start(self, t, key, state_hash):
        """Begin (or restart, after an escape) at step t with the given state hash."""
        self.t0 = t
        self.phase = state_hash % self.w
        self.ring = [None, None]        # full keys at t-1, t-2
        self.min_key = key              # running minimum of the current window
        self.prev_min = key             # minimum of the previous window
        self.ckpt_key = key             # key at the previous boundary
        self.escape_key = None

    # -- one step ------------------------------------------------------------
    def observe(self, t, key):
        """Feed K(t).  Returns None, 'two', 'ckpt', or 'abort'."""
        fired = None
        f = fp(key, self.bits)
        f2 = self.ring[1]
        if f2 is not None and f == fp(f2, self.bits):
            if key == f2:               # full-key confirmation
                fired = 'two'
            else:
                self.confirm_rejects += 1
        else:
            # keep the full keys so the confirmation above is possible
            pass
        self.ring[1] = self.ring[0]
        self.ring[0] = key
        if self.min_key is None or key < self.min_key:
            self.min_key = key
        if (t - self.t0) % self.w == 0 and t != self.t0:
            if key == self.ckpt_key or self.min_key == self.prev_min:
                if fired is None:
                    fired = 'ckpt'
                self.escape_key = self.min_key if self.min_key is not None else key
            self.prev_min = self.min_key
            self.ckpt_key = key
            self.min_key = key
        if fired is not None:
            self.escape_key = key          # CLASS-ONLY action: the current class key
            if fired == 'two':
                self.fires_two += 1
            else:
                self.fires_ckpt += 1
            if self.escapes >= self.esc_cap:
                self.aborted = True
                return 'abort'
            self.escapes += 1
        return fired
