#!/usr/bin/env python3
"""End-to-end toy: the production walk on toy towers of exactly the production shape.

Scale changes four numbers only: the field modulus, the class count N, the
distinguished density (2^-7 instead of 2^-80, because a toy trail must be shorter than
sqrt(N)) and the escape budget.  Everything else is the production algorithm:

  * the class rule: kappa(y) = the minimum of the twelve conjugate keys of the held
    w-coefficient, together with the position k of that minimum, which is the Frobenius
    power that carries y to the canonical representative,
  * the EQUIVARIANT step rule from classwalk.py: the table entry is addressed at
    (routing twist + k) mod 12, i.e. sigma^(tw+k)(m0).  This is the v6 repair of the v5
    defect: with the index frozen at a per-trail value the step was not a class function
    and the walk effectively ran on r = 12N elements,
  * the step mixture, thresholded on the routing word by the same f the budget uses,
  * the detector from detector.py (imported, not reimplemented): two-cycle rule on
    fp(K) vs fp(K-2) with full-key confirmation, checkpoint rule on the windowed
    minimum, escape = double traversal-min with history clear and phase re-derivation,
  * the metered trail cutoff and the fixed report cap with a refusal flag,
  * replay WITH the detector re-executed (the detector is history dependent, so a
    replay that skipped it would not reproduce the report), the Frobenius alignment,
    the exponent solve and the verification against z.

Two probes are recorded besides the end-to-end run:

  * `equivariance_probe` -- on real tower elements, take y and its twelve Frobenius
    images and check that a step taken from each of them lands in ONE class.  This is
    the property the v5 step lacked (59 per cent disagreement) and the one the claim
    needs; it is measured, not asserted.
  * `cycle_regression` (2..8) and `cycle_adversarial` ({2,4,8,16,32,64,67}) -- synthetic
    closed cycles at every checkpoint phase, so the detector's fire rule is exercised
    in the regime where a fixed-phase escape re-enters the same trap.
"""
import json, random, sys, time
import tower
import routing
import classwalk
import detector as D
import constants as K

TOY_DP_BITS = 7           # toy distinguished density 2^-7 (production: 2^-80)
MAX_SEG = 4096            # hard segment cap (never reached: the metered cutoff binds)
R_TABLE = 16
TOY_X = -56               # BLS12 parameter, same parameterisation as BLS12-381


def primes(n):
    s = [True]*(n+1); s[0] = s[1] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]:
            for j in range(i*i, n+1, i):
                s[j] = False
    return [i for i in range(n+1) if s[i]]


SMALL = primes(1 << 16)
is_prime = tower.is_prime


def find_toy():
    """A BLS12 tower of production shape: r = x^4-x^2+1, p = (x^6-2x^5+2x^3+x+1)/3."""
    x = TOY_X
    r = x**4 - x**2 + 1
    p = (x**6 - 2*x**5 + 2*x**3 + x + 1)//3
    assert is_prime(p) and is_prime(r), (p, r)
    assert p % 12 == 7, p % 12
    assert (r - 1) % 12 == 0
    assert all(pow(p, k, r) != 1 for k in (1, 2, 3, 4, 6)) and pow(p, 12, r) == 1
    return p, r, (r-1)//12


def build_table(T, r, g, z, rng, R=R_TABLE):
    """R base multipliers, each with its twelve Frobenius twists and the pre-twisted
    exponents, so a step needs no exponent conversion."""
    muls = []
    for _ in range(R):
        u = rng.randrange(1, r); v = rng.randrange(1, r)
        el = T.mul(T.pow(g, u), T.pow(z, v))
        twists = []
        for i in range(12):
            twists.append((el, (u*pow(T.p, i, r)) % r, (v*pow(T.p, i, r)) % r))
            el = T.frob(el)
        muls.append(twists)
    return muls


def key_and_twist(T, y):
    """Class key and held Frobenius offset: the argmin over the twelve conjugate keys.

    Both come out of one pass over the twelve keys (tower.kappa), which is the fold row
    already charged in ledger.py, so the corrected rule costs no extra field arithmetic.
    """
    return T.kappa(y)


def step(T, y, muls, r, a, b, key, held):
    """One production step from the held element y, given its (key, held twist).

    Returns (y', a', b', key', held').  Exactly ONE kappa fold is paid here -- the fold
    of the new state -- and its key is carried to the next step, which is the production
    accounting (the fold of a state is simultaneously its routing key, its distinguished
    test and its next twist).
    """
    tbl, tw, dp, st = routing.fields(key)
    if st < K.SQR_THRESHOLD:
        ny = T.cyclo_sqr_fast(y); na = (2*a) % r; nb = (2*b) % r
    else:
        el, uu, vv = muls[tbl % R_TABLE][classwalk.twist_index(tw, held)]
        ny = T.mul(y, el); na = (a + uu) % r; nb = (b + vv) % r
    nkey, nheld = key_and_twist(T, ny)
    return ny, na, nb, nkey, nheld


def escape_step(T, y, muls, r, a, b, escape_key, key, held):
    """Double traversal-min: two table steps at the multiplier indexed by escape_key.

    The index is (escape twist - held twist), exactly as in step(), so the escape is a
    class map: two trails sitting at the common class perform the identical action.
    Returns (y', a', b', key', held') with one fold for the final state.
    """
    tbl, tw, _, _ = routing.fields(escape_key)
    el, uu, vv = muls[tbl % R_TABLE][classwalk.twist_index(tw, held)]
    y = T.mul(y, el); a = (a + uu) % r; b = (b + vv) % r
    y = T.mul(y, el); a = (a + uu) % r; b = (b + vv) % r
    k, h = key_and_twist(T, y)
    return y, a, b, k, h


def trail_iter(T, muls, r, trail, limit, dp_bits, det=None):
    """Yield ('step', t, y, a, b, key) for each step, applying the detector; stops at
    the first report or at the metered cutoff."""
    y, a, b = trail
    if det is None:
        det = D.Detector()
        k0, _ = key_and_twist(T, y)
        det.start(0, k0, routing.mix(k0))
    k, held = key_and_twist(T, y)
    for t in range(1, limit+1):
        y, a, b, k, held = step(T, y, muls, r, a, b, k, held)
        ev = det.observe(t, k)
        if ev == 'abort':
            yield ('abort', t, y, a, b, k)
            return
        if ev is not None:
            y, a, b, k, held = escape_step(T, y, muls, r, a, b, det.escape_key, k, held)
            det.start(t, k, routing.mix(k))
        yield ('step', t, y, a, b, k)


def run_group(T, r, g, z, muls, rng, cap, segments, dp_bits, w_trail):
    """Independent trails, each ending at one report; the ordered store detects a repeat."""
    store = {}
    loss = False
    trails = []
    reports = 0
    aborts = 0
    escapes = 0
    for seg in range(segments):
        a0 = rng.randrange(1, r); b0 = rng.randrange(1, r)
        y0 = T.mul(T.pow(g, a0), T.pow(z, b0))
        trails.append((y0, a0, b0))
        k0, _ = key_and_twist(T, y0)
        det = D.Detector()
        det.start(0, k0, routing.mix(k0))
        hit = None
        for ev in trail_iter(T, muls, r, trails[-1], w_trail, dp_bits, det):
            if ev[0] == 'abort':
                aborts += 1
                break
            _, t, y, a, b, k = ev
            if (routing.fields(k)[2] & ((1 << dp_bits)-1)) != 0:
                continue
            reports += 1
            if k in store:
                hit = (store[k], (seg, t))
                break
            if len(store) >= cap:
                loss = True
                break
            store[k] = (seg, t)
        escapes += det.escapes
        if hit is not None:
            return dict(hit=True, pair=hit, loss=loss, reports=reports, trails=trails,
                        aborts=aborts, escapes=escapes)
        if loss:
            break
    return dict(hit=False, loss=loss, reports=reports, trails=trails, size=len(store),
                aborts=aborts, escapes=escapes)


def replay(T, muls, r, trail, steps, dp_bits):
    """Re-execute the trail with the detector: the detector is history dependent, so a
    replay that skipped it could not reproduce the report."""
    k0, _ = key_and_twist(T, trail[0])
    det = D.Detector()
    det.start(0, k0, routing.mix(k0))
    last = None
    for ev in trail_iter(T, muls, r, trail, steps, dp_bits, det):
        if ev[1] == steps:
            last = ev
            break
    return last[2], last[3], last[4], last[5]


def group(T, r, g, z, muls, rng, secret, cap, segments, dp_bits, w_trail):
    out = run_group(T, r, g, z, muls, rng, cap, segments, dp_bits, w_trail)
    if not out['hit']:
        return dict(recovered=False, hit=False, loss=out['loss'],
                    aborts=out['aborts'], escapes=out['escapes'])
    (t1, s1), (t2, s2) = out['pair']
    y1, a1, b1, k1 = replay(T, muls, r, out['trails'][t1], s1, dp_bits)
    y2, a2, b2, k2 = replay(T, muls, r, out['trails'][t2], s2, dp_bits)
    if k1 != k2:
        return dict(recovered=False, hit=True, replay_mismatch=True,
                    aborts=out['aborts'], escapes=out['escapes'])
    j = None
    for kk in range(12):
        if y1 == T.pow(y2, pow(T.p, kk, r)):
            j = kk
            break
    if j is None:
        return dict(recovered=False, hit=True, aligned=False,
                    aborts=out['aborts'], escapes=out['escapes'])
    pj = pow(T.p, j, r)
    num = (a2*pj - a1) % r
    den = (b1 - b2*pj) % r
    if den == 0:
        return dict(recovered=False, hit=True, aligned=True, singular=True,
                    aborts=out['aborts'], escapes=out['escapes'])
    x = num*pow(den, r-2, r) % r
    ok = (T.pow(g, x) == z)
    return dict(recovered=True, hit=True, aligned=True, verified=ok,
                secret_match=(x == secret), aborts=out['aborts'], escapes=out['escapes'])


def equivariance_probe(T, r, g, z, rng, samples=60):
    """Measure the class-function property on real tower elements.

    For each sample y, walk the twelve conjugates sigma^k(y) and take one step from each
    with the same multiplier table; all twelve next elements must lie in ONE class.  The
    v5 rule failed this 59 per cent of the time; the v6 rule must fail it zero times.
    """
    R = 16
    muls = build_table(T, r, g, z, rng, R)
    probes = 0
    agree = 0
    for _ in range(samples):
        a = rng.randrange(1, r); b = rng.randrange(1, r)
        y = T.mul(T.pow(g, a), T.pow(z, b))
        y0 = y
        keys = []
        kk, hh = key_and_twist(T, y)
        for k in range(12):
            nxt, _, _, nkey, _ = step(T, y, muls, r, a, b, kk, hh)
            keys.append(nkey)                          # key of the NEXT element
            y = T.frob(y)
            kk, hh = key_and_twist(T, y)
        probes += 1
        if len(set(keys)) == 1:
            agree += 1
    return dict(samples=probes, all_twelve_agree=agree,
                fraction=agree/probes if probes else 0.0)


def cycle_regression(lengths, w=64):
    """Adversarial short-cycle regression, exhaustive over cycle length and phase.

    A synthetic cycle of length L with distinct keys is walked for 3*w steps at every
    phase in [0, w).  The checkpoint rule must fire for every (L, phase) pair for
    L <= w -- this is the regression for the documented fixed-window trap, where a
    4-cycle A->B->C->D->A with the boundary at B never escapes.  Lengths above w are
    outside the rule's guarantee (a window of w never contains the whole cycle) and are
    recorded rather than required.
    """
    out = dict(lengths=lengths, window=w, phases=w, all_phases_fire=True,
               trapped=0, fired=0, phases_by_length={}, within_window=list(
                   x for x in lengths if x <= w), beyond_window=list(
                   x for x in lengths if x > w))
    for L in lengths:
        keys = [1000 + 7*i for i in range(L)]      # distinct, increasing
        fired_here = 0
        for phase in range(w):
            det = D.Detector()
            det.start(0, keys[0], phase)           # phase = the per-trail offset
            t = 0
            fired = 0
            for _ in range(3*w):
                t += 1
                k = keys[(t % L)]
                ev = det.observe(t, k)
                if ev is not None and ev != 'abort':
                    fired += 1
                    break
            if fired:
                fired_here += 1
        out['phases_by_length'][str(L)] = fired_here
        out['fired'] += fired_here
        if L <= w and fired_here < w:
            out['all_phases_fire'] = False
            out['trapped'] += w - fired_here
    return out


def main():
    t0 = time.time()
    p, r, N = find_toy()
    T = tower.Tower(p)
    rng = random.Random(424242)
    e = (p**12 - 1)//r
    g = None
    while g is None:
        x = tuple(tuple((rng.randrange(p), rng.randrange(p)) for _ in range(3)) for _ in range(2))
        h = T.pow(x, e)
        if h != T.one12() and T.pow(h, r) == T.one12():
            g = h
    S = K.step_budget(n=N)
    theta_bits = TOY_DP_BITS
    # metered cutoff: 16 DP lengths, and never past the rho scale of this field, so a
    # trail cannot reach a cycle of the class map (production: 16*2^80 << 2^125.6)
    w_trail = min(K.TRAIL_METERS << theta_bits, int((3.14159*N/2) ** 0.5))
    cap = (11*S) // (10 << theta_bits) + K.CAP_SLACK     # Chernoff report cap
    groups, segments = 40, 80
    res = dict(groups=0, hits=0, recovered=0, verified=0, secret_match=0, cap_losses=0,
               aligned_fail=0, singular=0, replay_mismatch=0, aborts=0, escapes=0,
               fires_two=0, fires_ckpt=0, confirm_rejects=0,
               p=p, r=r, N=N, dp_bits=theta_bits, segments=segments, store_cap=cap,
               w_trail=w_trail, toy_x=TOY_X, S=int(S), sqr_threshold=K.SQR_THRESHOLD,
               trail_meters=K.TRAIL_METERS, rule="v6-equivariant")
    for gi in range(groups):
        secret = rng.randrange(1, r)
        z = T.pow(g, secret)
        muls = build_table(T, r, g, z, rng)
        o = group(T, r, g, z, muls, rng, secret, cap, segments, theta_bits, w_trail)
        res['groups'] += 1
        for k in ('aborts', 'escapes'):
            res[k] += o.get(k, 0)
        if o.get('hit'):
            res['hits'] += 1
        if o.get('loss'):
            res['cap_losses'] += 1
        if o.get('recovered'):
            res['recovered'] += 1
            if o.get('verified'):
                res['verified'] += 1
            if o.get('secret_match'):
                res['secret_match'] += 1
        if o.get('hit') and o.get('aligned') is False:
            res['aligned_fail'] += 1
        if o.get('replay_mismatch'):
            res['replay_mismatch'] += 1
        if o.get('singular'):
            res['singular'] += 1
    res['equivariance_probe'] = equivariance_probe(T, r, g, T.pow(g, 12345), rng)
    # fixed-cap refusal: a store that is too small must refuse and set the loss flag
    secret = rng.randrange(1, r)
    z = T.pow(g, secret)
    muls = build_table(T, r, g, z, rng)
    small = run_group(T, r, g, z, muls, rng, 1, segments, theta_bits, w_trail)
    res['cap_probe_refused'] = bool(small['loss'])
    res['cap_probe_hit'] = bool(small['hit'])
    res['detector_cycle_test'] = cycle_regression(list(range(2, 9)))
    res['detector_cycle_adversarial'] = cycle_regression([2, 4, 8, 16, 32, 64, 67])
    res['seconds'] = round(time.time()-t0, 1)
    res['all_verified'] = (res['recovered'] == res['verified'] == res['secret_match']
                           and res['recovered'] > 0 and res['cap_probe_refused']
                           and res['detector_cycle_test']['all_phases_fire']
                           and res['detector_cycle_adversarial']['all_phases_fire']
                           and res['equivariance_probe']['all_twelve_agree']
                           == res['equivariance_probe']['samples']
                           and res['replay_mismatch'] == 0)
    json.dump(res, open("toy-out.json", "w"), indent=1)
    print("toy p=%d r=%d N=%d dp_bits=%d segments=%d w_trail=%d cap=%d"
          % (p, r, N, theta_bits, segments, w_trail, cap))
    print("groups %d hits %d recovered %d verified %d secret_match %d cap_loss %d "
          "aligned_fail %d replay_mismatch %d"
          % (res['groups'], res['hits'], res['recovered'], res['verified'],
             res['secret_match'], res['cap_losses'], res['aligned_fail'],
             res['replay_mismatch']))
    print("escapes %d aborts %d" % (res['escapes'], res['aborts']))
    print("equivariance: %d/%d samples put all twelve images in one next class"
          % (res['equivariance_probe']['all_twelve_agree'],
             res['equivariance_probe']['samples']))
    print("cap probe: refused %s  cycle 2..8 all phases fire: %s  adversarial %s"
          % (res['cap_probe_refused'], res['detector_cycle_test']['all_phases_fire'],
             res['detector_cycle_adversarial']['phases_by_length']))
    print("TOY", "PASS" if res['all_verified'] else "FAIL", "%.1fs" % res['seconds'])
    return 0 if res['all_verified'] else 1


if __name__ == "__main__":
    sys.exit(main())
