#!/usr/bin/env python3
"""sigma_set.py -- EXHAUSTIVE answer to the question the v14 dispatch posed:

    does a nonempty subset S of the twelve column exponents Z/12 exist that is
    (a) closed under the routing action k -> k+1 that sigma induces on the
        column index, and (b) free of inverse pairs {k, k+6}?

Why it is the right question.  The multiplier of a multiply step is
m = TBL[t][kappa(y)] with TBL[t][k] = sigma^{-k}(M_t).  Merge-safety requires
m(sigma^j y) = sigma^j m(y) for every j (that is what makes the successor CLASS a
function of the class), and kappa(sigma^j y) = kappa(y) - j.  If the table offers
columns only from a set S, the routing must map every residue into S, so S has to
be closed under k -> k+1.  If S is closed under +1 it is closed under the whole
transitive Z/12 action, so it is empty or everything.

Consequence, proved twice below (once by enumeration over all 2^12 subsets, once by
the transitivity argument, which the enumeration also confirms): the largest
inversion-free closed set is EMPTY.  A six-column lower-half table (v12b's, which
the v12b reviewer praised as 'inverse-free') is NOT closed: on it kappa and
kappa+6 select the same column, so m(sigma^6 y) = m(y) while sigma^6 m(y) =
m(y)^{-1}, forcing m^2 = 1 for EVERY table entry, and a nonidentity GT element has
no order-2 point when r is odd.  That is exactly the equivariance failure v12b's
own verdict records, now derived rather than observed.  Hence v14 ships the
routing GUARD (classkey.Guard), not a column-subset table.
"""
import itertools
import json

N = 12
PAIR = {k: (k + 6) % 12 for k in range(N)}


def closed(S):
    """closed under the routing action k -> k+1 (mod 12)?"""
    return all(((k + 1) % 12) in S for k in S)


def inversion_free(S):
    """no k with both k and k+6 present?"""
    return not any((k in S) and (PAIR[k] in S) for k in S)


def enumerate_subsets():
    closed_sets, both, max_if_closed, hist = [], [], None, {}
    for mask in range(1 << N):
        S = frozenset(k for k in range(N) if mask >> k & 1)
        c, f = closed(S), inversion_free(S)
        if c:
            closed_sets.append(sorted(S))
        if c and f:
            both.append(sorted(S))
        if c and f and (max_if_closed is None or len(S) > len(max_if_closed)):
            max_if_closed = sorted(S)
        hist[(c, f)] = hist.get((c, f), 0) + 1
    return dict(subsets=1 << N, closed=[sorted(s) for s in closed_sets],
                closed_and_inversion_free=both,
                largest_closed_inversion_free=(max_if_closed if max_if_closed else []),
                histogram={'closed_and_free': hist.get((True, True), 0),
                           'closed_not_free': hist.get((True, False), 0),
                           'not_closed_free': hist.get((False, True), 0),
                           'not_closed_not_free': hist.get((False, False), 0)})


def six_column_refutation():
    """The lower-half table is inverse-free but not closed, and closure fails by a
    forced involution: the two exponents k and k+6 collapse to one column."""
    S = list(range(6))
    return dict(set=S, closed=closed(set(S)), inversion_free=inversion_free(set(S)),
                collapse_pairs=[(k, (k + 6) % 12) for k in S],
                forced_relation='m(sigma^6 y) = m(y) and sigma^6 m(y) = m(y)^-1 '
                                '=> m^2 = 1 for every entry; r is odd so the only '
                                'such GT element is the identity')


def transitivity_argument():
    """One line: +1 mod 12 generates Z/12, so a nonempty closed subset is Z/12 and
    Z/12 contains {0,6}."""
    orbit = set()
    k = 0
    for _ in range(N):
        orbit.add(k)
        k = (k + 1) % 12
    return dict(orbit_of_0=sorted(orbit), is_transitive=(len(orbit) == N),
                full_set_has_inverse_pair=(0 in orbit and 6 in orbit))


if __name__ == '__main__':
    out = dict(enumeration=enumerate_subsets(), lower_half=six_column_refutation(),
               transitivity=transitivity_argument(),
               conclusion='no nonempty sigma-stable inversion-free column set exists; '
                          'the routing guard is the only remaining structural route')
    json.dump(out, open('evidence/sigma-set.json', 'w'), indent=1, sort_keys=True)
    e = out['enumeration']
    print('subsets scanned %d' % e['subsets'])
    print('closed sets            %s' % e['closed'])
    print('closed AND inv-free    %s' % e['closed_and_inversion_free'])
    print('largest closed inv-free %s' % e['largest_closed_inversion_free'])
    print('lower half: closed %s  inversion_free %s'
          % (out['lower_half']['closed'], out['lower_half']['inversion_free']))
    print('transitive %s' % out['transitivity']['is_transitive'])
