#!/usr/bin/env python3
"""DP-free regression, written against the v10 implementation before repair."""
import unittest
from retention import run_walk

class SynthesisTest(unittest.TestCase):
    def test_dp_free_two_cycle_fires_charged_escape(self):
        def advance(st):
            return (1-st[0],1-st[0],0) if st[0]<2 else (3,0,1)
        out=run_walk(advance,lambda st:(2,2,0),(0,0,0),
                     lambda i:False,8,30,4,7)
        self.assertGreater(out['escapes'],0,'DP-free cycle must reach an escape hook')
        self.assertEqual(out['escape_R_adds'],out['escapes'])
        self.assertEqual(out['deletions'],0)

    def test_pool_overflow_aborts_without_any_deletion(self):
        from retention import FixedPool
        pool=FixedPool(2)
        pool.insert(1,(1,2,3));pool.insert(2,(4,5,6))
        before=dict(pool.records);pool.insert(3,(7,8,9))
        self.assertTrue(pool.full)
        self.assertEqual(pool.records,before)
        self.assertEqual(pool.deletions,0)
        self.assertEqual(pool.insert(1,(0,0,0)),(1,2,3))
        self.assertEqual(pool.records,before)

    def test_pre_escape_record_recovers_after_detector_reset(self):
        def advance(st):
            return (0,4,1) if st[0]==2 else (1-st[0],1-st[0],0)
        out=run_walk(advance,lambda st:(2,2,0),(0,0,0),
                     lambda i:True,8,30,4,7)
        self.assertEqual(out['recovered_scalar'],3)
        self.assertEqual(out['cross_escape_recovery'],1)
        self.assertEqual(out['escapes'],1)
        self.assertEqual(out['escape_retention_failures'],0)

    def test_allowance_exhaustion_fails_closed(self):
        out=run_walk(lambda st:(1-st[0],1-st[0],0),lambda st:st,
                     (0,0,0),lambda i:False,8,100,4,7)
        self.assertEqual(out['escapes'],4)
        self.assertEqual(out['fruitless_meetings'],5)
        self.assertEqual(out['allowance_failure'],1)
        self.assertEqual(out['deletions'],0)

    def test_escapes_unbounded_per_trail_but_charged(self):
        out=run_walk(lambda st:(0,st[1],st[2]),lambda st:st,
                     (0,0,0),lambda i:False,8,100,1<<40,7)
        self.assertEqual(out['ordinary_steps'],99)
        self.assertGreater(out['escapes'],4)
        self.assertEqual(out['escape_R_adds'],out['escapes'])
        self.assertEqual(out['allowance_failure'],0)
        self.assertEqual(out['escape_retention_failures'],0)

    def test_report_cap_stops_without_deleting(self):
        out=run_walk(lambda st:(st[0]+1,st[1]+1,0),lambda st:st,
                     (0,0,0),lambda i:True,8,30,4,101,report_cap=2)
        self.assertEqual(out['report_cap_failure'],1)
        self.assertEqual(out['reports'],2)
        self.assertEqual(out['pool_final'],2)
        self.assertEqual(out['deletions'],0)

    def test_exhaustive_period_four_all_DP_partitions(self):
        from itertools import product
        checked=blind=0
        for n in range(1,5):
            for dest in product(range(n),repeat=n):
                for start in range(n):
                    seen={};path=[];u=start
                    while u not in seen:
                        seen[u]=len(path);path.append(u);u=dest[u]
                    mu=seen[u];cycle=path[mu:];ell=len(cycle)
                    t=h=start
                    for tick in range(1,mu+ell+1):
                        t=dest[t];h=dest[dest[h]]
                        if t==h:break
                    self.assertEqual(t,h)
                    self.assertLessEqual(tick,mu+ell)
                    for mask in range(1<<n):
                        checked+=1
                        if not any(mask&(1<<i) for i in cycle):
                            blind+=1
                            self.assertEqual(t,h)
        self.assertGreater(checked,16000)
        self.assertGreater(blind,1000)

if __name__=='__main__':unittest.main()
