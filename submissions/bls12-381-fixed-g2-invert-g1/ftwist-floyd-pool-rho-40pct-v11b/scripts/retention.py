#!/usr/bin/env python3
"""Floyd is a stall oracle; a zero-delete DP pool retains emitted credit.

Toy state: (full canonical point identifier, A, D). Production raw-state
alignment is equivalent, with full Frobenius confirmation at equality.
"""
class FixedPool:
    def __init__(self,capacity):
        if capacity<=0:raise ValueError('capacity must be positive')
        self.capacity=capacity;self.records={};self.insertions=0
        self.deletions=0;self.escapes=0;self.peak=0;self.full=False

    def insert(self,key,record):
        if key in self.records:return self.records[key]
        if len(self.records)>=self.capacity:
            self.full=True
            return None  # fail closed; never evict, overwrite, or continue
        self.records[key]=record;self.insertions+=1
        self.peak=max(self.peak,len(self.records))
        return None

    def escape_event(self):
        self.escapes+=1


def run_walk(advance,escape,start,is_dp,capacity,ordinary_cap,escape_allowance,
             modulus,report_cap=1<<84,observer=None):
    """Three evaluations/tick; report all emitted DPs, including hare hops.

    A fruitless confirmed Floyd meeting pays one in-place offset at the
    current tortoise. Only detector counters/hare initialization change.
    Scalar labels, current point history and every pool record survive.
    Per-trail escapes are UNBOUNDED but individually charged: each one costs
    one dense product plus the ESCAPE_ROWS instructions in the ledger, and
    the trail itself terminates at ordinary_cap. escape_allowance is the
    globally priced first-moment budget, not a per-trail cap; exhausting it
    fails closed, and the ledger prices that event by
    ESCAPE_MEAN_CEILING/ESCAPE_CAP. observer is read-only evidence
    instrumentation: it never changes the transition or escape decision.
    """
    pool=FixedPool(capacity);tort=hare=start
    steps=reports=hits=epoch=cross=meetings=fruitless=retention_failures=0
    recovered=None;allowance_failure=report_failure=pool_failure=0
    first_meeting_steps=0

    def report(state):
        nonlocal reports,hits,recovered,cross,report_failure,pool_failure
        if not is_dp(state[0]):return False
        if reports>=report_cap:
            report_failure=1;return True
        reports+=1
        old=pool.insert(state[0],(state,epoch,steps))
        if pool.full:
            pool_failure=1;return True
        if old is None:return False
        hits+=1;prior,prior_epoch,_=old
        den=(state[2]-prior[2])%modulus
        num=(prior[1]-state[1])%modulus
        if den:
            recovered=num*pow(den,-1,modulus)%modulus
            cross=int(prior_epoch!=epoch)
            return True
        assert num==0,'equal full points must have consistent zero labels'
        return False

    stopped=report(start)
    while not stopped and steps+3<=ordinary_cap:
        tort=advance(tort);steps+=1
        if report(tort):break
        hare=advance(hare);steps+=1
        if report(hare):break
        hare=advance(hare);steps+=1
        if report(hare):break
        if tort[0]!=hare[0]:continue
        meetings+=1
        if not first_meeting_steps:first_meeting_steps=steps
        den=(tort[2]-hare[2])%modulus
        num=(hare[1]-tort[1])%modulus
        if den:
            recovered=num*pow(den,-1,modulus)%modulus;break
        assert num==0,'confirmed meeting has inconsistent labels'
        fruitless+=1
        if observer is not None:observer(tort,hare,steps)
        # No per-trail escape cap: charge this offset and continue. Only the
        # globally priced allowance can terminate the trail early.
        if pool.escapes>=escape_allowance:
            allowance_failure=1;break
        # Counts plus immutable insertion-only semantics prove retention
        # without copying an exponential pool at each escape.
        before=(len(pool.records),pool.insertions,pool.deletions)
        pool.escape_event();tort=escape(tort);epoch+=1
        retention_failures+=int(before!=(len(pool.records),pool.insertions,pool.deletions))
        hare=tort  # detector reset, NOT trail/point/label/pool reset
        if report(tort):break
    return dict(success=int(recovered is not None),recovered_scalar=recovered,
        ordinary_steps=steps,escapes=pool.escapes,escape_R_adds=pool.escapes,
        reports=reports,hits=hits,cross_escape_recovery=cross,
        floyd_meetings=meetings,fruitless_meetings=fruitless,
        first_meeting_steps=first_meeting_steps,
        escape_retention_failures=retention_failures,deletions=pool.deletions,
        insertions=pool.insertions,pool_peak=pool.peak,pool_final=len(pool.records),
        allowance_failure=allowance_failure,report_cap_failure=report_failure,
        pool_cap_failure=pool_failure)
