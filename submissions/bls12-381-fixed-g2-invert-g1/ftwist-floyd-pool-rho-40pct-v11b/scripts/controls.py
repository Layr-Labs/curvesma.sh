#!/usr/bin/env python3
"""Fixed-capacity append-only-until-full DP pool; escapes cannot erase it."""
class FixedPool:
    def __init__(self, capacity):
        if capacity <= 0:
            raise ValueError('capacity must be positive')
        self.capacity=capacity
        self.records={}
        self.deletions=0
        self.escapes=0
        self.insertions=0
        self.peak=0

    def insert(self, key, record):
        if key in self.records:
            return self.records[key]
        if len(self.records) >= self.capacity:
            # FIFO replacement is the ONLY deletion path. A duplicate lookup
            # never overwrites a record and an escape has no pool side effect.
            oldest=next(iter(self.records))
            del self.records[oldest]
            self.deletions+=1
        self.records[key]=record
        self.insertions+=1
        self.peak=max(self.peak,len(self.records))
        return None

    def escape_event(self):
        self.escapes+=1


def run_walk(advance, escape, start, is_dp, capacity, ordinary_cap,
             escape_cap, modulus, report_cap=1<<84):
    """One trail. Only a DP same-label hit may request an in-place jump.

    Distinct labels on any full-point match can recover, including a useful
    self-visit. This is optimistic relative to requiring a different worker.
    Pool records are immutable until FIFO replacement on a full pool.
    No hidden all-state detector, timeout jump, restart or DP deepening exists.
    """
    pool=FixedPool(capacity)
    state=start
    steps=epoch=reports=hits=cross=retention_failures=0
    recovered=None
    cap_failure=report_failure=0
    while steps < ordinary_cap:
        i,A,D=state
        if is_dp(i):
            if reports>=report_cap:
                report_failure=1
                break
            reports+=1
            old=pool.insert(i,(state,epoch,steps))
            if old is not None:
                hits+=1
                prior,prior_epoch,_=old
                den=(D-prior[2])%modulus
                num=(prior[1]-A)%modulus
                if den:
                    recovered=num*pow(den,-1,modulus)%modulus
                    cross=int(prior_epoch != epoch)
                    break
                assert num==0,'labels disagree on an exactly equal point'
                if pool.escapes >= escape_cap:
                    cap_failure=1
                    break
                snapshot=dict(pool.records)
                pool.escape_event()
                state=escape(state)
                epoch+=1
                retention_failures+=int(pool.records != snapshot)
                continue
        state=advance(state)
        steps+=1
    return dict(success=int(recovered is not None),recovered_scalar=recovered,
        ordinary_steps=steps,escapes=pool.escapes,escape_R_adds=pool.escapes,
        reports=reports,hits=hits,cross_escape_recovery=cross,
        escape_retention_failures=retention_failures,deletions=pool.deletions,
        insertions=pool.insertions,pool_peak=pool.peak,pool_final=len(pool.records),
        escape_cap_failure=cap_failure,report_cap_failure=report_failure)


def old_floyd(advance,escape,start,cap,r):
    """The inherited v9 detector semantics, two offsets and BOTH states reset."""
    tort=hare=start;esc=steps=0
    for _ in range(cap):
        tort=advance(tort);hare=advance(advance(hare));steps+=3
        if tort[0]!=hare[0]:continue
        den=(tort[2]-hare[2])%r;num=(hare[1]-tort[1])%r
        if den:return num*pow(den,-1,r)%r,esc,steps
        assert num==0
        tort=hare=escape(escape(tort));esc+=1
    return None,esc,steps

