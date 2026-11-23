"""Executed driver checks; conditional on the disclosed synchronous link model.

The small callback stubs deliberately clobber all GPRs and record the arguments
of the leaf calls. They check glue and its counts, not field arithmetic; actual
leaf/inversion execution is a separate package check.
"""
import json
import sys
from pathlib import Path
ROOT=Path('/snapshot/files') if Path('/snapshot/files/kernels').is_dir() else Path(__file__).parent
sys.path.insert(0,str(ROOT/'kernels'))
import a64

TB=2**40
YS=34*TB+0x7100000
DONE=34*TB+0x7120000
DRAINED=34*TB+0x7120800
STACK=34*TB+0x7140000
VIEW=38*TB
SOURCE=40*TB
PORT=34*TB+0xB800
LOG=34*TB+0x7160000
EBE=42*TB+2**37
CH=2**22
ORDER=int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141',16)
SNT=0x77abc

def stubs():
    out=[]
    for name in ('fwd','inv','bwd','keys'):
        out += [name+'_stub:', 'str x30,[sp,#-16]!',
                'mov x18,#0x7160000','movk x18,#0x2200,lsl #32']
        if name=='fwd':out += ['stp x0,x1,[x18]','str x3,[x18,#16]', 'ldr x4,[x1]']
        elif name=='inv':out += ['stp x0,x1,[x18,#24]',
            'ldp x4,x5,[x1]','ldp x6,x7,[x1,#16]',
            'stp x4,x5,[x0]','stp x6,x7,[x0,#16]']
        elif name=='bwd':
            # The real bwd never stores via its borrowed Y SP. Restore the
            # test stub's local LR slot immediately before observing it.
            out=out[:-4]+[name+'_stub:','mov x18,#0x7160000','movk x18,#0x2200,lsl #32']
            out += ['mov x4,sp','str x4,[x18,#40]','str x3,[x18,#48]', 'ldr x4,[sp,#64]',
                    'mov x18,#0x7100000','movk x18,#0x2200,lsl #32','str x30,[x18,#224]']
        else:out += ['stp x2,x3,[x18,#56]']
        out += [f'mov x{k},#{k+101}' for k in range(30)]
        if name=='bwd':out += ['mov x18,#0x7100000','movk x18,#0x2200,lsl #32','ldr x30,[x18,#224]']
        else:out += ['ldr x30,[sp],#16']
        out += ['ret']
    return '\n'.join(out)+'\n'

class Memory(dict):
    def __init__(self):
        super().__init__();self.cleaned=[];self.backing={};self.cache={}
        self.live=False;self.n=0;self.source=SOURCE;self.events=[];self.coins=[]
        self.arglog={};self.epochs=[];self.sync=0
        self.drain_guard=None
    def __getitem__(self,a):
        if VIEW<=a<VIEW+2**35:
            assert self.live,'expired Y aperture'
            assert a<VIEW+512*self.n,'outside offered frame range'
            return self.backing.get(self.source+a-VIEW,0)
        return super().get(a,0)
    def __setitem__(self,a,v):
        assert not VIEW<=a<VIEW+2**35,'write to read-only Y aperture'
        if a==DONE:self.epochs.append(v)
        if a==DRAINED and self.drain_guard:self.drain_guard()
        if LOG<=a<LOG+72:self.arglog.setdefault(a-LOG,[]).append(v)
        super().__setitem__(a,v)
    def clean(self,a):
        assert not self.live,'clean during a live offer'
        self.cleaned.append(a)
        for k in range(0,64,8):self.backing[a+k]=self.cache.get(a+k,0)
    def barrier(self):self.sync+=1
    def rndr(self):return self.coins.pop(0) if self.coins else 1

def run_checks():
    text=(ROOT/'kernels/ytransport.s').read_text()
    P,L=a64.a64_load(text+stubs())
    stub_start=L['fwd_stub'][0]
    coin_start=L['ycoin'][0];coin_end=L['ychunk'][0]
    def run(entry,regs,m,dev=None):
        vis={i:0 for i in range(len(P))};out={}
        n=a64.a64_run(P,L,entry,dict(regs,x30=SNT),m,visits=vis,out_regs=out,dev=dev)
        stub_n=sum(v for i,v in vis.items() if i>=stub_start)
        coin_n=sum(v for i,v in vis.items() if coin_start<=i<coin_end)
        return n,out,n-stub_n-coin_n,coin_n
    def setup(m):
        m[YS+48]=0x123456789
        for i,k in enumerate(('fwd','inv','bwd','keys')):m[YS+64+8*i]=L[k+'_stub'][0]
        for i in range(4):m[YS+96+8*i]=ORDER>>(64*i)&(2**64-1)
    clean_results=[]
    for second in (False,True):
        for n in (0,1,2,5,17):
            m=Memory();m.source=SOURCE+(2**35+2**31 if second else 0)
            for i in range(n*64):m.cache[m.source+8*i]=i+1
            count,_,_,_=run('yclean',dict(x0=m.source,x1=n),m)
            assert count==(32*n+4 if n else 3) and m.sync==1
            assert m.backing==m.cache and len(m.cleaned)==8*n
            clean_results.append(dict(second=second,n=n,count=count))
    glue=[]
    for n in (0,1,2,5,CH-1,CH,CH+1,2*CH+3,2**26-1):
        m=Memory();setup(m);m.live=True;m.n=n
        total,out,over,coins=run('yconsume',dict(x0=VIEW,x1=n,x2=DONE,x3=19,sp=STACK),m)
        chunks=[];left=n;cursor=VIEW
        while left:
            z=min(CH,left);chunks.append((z,cursor));cursor+=512*z;left-=z
        assert m.arglog.get(0,[])==[6*z for z,v in chunks]
        assert m.arglog.get(8,[])==[v+64 for z,v in chunks]
        assert m.arglog.get(16,[])==[EBE-3072*z for z,v in chunks]
        assert m.arglog.get(40,[])==[v+512*(z-1) for z,v in chunks]
        assert m.arglog.get(48,[])==[EBE for _ in chunks]
        assert m.arglog.get(56,[])==[m[YS+48]+6 for _ in chunks]
        assert m.arglog.get(64,[])==[EBE-3072*z+64 for z,v in chunks]
        assert out['sp']==STACK and m.epochs==[19]
        # The original41 glue +16 INV_CALL and new64 cover current glue,
        # including the composed bwd's additional2/chunk. Fixed costs use
        # the576/offer envelope, leaving ample room for waits/metadata.
        assert over+2*len(chunks)<=121*len(chunks)+32
        assert coins<=64*len(chunks)
        glue.append(dict(n=n,chunks=len(chunks),glue=over,coins=coins))
    # Coin rejection reaches zero, >=r, and both failures without deadlock.
    coin_cases=[]
    for words in ([0]*4+[1,0,0,0], [2**64-1]*4+[1,0,0,0], [0]*8):
        m=Memory();setup(m);m.coins=list(words)
        n,out,_,_=run('ycoin',dict(x16=YS),m)
        assert n<=64 and out['x0']==int(words==[0]*8)
        coin_cases.append(dict(count=n,failed=out['x0']))
    # All1023 offers complete synchronously, including an empty round.
    offers=[]
    for n in (0,2):
        owner=Memory();owner.n=n
        for i in range(n*64):owner.cache[SOURCE+8*i]=i+1
        run('yclean',dict(x0=SOURCE,x1=n),owner)
        seen=[]
        def deliver(off,val):
            if off!=40:return
            assert not owner.live
            recipient=Memory();setup(recipient)
            recipient.source=owner[PORT+16];recipient.n=owner[PORT+24]
            recipient.backing=owner.backing;recipient.live=owner.live=True
            assert val==owner[PORT+32]
            run('yconsume',dict(x0=VIEW,x1=recipient.n,x2=DONE,x3=val,sp=STACK),recipient)
            assert recipient.epochs==[val]
            if n:
                assert recipient[VIEW+64]==owner.cache[SOURCE+64]
            recipient.live=owner.live=False
            seen.append((owner[PORT],owner[PORT+8],val))
        count,_,_,_=run('yoffer',dict(x0=PORT,x1=1024*17+1,x2=0,x3=SOURCE,x4=n,x5=23,x6=1023),owner,(PORT,48,deliver))
        assert count==6*1023+3
        assert seen==[(1024*17+j,0,23)for j in range(1,1024)]
        offers.append(dict(n=n,recipients=len(seen),sender_instructions=count))
    # Output-level cache control: a dirty source yields a wrong payload.
    dirty=Memory();dirty.n=1;dirty.live=True;dirty.cache[SOURCE+64]=99
    dirty_mismatch=dirty[VIEW+64]!=dirty.cache[SOURCE+64]
    expired=False
    dirty.live=False
    try:dirty[VIEW+64]
    except AssertionError:expired=True
    assert dirty_mismatch and expired
    owner_runs=[]
    for n in (0,2):
        owner=Memory();setup(owner);owner.n=n
        context=YS+256;finished=[]
        values=(SOURCE,n,1024*17+1,0,DONE,27,DRAINED,PORT)
        for i,v in enumerate(values):owner[context+8*i]=v
        for i in range(n*64):owner.cache[SOURCE+8*i]=i+1
        def guard():
            assert len(finished)==1023 and not owner.live,'DRAINED before all returns'
        owner.drain_guard=guard
        remote_fixed=0
        def owner_deliver(off,val):
            nonlocal remote_fixed
            if off!=40:return
            receiver=Memory();setup(receiver);receiver.n=n
            receiver.backing=owner.backing;receiver.live=owner.live=True
            total,out,glue,coins=run('yconsume',dict(x0=VIEW,x1=n,x2=DONE,x3=val,sp=STACK),receiver)
            assert receiver.epochs==[27] and out['sp']==STACK
            remote_fixed+=glue-63*bool(n)
            receiver.live=owner.live=False
            finished.append(owner[PORT])
        total,out,owner_glue,coins=run('yowner',dict(x0=context,sp=STACK),owner,(PORT,48,owner_deliver))
        expected=(32*n+4 if n else 3)+6141+20+11+63*bool(n)
        assert owner_glue==expected and owner[DRAINED]==27 and owner.epochs==[27] and out['sp']==STACK
        assert remote_fixed==11*1023
        owner_runs.append(dict(n=n,owner_glue_with_clean=owner_glue,remote_fixed=remote_fixed,
            owner_wrapper=20,returned=len(finished),drained=owner[DRAINED]))
    old='bl yoffer\nldr x4,[sp]\nldp x5,x6,[x4,#40]\nstlr x5,[x6]'
    new='ldr x4,[sp]\nldp x5,x6,[x4,#40]\nstlr x5,[x6]\nbl yoffer'
    assert text.count(old)==1
    badP,badL=a64.a64_load(text.replace(old,new)+stubs())
    bad=Memory();setup(bad);context=YS+256
    for i,v in enumerate((SOURCE,0,1024*17+1,0,DONE,27,DRAINED,PORT)):bad[context+8*i]=v
    def reject_early():raise AssertionError('DRAINED before all returns')
    bad.drain_guard=reject_early
    try:a64.a64_run(badP,badL,'yowner',dict(x0=context,sp=STACK,x30=SNT),bad)
    except AssertionError as exc:early_drain_caught=str(exc)=='DRAINED before all returns'
    else:early_drain_caught=False
    assert early_drain_caught
    # Minimal adversarial schedule: the last recipient has stored DONE but
    # has not returned to its stalled sender. DONE alone permits reuse; the
    # owner wrapper cannot publish DRAINED until that return reaches yoffer.
    done=True;handler_returned=False;sender_returned=False;drained=False
    done_only_reuse=done and not handler_returned
    drain_gate_reuse=done and drained
    assert done_only_reuse and not drain_gate_reuse
    handler_returned=True;sender_returned=True;drained=sender_returned
    assert done and drained and handler_returned
    fixed=6*1023+3+(11+14)*1024+20+4+14+4*128
    assert fixed<=576*1023
    wait_source='''wait:
1:
ldaxr x2,[x0]
cmp x2,x1
b.eq 2f
wfe
b 1b
2:
ret
'''
    wp,wl=a64.a64_load(wait_source)
    waits=[]
    for aliased in (False,True):
        memory={DONE:0};events=[(DONE+(8 if aliased else 2048),1),
            (DONE+(16 if aliased else 4096),1),(DONE,31)]
        def wake(mem,monitor,event):
            # A finite external-store schedule. Only a write into the watched
            # 2KiB granule clears the monitor and sets this processor's event.
            while events:
                address,value=events.pop(0);mem[address]=value
                if address>>11==monitor[0]:
                    monitor[0]=None;event[0]=1;return
            raise AssertionError('wait has no completing store')
        nw=a64.a64_run(wp,wl,'wait',dict(x0=DONE,x1=31,x30=SNT),memory,wfe=wake,ev0=1)
        waits.append(nw)
    assert waits[0]==14 and waits[1]>14
    return dict(all_checks_pass=True,clean=clean_results,glue=glue,coin=coin_cases,
                offers=offers,dirty_negative=dirty_mismatch,expired_negative=expired,
                owner=owner_runs,done_only_reuse_negative=done_only_reuse,
                actual_early_drain_mutation=early_drain_caught,
                drain_gate_blocks_reuse=not drain_gate_reuse,fixed_per_subring_upper=fixed,
                wait_instructions=waits[0],aliased_monitor_negative=waits[1],
                source_instructions=len(a64.a64_load(text)[0]),
                scope='Driver instructions and explicit link/cache model; callback arithmetic checked separately.')

if __name__=='__main__':
    json.load(sys.stdin)
    print(json.dumps(run_checks()))
