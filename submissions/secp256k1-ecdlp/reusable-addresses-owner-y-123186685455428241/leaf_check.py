"""Compact current-leaf check and prepared-memory fixture.

Import is inert. Call only in the credential-free verifier with an initialized
reference experiment module E and certificate K. No lcert/import time is added.
leaf_check uses independently derived small table/record values. run_cold
executes the current table builder and Y producer separately in the same suite.
"""
def leaf_check(E,K,leaf_source,pointer_init_source,route_pair=None,n=4):
    import time
    started=time.time()
    rp,rl=route_pair or E.kload('route')
    orig=(E.TLC,E.popk,E.offs7,E.arun,E.CE)
    OldMemory,old_popk,old_offs,arun,_=orig
    E.CE=K
    view=38*(1<<40)
    pair=(view,view+(1<<35))
    add_dev=pair not in E.a64.DEV
    if add_dev:E.a64.DEV.append(pair)
    init=pointer_init_source.replace('mov x0,#0x2a1d00000000\n','').replace('mov x3,#4194304\n','')
    ip,il=E.a64.a64_load(init)
    I,L=E.a64.a64_load(leaf_source)
    saved={};memories=[]
    class ReadOnlyViewWrite(AssertionError):pass
    class Guarded(OldMemory):
        def __init__(self,*args):
            super().__init__(*args);self.protect=False
            self.phase=None;self.view=(view,view+512*n)
            self.view_reads=set();self.view_writes=[]
            self.coins=iter(())
            memories.append(self)
        def __getitem__(self,a):
            if self.protect and self.view[0]<=a<self.view[1]:self.view_reads.add(a)
            return super().__getitem__(a)
        def __setitem__(self,a,v):
            if self.protect and self.view[0]-64<=a<self.view[1]:
                self.view_writes.append(a);raise ReadOnlyViewWrite(a)
            return super().__setitem__(a,v)
        def rndr(self):return next(self.coins)
        def sync(self,*args):pass
    def linked_offs(v,i):
        out=old_offs(v,i);out[0]+=64;return out
    def linked_popk(mem,xs,b):
        old_popk(mem,xs,b)
        for i,v in enumerate((b,pow(b,3,E.p),3*b*b%E.p)):
            for a in old_offs(v,i):mem[a+120]=a;mem[a+248]=a+128
    def run(prog,labs,entry,regs,N,*a,**kw):
        base=N
        while hasattr(base,'m'):base=base.m
        hot=isinstance(base,Guarded) and entry in ('fwd','bwd','keys')
        if hot and entry=='fwd':
            sb=regs['x3'];ny=regs['x0']//6
            saved['initializer']=arun(ip,il,'card44_init',dict(x0=sb,x3=ny),base)
            addresses=[]
            for g in range(ny):
                e=sb+3072*g
                addresses += [e+j*512+176 for j in range(6)]+[e+128]
                for j in (0,3):addresses += [e+(j+1)*512+128,e+j*512+144,e+j*512+160]
            saved['pointers']={a:base[a] for p in addresses for a in (p,p+8)}
            saved['initial']=dict(base);saved['fwd_regs']=dict(regs);saved['sb']=sb
        old_protect=base.protect if isinstance(base,Guarded) else False
        if hot:base.protect=True
        try:return arun(prog,labs,entry,regs,N,*a,**kw)
        finally:
            if hot:base.protect=old_protect
    try:
        E.TLC=Guarded;E.popk=linked_popk;E.offs7=linked_offs;E.arun=run
        rng=E.Rd(20261123)
        xs=E.mbx(rng)
        ys=[E.PROD.mul(rng.rr(1,E.r),E.GP) for _ in range(n)]
        rho=7
        d=E.lfb(I,L,xs,ys,6*n,None,K,rprog=rp,rlab=rl,ybase=view,rho=rho)
        reference=memories[-1]
        pointers=saved['pointers']
        df=d['fwd']-d['fixup_extra']['fwd']-K.fwd_count(6*n)
        db=d['bwd']-d['fixup_extra']['bwd']-K.bwd_count(6*n)
        good=(df==-n and db==-24*n+2 and saved['initializer']==47*n+1
            and d['frac_mis']==0 and d['acc_ok'] and d['records_ok']
            and d['sp_kept'] and d['keys_formula'] and d['final_R_times_rho_is_one']
            and not reference.view_writes and all(reference[a]==v for a,v in pointers.items()))
        assert good,d
        mem=Guarded(rp,rl);mem.update(saved['initial'])
        mem.coins=iter(E.limbs(rho));mem.protect=True
        sb=saved['sb']
        expected_eb={a:reference[a] for a in range(sb,E.EBE,8) if a in reference}
        summary={'ok':True,'n_y':n,'fwd_delta':df,'bwd_delta':db,
                 'initializer':saved['initializer'],'static_pointer_words':len(pointers),
                 'fraction_mismatches':d['frac_mis'],'records_ok':d['records_ok'],
                 'stack_restored':d['sp_kept'],'read_only':not reference.view_writes,
                 'seconds':round(time.time()-started,3)}
        # Program-specific scratch words may differ. Compare the semantic EB
        # output and records; no instruction counts in the fixture are claims.
        return summary,dict(memory=mem,reference=reference,expected_eb=expected_eb,
            pointers=pointers,view=view,n=n,rho=rho,sb=sb,xs=xs,ys=ys,
            fwd_regs=saved['fwd_regs'],leaf_program=(I,L),route_program=(rp,rl),
            result=d,read_only_exception=ReadOnlyViewWrite)
    finally:
        E.TLC,E.popk,E.offs7,E.arun,E.CE=orig
        # Keep this Device mapping active for the returned fixture. It is a
        # local interpreter assertion; it does not create a memory-map claim.

def run_cold(E,K,field_source):
    """Execute actual tbldk recurrences and normal/mining yround branches."""
    import time
    started=time.time()
    fp,fl=E.a64.a64_load(field_source)
    old_offs,old_ce=E.offs7,E.CE
    E.CE=K
    rng=E.Rd(20261123)
    PB,SK,AB,DST=0x600000,0x500000,0x300000,0x700000
    tables=[]
    try:
        for cl in K.EIS:
            mem={};rows=3
            for i,v in enumerate(E.eiw(K.EIS[cl])):mem[PB+8*i]=v
            for a in range(SK-4096,SK,8):mem[a]=0
            steps=[rng.rr(E.p) for _ in range(6)]
            starts=[rng.rr(E.p) for _ in range(6)]
            for j in range(6):
                E.put(mem,AB+64*j,steps[j]);E.put(mem,AB+64*j+32,starts[j])
            count=E.arun(fp,fl,'tbldk',dict(x0=DST,x1=AB,x2=rows,x3=PB,sp=SK),mem)
            coeff=links=True
            for row in range(rows):
                a=DST+256*row
                links &= mem[a+120]==a and mem[a+248]==a+128
                for j in range(6):
                    c1,c2=K.eis_red(K.EIS[cl],(starts[j]+row*steps[j])%E.p)
                    coeff &= E.get(mem,a+E.KOFF[j],3)==c1%2**192 and E.get(mem,a+E.KOFF[j]+24,2)==c2
            tables.append(dict(class_name=cl,rows=rows,count=count,
                extra=count-K.tbldk_count(rows,cl),coefficients_ok=bool(coeff),links_ok=bool(links)))
        def linked_offs(v,i):
            out=old_offs(v,i);out[0]+=64;return out
        E.offs7=linked_offs
        class Costs:
            YR_BLK=K.YR_BLK+99
            def __getattr__(self,name):return getattr(K,name)
        yr,mem,*chains=E.ygen(fp,fl,rng,Costs())
        ok=(all(x['extra']==x['rows'] and x['coefficients_ok'] and x['links_ok'] for x in tables)
            and yr['centres_ok'] and yr['counts_ok'] and yr['mismatches']==0)
        assert ok,(tables,yr)
        return dict(ok=True,tables=tables,yround=yr,
            normal_y_extra=3,mining_y_extra=0,seconds=round(time.time()-started,3))
    finally:E.offs7,E.CE=old_offs,old_ce
