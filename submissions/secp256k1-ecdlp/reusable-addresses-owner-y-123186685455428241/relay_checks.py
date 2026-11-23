"""Reusable K256 relay checks. Candidate execution belongs in llm-verifier."""
import re,time
TB=1<<40;MASK=(1<<64)-1;K=256;STRIDE=12288
UP=(1<<64)-(1<<55);UPN=2932031007402;UPEND=UP+UPN*STRIDE
LOW=18013298997846016;LOWN=(1<<42)-UPN;LOWEND=LOW+LOWN*STRIDE
RCB=256*TB+65536;RCAP=277059750328320;REND=RCB+64*RCAP
TABLE=(1<<48)-1024


def block(d):return UP+d*STRIDE if d<UPN else LOW+(d-UPN)*STRIDE


def reference_geometry(source):
    lo,hi=source.index('.globl h2s0'),source.index('.macro LNONE')
    source=source[:lo]+source[lo:hi].replace('mov x18,#4096','mov x18,#8192')+source[hi:]
    lo,hi=source.index('.macro RSTUB'),source.index('RSTUB x20,x21')
    source=source[:lo]+source[lo:hi].replace('sub x9,x11,#4064','sub x9,x11,#8192\nadd x9,x9,#32')+source[hi:]
    lo,hi=source.index('.globl fin\n'),source.index('.macro H2B')
    source=source[:lo]+source[lo:hi].replace('add x8,x8,#0x2000','add x8,x8,#0x3000')+source[hi:]
    return re.sub(r'0x(?:5000|10000)\b','0x20000',source)


def run(E,source,baseline,mode='compact',export=False):
    A=E.a64;M=E.machine;started=time.time()
    if 'mov x18,#4096'in baseline[baseline.index('.globl h2s0'):baseline.index('.macro LNONE')]:baseline=reference_geometry(baseline)
    ops,labs=A.a64_load(source);bi,bl=A.a64_load(baseline)
    assert UPN+LOWN==1<<42 and UPEND==(1<<64)-8192 and LOWEND==M.LINK
    assert REND%65536==0 and REND<=LOW and LOW-REND==57344
    assert len({block(d) for d in (0,UPN-1,UPN,(1<<42)-1)})==4
    class Guard(A.S1Mem):
        def ro(self,a):
            for lo,hi in self.rg:
                if lo<=a<hi:return(a-lo)%STRIDE>=8192
            return False
    class Port(M.Port1):
        def __init__(self,mem):
            super().__init__(mem);self.msgs=[]
        def store(self,a,v):
            if a==M.MB2+40:self.send(a,M.View(self.m,M.MB2,6,0),self.msgs)
            elif any(lo<=a<hi and(a-lo)%STRIDE==0 for lo,hi in ((UP,UPEND),(LOW,LOWEND))):
                self.send(a,M.View(self.m,a,4+4*K),self.blocks)
    def init_mem(destinations,changed,cap,remaining=2):
        mem=Guard(M.dd(int),((UP,UPEND),(LOW,LOWEND)),M.AL,'vbarr')
        init=E.arun(ops,labs,'relay_source_init',{},mem)
        assert init==262
        for d in destinations:
            mem[M.LC2+8*d]=block(d)+32
            mem[(M.LC2+8*d)|M.RC2]=M.Node.R2
        mem[M.PW]=M.Node.E2
        reg={};pr,lb=(ops,labs)if changed else(bi,bl)
        E.arun(pr,lb,'h2s0',dict(x12=cap,x0=0,x8=remaining),mem,out_regs=reg)
        return mem,reg
    def rows(n,d,seed=0):return[(d|((seed+j)<<42),0x1234567890+seed+j,0xABC000+seed+j,0xDEF000+seed+j)for j in range(n)]
    def send(mem,reg,records,changed):
        pr,lb=(ops,labs)if changed else(bi,bl);port=Port(mem);out={};n=len(records)
        msg=[0,32*n,1,0]+[v for row in records for v in row]+[M.STL]*4*(128-n)
        count=E.arun(pr,lb,'h2s',reg,M.LV(mem,msg,[]),out_regs=out,dev=port.dev,
            ls64=([M.WIN,(TABLE,TABLE+512)],[]))
        state=(out['x14']<<64|out['x13'],out['x8'if changed else'x5']<<64|out['x15'],out['x28'],mem[M.DONE])
        if changed:assert out['sp']==TABLE and out['x16']==M.PW and out['x0']==M.HSP[0]
        return out,count,state,port
    def fault_case(d,hit,drop,n=128,rbase=None):
        rbase=M.Node.R2 if rbase is None else rbase
        cap=(rbase+(0 if drop else 65536))&MASK;results=[]
        for changed in (False,True):
            mem,reg=init_mem([d],changed,cap)
            # ExactRC2address96TiB: machine.RC2 is64TiB ORmask, notbase.
            rc=(M.LC2+8*d)|M.RC2
            mem[rc]=rbase
            assert mem[rc]==rbase
            start=255-hit
            for j,row in enumerate(rows(start,d,10000)):
                for w,v in enumerate(row):mem[block(d)+32+32*j+8*w]=v
            mem[M.LC2+8*d]=block(d)+32+32*start
            rec=rows(n,d)
            out,count,state,pt=send(mem,reg,rec,changed)
            assert len(mem.log)==1
            el,esr,far,spsr=mem.log[0]
            assert far==block(d)+8192 and esr==0x9600004F
            assert len(pt.blocks)==int(not drop),('send/drop not exercised',d,hit,drop,len(pt.blocks))
            if not drop:
                expected=rows(start,d,10000)+rec[:hit+1]
                assert pt.blocks[0][1:4]==[rbase,8192,M.Node.E2]
                assert pt.blocks[0][4:]==[v for r in expected for v in r]
            assert mem[rc]==cap if drop else mem[rc]==(rbase+8192)&MASK
            assert mem[M.LC2+8*d]==block(d)+32+32*(n-hit-1)
            results.append((count,state,dict(mem.m),pt.blocks,pt.msgs))
        a,b=results;assert a[1:]==b[1:]
        assert b[0]-a[0]==(-43 if n==128 else 0)
        return dict(d=d,hit=hit,drop=drop,n=n,rbase=rbase,delta=b[0]-a[0],sent=int(not drop))
    plain=[]
    ns=range(129)if mode=='plain'else (0,1,2,63,64,127,128)
    for n in ns:
        for remaining in (1,2):
            z=[]
            for changed in (False,True):
                mem,reg=init_mem([0],changed,M.Node.R2+65536,remaining)
                out,count,state,pt=send(mem,reg,rows(n,0),changed)
                assert not mem.log and not pt.blocks
                z.append((count,state,dict(mem.m),pt.msgs))
            a,b=z;assert a[1:]==b[1:]
            expect=-43 if n==128 else int(remaining==1)
            assert b[0]-a[0]==expect
            plain.append(dict(n=n,remaining=remaining,delta=expect))
    faults=[]
    ds=(UPN-1,)if mode=='upper'else(UPN,)if mode=='lower'else(0,UPN)
    hs=range(128)if mode in('upper','lower')else(0,1,63,64,126,127)
    for d in ds:
        for hit in hs:
            for drop in (False,True):
                for n in ((128,127)if hit<127 else(128,)):
                    faults.append(fault_case(d,hit,drop,n))
    for d in (0,UPN):
        for drop in (False,True):faults.append(fault_case(d,127,drop,rbase=(1<<64)-8192))
    # fin must preserveRC2'scompletedblockcursor andLC2'spartialcursor.
    finals=[]
    for d,rbase in ((UPN-1,M.Node.R2),(UPN,M.Node.R2),(UPN-1,(1<<64)-8192),(UPN,(1<<64)-8192)):
        for changed in (False,True):
            cap=(rbase+65536)&MASK;mem,reg=init_mem([d],changed,cap,100)
            mem[(M.LC2+8*d)|M.RC2]=rbase
            allrows=[];sent=[];watermarks=[];counts=[]
            for group in range(2):
                new=rows(300,d,group*300);allrows+=new
                for j in range(0,300,128):
                    reg,c,_,pt=send(mem,reg,new[j:j+128],changed);sent+=pt.blocks
                before=(mem[M.LC2+8*d],mem[(M.LC2+8*d)|M.RC2])
                pt=Port(mem);out={};pr,lb=(ops,labs)if changed else(bi,bl)
                c=E.arun(pr,lb,'fin',dict(x0=M.MB2,x7=M.LC2+8*d,x3=(M.LC2+8*d)|M.RC2,
                     x8=block(d)+32,x1=1,x2=M.Node.E2,x6=cap,x11=M.Node.E2),mem,dev=pt.dev,out_regs=out)
                assert before==(mem[M.LC2+8*d],mem[(M.LC2+8*d)|M.RC2])
                actual=(rbase+len(allrows)*32)&MASK
                assert pt.msgs==[[M.LC2+8*d,M.Node.E2,actual,0,0,0]]
                assert len(pt.blocks)==1
                q=len(allrows)//K*K;tail=allrows[q:]
                assert pt.blocks[0][1:4]==[(rbase+32*q)&MASK,32*len(tail),M.Node.E2],pt.blocks[0][:4]
                assert pt.blocks[0][4:4+4*len(tail)]==[v for r in tail for v in r]
                watermarks.append(actual);counts.append(c)
            assert len(sent)==2
            assert sent[1][4:]==[v for r in allrows[256:512]for v in r]
            assert mem[(M.LC2+8*d)|M.RC2]==(rbase+512*32)&MASK
            assert mem[M.LC2+8*d]==block(d)+32+88*32
            finals.append(dict(d=d,changed=changed,rbase=rbase,watermarks=watermarks,fin_counts=counts,old_tail_preserved=True))
    # Explicittwo-runinit3+finatthesplit, withoneupperandtwolowerdestinations.
    mem=Guard(M.dd(int),((UP,UPEND),(LOW,LOWEND)),M.AL,'vbarr');cur=M.LC2+8*(UPN-1);init_counts=[]
    for count,base in ((1,block(UPN-1)),(2,LOW)):
        out={};z=E.arun(ops,labs,'init3',dict(x0=count,x1=cur,x2=base+32,x3=STRIDE),mem,out_regs=out)
        init_counts.append(z);cur=out['x1']
    assert [mem[M.LC2+8*d]for d in (UPN-1,UPN,UPN+1)]==[block(d)+32 for d in (UPN-1,UPN,UPN+1)]
    commoncap=M.Node.R2+65536
    splitdest=(UPN-1,UPN,UPN+1);splitn=(0,1,255)
    for j,(d,n)in enumerate(zip(splitdest,splitn)):
        mem[M.LC2+8*d]=block(d)+32+32*n
        mem[(M.LC2+8*d)|M.RC2]=commoncap if j==2 else M.Node.R2+j*8192
        for t,row in enumerate(rows(n,d)):
            for w,v in enumerate(row):mem[block(d)+32+32*t+8*w]=v
    oldcursors=[(mem[M.LC2+8*d],mem[(M.LC2+8*d)|M.RC2])for d in splitdest]
    fcur=M.LC2+8*(UPN-1);rcur=fcur|M.RC2;split_fin_counts=[];split_msgs=[];split_blocks=[]
    for n,base in ((1,block(UPN-1)),(2,LOW)):
        pt=Port(mem);out={}
        z=E.arun(ops,labs,'fin',dict(x0=M.MB2,x7=fcur,x3=rcur,x8=base+32,x1=n,
            x2=M.Node.E2,x6=commoncap,x11=M.Node.E2),mem,dev=pt.dev,out_regs=out)
        split_fin_counts.append(z);split_msgs+=pt.msgs;split_blocks+=pt.blocks
        fcur,rcur=out['x7'],out['x3']
    assert(fcur,rcur)==(M.LC2+8*(UPN+2),(M.LC2+8*(UPN+2))|M.RC2)
    assert oldcursors==[(mem[M.LC2+8*d],mem[(M.LC2+8*d)|M.RC2])for d in splitdest]
    assert [m[0]for m in split_msgs]==[M.LC2+8*d for d in splitdest]
    assert [m[2]for m in split_msgs]==[M.Node.R2,M.Node.R2+8192+32,commoncap]
    assert len(split_blocks)==1 and split_blocks[0][2]==32
    # Execute the actual constant-loading split wrappers at reduced populations.
    # Keep threeinstructions percount so fixed setup charges stay exact.
    small=source
    for reg in ('x0','x1'):
        upper=f'mov {reg},#43690\nmovk {reg},#43690,lsl #16\nmovk {reg},#682,lsl #32'
        lower=f'mov {reg},#21846\nmovk {reg},#21845,lsl #16\nmovk {reg},#341,lsl #32'
        assert source.count(upper)==source.count(lower)==1
        small=small.replace(upper,f'mov {reg},#1\nmov {reg},{reg}\nmov {reg},{reg}')
        small=small.replace(lower,f'mov {reg},#2\nmov {reg},{reg}\nmov {reg},{reg}')
    si,sl=A.a64_load(small);sm=Guard(M.dd(int),((UP,UPEND),(LOW,LOWEND)),M.AL,'vbarr')
    split_init_wrapper=E.arun(si,sl,'relay_blocks_init_split',{},sm)
    assert split_init_wrapper==33
    assert [sm[M.LC2+8*j]for j in range(3)]==[UP+32,LOW+32,LOW+STRIDE+32]
    for j in range(3):sm[(M.LC2+8*j)|M.RC2]=M.Node.R2
    sp=Port(sm);sr={}
    split_fin_wrapper=E.arun(si,sl,'relay_fin_split',dict(x0=M.MB2,x2=M.Node.E2,x6=commoncap,x11=M.Node.E2),sm,dev=sp.dev,out_regs=sr)
    assert split_fin_wrapper==69 and len(sp.msgs)==3 and not sp.blocks
    assert(sr['x7'],sr['x3'])==(M.LC2+24,(M.LC2+24)|M.RC2)
    # Guardfootprintandnegativecontrol: page1iswritable,page2RO; oldbit12
    # rulewouldfaultonrecord127ratherthan255, which thischeckrejects.
    assert all(not mem.ro(block(d)+32+32*j)for d in(0,UPN)for j in range(255))
    assert all(mem.ro(block(d)+8192)for d in(0,UPN))
    old_guard_wrong=any(bool((block(d)+32+32*j)&4096)!=mem.ro(block(d)+32+32*j)for d in(0,UPN)for j in range(256))
    assert old_guard_wrong
    # Realrequestcapboundary: one accepted thenone rejected request.
    mem2=M.dd(int);curaddr=256*TB+0x400;mem2[curaddr]=REND-64
    hdr=M.LINK+0x9000
    for j in range(6):mem2[hdr+8*j]=0x7770+j
    count=E.arun(ops,labs,'dlvq',dict(x0=hdr,x10=curaddr,x11=REND,x9=64,x29=M.CNT,x28=1),mem2)
    assert mem2[curaddr]==REND and[mem2[REND-64+8*j]for j in range(5)]==[0x7771,0x7772,0x7773,0x7774,0x7775]
    prior=dict(mem2)
    count2=E.arun(ops,labs,'dlvq',dict(x0=hdr,x10=curaddr,x11=REND,x9=64,x29=M.CNT,x28=1),mem2)
    assert mem2[curaddr]==REND and all(a<REND or a>=LOW for a in mem2 if a not in prior)
    assert REND not in mem2
    saved=ops,labs
    bad=source.replace('sub x9,x11,#8192\nadd x9,x9,#32','sub x9,x11,#4096\nadd x9,x9,#32')
    assert bad!=source
    ops,labs=A.a64_load(bad)
    try:fault_case(0,0,False)
    except (AssertionError,KeyError)as exc:stub_negative=str(exc)[:180]
    else:raise AssertionError('wrong guard displacement not detected')
    ops,labs=saved
    result=dict(ok=True,plain_cases=len(plain),fault_cases=len(faults),sent_fault_cases=sum(f['sent']for f in faults),
        dropped_fault_cases=sum(not f['sent']for f in faults),full_source_delta=-43,partial_last_final_extra=1,
        finals=finals,split_init_counts=init_counts,split_fin_counts=split_fin_counts,
        split_wrapper_counts=[split_init_wrapper,split_fin_wrapper],
        split_fin_cursor_preserved=True,request_cap=RCAP,request_cap_end=REND,
        request_counts=[count,count2],guard_negative_detected=old_guard_wrong,route_slots=len(ops),
        stub_negative_detected=stub_negative,route_bytes=4*len(ops),seconds=round(time.time()-started,3))
    if export:
        import base64,lzma
        targets=set();code=[];branches=('b','bl','cbz','cbnz','tbz','tbnz','cbhs','cbhi','cbeq','cbne','cblt','cbge')
        for p,(op,args0)in enumerate(ops):
            args=list(args0)
            if op in branches or op.startswith('b.'):
                t=args[-1];num=t[:-1].isdigit()and t[-1:]in('f','b');where=labs[t[:-1]if num else t]
                dest=next(q for q in where if q>p)if num and t[-1]=='f'else next(q for q in reversed(where)if q<=p)if num else where[0]
                targets.add(dest);args[-1]='label_'+str(dest)
            code.append(op+' '+','.join(args))
        flat=[]
        for p,line in enumerate(code):
            if p in targets:flat.append('label_'+str(p)+':')
            flat.append(line)
        if len(code)in targets:flat.append('label_'+str(len(code))+':')
        result['flat_lzma_base64']=base64.b64encode(lzma.compress(('\n'.join(flat)+'\n').encode())).decode()
    return result
