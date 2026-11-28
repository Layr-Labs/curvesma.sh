"""Finite Card66 archive/record boundary checks; no production-scale run claim.
Real emitted ystage routines are executed on sparse full-width memory. The sweep
probe substitutes only its initial owner counter with the last valid owner; it
executes every active bank of that owner, including begin/send/ACK/tail control.
Network ACKs are deterministic test fixtures, not a fabric-correctness proof.
"""
import json, pathlib, sys
from collections import defaultdict


def check(root=None):
    root = pathlib.Path(root) if root is not None else next(p for p in (pathlib.Path(__file__).resolve().parent, pathlib.Path('/snapshot/files'), pathlib.Path.cwd()) if (p/'kernels/ystage.s').is_file())
    sys.path.insert(0, str(root/'kernels'))
    import a64
    source=(root/'kernels/ystage.s').read_text()+'\n'+(root/'kernels/ytransport.s').read_text()
    code, labels=a64.a64_load(source)
    U=1<<62; META=U; INBOX=U+(1<<32); CTRL=U+(1<<33)
    CMDIN=INBOX+65536; CMDCTRL=CTRL+512
    G=34<<40; SRC=G+0x7200000; CMDSRC=G+0x7600000
    SP=G+0x7000000; CB=G+0x2000000; FRAME=40<<40
    SENT=0x51515151; CAP=1<<27; QMAX=67334001488
    TAG=1<<63; result={}; controls={}
    class Mem(defaultdict):
        def __init__(self):super().__init__(int)
    def run(entry,mem,regs=None,dev=None,stop=None,alternate=None,limit=2000000):
        rr={'sp':SP,'x30':SENT};rr.update(regs or {})
        co,la=alternate or (code,labels)
        n=a64.a64_run(co,la,entry,rr,mem,dev=dev,ls64=([(0,1<<64)],[(0,1<<64)]),out_regs=rr,stop_pc=stop,limit=limit)
        return n,rr
    def ack_device(mem,ctrl,save=None):
        def wr(off,value):
            if off==24 and value:
                if save is not None:save.append((mem[ctrl],mem[ctrl+8],[mem[ctrl+4096+8*i] for i in range(16)]))
                mem[ctrl+16]=1
        return (ctrl,8192,wr,lambda:None,lambda:True)
    def put_payload(mem,base,seed):
        for i in range(56):mem[base+64+8*i]=seed+i
    # Every operation executes with its production constants; sparse backing
    # allows high addresses without allocating an exponentially large memory.
    for n in (0,1,3,17):
        mm=Mem();steps,_=run('ybank_begin',mm,{'x0':(1<<44)+19,'x1':n})
        assert all(mm[META+16*i]==0 and mm[META+16*i+8]==TAG+(1<<44)+19+i for i in range(n))
        assert mm[CTRL+104]==0 and steps==4*n+9,(n,steps)
    result['bank_begin_per_slot']=4;result['bank_begin_fixed']=9
    mm=Mem(); steps,_=run('yboot_bank',mm,{'x0':(1<<58)-1,'x1':0})
    assert mm[CTRL+224]==(1<<58)-1 and mm[CTRL+232]==0 and mm[CTRL+104]==2
    result['bank_boot']=steps
    # The upper slot and former alias boundary both traverse real receiver code.
    receiver_steps=[]
    slots=(0,(1<<19)-1,1<<19,(1<<27)-1,QMAX%CAP)
    for slot in slots:
        mm=Mem(); tag=TAG+(1<<44)+slot
        mm[CTRL+24]=1;mm[CTRL+96]=0;mm[META+16*slot+8]=tag
        mm[INBOX]=slot;mm[INBOX+8]=tag;mm[INBOX+16]=7;mm[INBOX+24]=0
        put_payload(mm,INBOX,0xA100+slot)
        steps,_=run('yreceive',mm);receiver_steps.append(steps)
        assert mm[CTRL+48]==1 and mm[META+16*slot]==1
        target=slot<<35
        assert target+512<=1<<62
        assert all(mm[target+64+8*i]==0xA100+slot+i for i in range(56))
    for count in ((1<<26)-2,(1<<26)-1):
        mm=Mem();slot=CAP-1;tag=TAG+91
        mm[CTRL+24]=1;mm[META+16*slot]=count;mm[META+16*slot+8]=tag
        mm[INBOX]=slot;mm[INBOX+8]=tag;put_payload(mm,INBOX,0xE500)
        steps,_=run('yreceive',mm);receiver_steps.append(steps)
        target=(slot<<35)+(count<<9)
        assert target+512<=1<<62
        assert mm[CTRL+48]==1 and mm[META+16*slot]==min(count+1,(1<<26)-1)
        assert all(mm[target+64+8*i]==0xE500+i for i in range(56))
    result['receiver_cases']=len(receiver_steps);result['receiver_insns']=receiver_steps
    for slot,wrong,sealed in [(1<<27,False,False),(1<<26,True,False),(1<<26,False,True)]:
        mm=Mem();tag=TAG+42;mm[CTRL+24]=1;mm[CTRL+104]=int(sealed)
        mm[INBOX]=slot;mm[INBOX+8]=tag;mm[META+slot*16+8]=tag+int(wrong)
        run('yreceive',mm);assert mm[CTRL+48]==2
    controls['receiver_rejects_oob_stale_sealed']=True
    old=a64.a64_load(source.replace('lsr x8,x19,#27','lsr x8,x19,#19'))
    mm=Mem();slot=1<<19;mm[CTRL+24]=1;mm[INBOX]=slot;mm[INBOX+8]=TAG+7;mm[META+16*slot+8]=TAG+7
    run('yreceive',mm,alternate=old);assert mm[CTRL+48]==2
    controls['old_receiver_rejects_valid_new_slot']=True
    # Real yproduce builds the packet and sends it, including the high-k w3.
    producer_steps=[]
    for rel in (0,(1<<19),1<<29,QMAX-1,(1<<37)-1):
        mm=Mem();first=123456;k=first+rel;c=(1<<48)-1;u=(1<<21)-1
        mm[CB+2240]=FRAME+512;mm[CB+2248]=k;mm[CB+2256]=c
        mm[CB+2264]=u;mm[CB+2272]=first;mm[CB+2280]=(1<<58)-1
        mm[FRAME+56]=c>>15;put_payload(mm,FRAME,0xB200)
        packets=[];steps,_=run('yproduce',mm,{'x29':CB},dev=ack_device(mm,SRC,packets));producer_steps.append(steps)
        expected=(rel>>27,rel&(CAP-1))
        assert packets[0][:2]==(c,expected[0])
        assert packets[0][2][0:2]==[expected[1],TAG|k]
        expectw3=((k-1)>>36)|((c>>15)<<9)|(u<<42)
        assert mm[FRAME+360]==expectw3 and mm[SRC+4096+360]==expectw3
        assert mm[CB+2264]==u+1
    result['producer_cases']=5;result['producer_insns']=producer_steps
    old=a64.a64_load(source.replace('ubfx x9,x8,#27,#10','ubfx x9,x8,#19,#10').replace('and x8,x8,#0x7ffffff','and x8,x8,#0x7ffff'))
    mm=Mem();mm[CB+2240]=FRAME+512;mm[CB+2248]=(1<<29)+1;mm[CB+2272]=1
    packets=[];run('yproduce',mm,{'x29':CB},dev=ack_device(mm,SRC,packets),alternate=old)
    assert packets[0][1]==0 and packets[0][2][0]==0
    controls['old_producer_aliases_2p29_with_zero']=True
    # Replay at the highest valid slot with real sender and one full payload.
    mm=Mem();slot=CAP-1;tag=TAG+(1<<44)+17;owner=(1<<48)-1;target=slot<<35
    mm[CTRL+104]=1;mm[CTRL+224]=(owner<<10)|1023;mm[CTRL+232]=0
    mm[META+16*slot]=1;mm[META+16*slot+8]=tag;put_payload(mm,target,0xC300)
    packets=[];steps,_=run('yreplay_request',mm,{'x0':slot,'x1':tag,'x2':owner},dev=ack_device(mm,CTRL,packets))
    assert packets[0][:2]==(owner,0) and packets[0][2][:2]==[0,tag]
    assert all(mm[CTRL+4096+64+8*i]==0xC300+i for i in range(56))
    assert mm[CTRL+128]==target+512 and mm[CTRL+56]==1 and mm[CTRL+192]==0
    result['replay_one_payload_insns']=steps
    # Actual command receive->begin->ACK path; a short population checks the
    # body while the sweep test below checks production-size command arguments.
    mm=Mem();mm[CMDCTRL+24]=1;mm[CTRL+104]=2
    for off,v in [(0,1),(8,1<<44),(16,3),(24,99),(32,1<<58),(40,0)]:mm[CMDIN+off]=v
    steps,_=run('ycmd_bank',mm)
    assert mm[CTRL+112]==99 and mm[CMDCTRL+48]==1
    assert [mm[META+16*i+8] for i in range(3)]==[TAG+(1<<44)+i for i in range(3)]
    result['command_begin3_insns']=steps
    # Owner's request construction uses same quotient/remainder as producer.
    a,b=labels['yowner_round'][0],labels['yowner_group'][0]
    stop=next(i for i in range(a,b) if code[i]==('ldr',['x3','[sp,#16]']))
    for rel in (1<<29,QMAX-1):
        mm=Mem();packets=[];steps,_=run('yowner_round',mm,{'x0':123+rel,'x1':123,'x2':(1<<48)-1},dev=ack_device(mm,CMDSRC,packets),stop=stop)
        expected_bank,expected_slot=divmod(rel,CAP)
        assert packets[0][:2]==((((1<<48)-1)<<10)|expected_bank,0xc0000)
        assert packets[0][2][0:4]==[4,expected_slot,TAG|(123+rel),(1<<48)-1]
    result['owner_request_prefix_insns']=steps
    # Full bank-count sweep for one boundary owner. Only initial owner index is
    # replaced; all998/1024 bank route/slot-count computations execute unedited.
    class SweepMem(Mem):
        def instruction(self,pc,op,args,regs,event):
            if pc==labels['ybs_owner'][0]:self[regs['sp']+48]=(1<<48)-1
    for q in (QMAX,1<<37):
        mm=SweepMem();packets=[];first=(1<<44)+1
        steps,_=run('ybank_sweep',mm,{'x0':1,'x1':first,'x2':first+q-1},dev=ack_device(mm,CMDSRC,packets))
        nb=(q+CAP-1)//CAP
        assert len(packets)==nb
        for i,(lo,hi,packet) in enumerate(packets):
            assert (lo,hi)==((((1<<48)-1)<<10)|i,0xc0000)
            assert packet[:4]==[1,first+i*CAP,min(CAP,q-i*CAP),first]
            assert packet[4:6]==[1<<58,0]
        result['sweep_%s_banks_insns'%nb]=steps
    mm=Mem();trapped=False
    try:run('ybank_sweep',mm,{'x0':1,'x1':1,'x2':(1<<37)+1})
    except AssertionError as e:trapped='UDF' in str(e)
    assert trapped;controls['sweep_rejects_over_2p37']=True
    # Both forwarding entrypoints are extracted unchanged with their exact
    # GENERIC_COMMIT/LINE64END macros. Only generic_send is a packet-capture stub.
    route=(root/'kernels/route.s').read_text()
    macro=route[:route.index('.endm')+len('.endm')]
    start=route.index('.macro LINE64END');end=route.index('.endm',start)+len('.endm')
    body=route[route.index('.globl fwdq\n'):route.index('.globl sendr\n')]
    stub="""\n.text\ngeneric_send:\nmov x10,#0x10000\nldp x11,x12,[x0]\nstp x11,x12,[x10]\nldp x11,x12,[x0,#48]\nstp x11,x12,[x10,#16]\nret\n"""
    forwarding=macro+'\n'+route[start:end]+'\n'+body+stub
    fq=a64.a64_load(forwarding)
    def forward(entry,w2,w3,alt=fq):
        mm=Mem();mm[0x20008]=w2;mm[0x20010]=w3
        run(entry,mm,{'x0':0x30000,'x1':0x20000,'x6':0x20040,'x7':42},alternate=alt)
        assert mm[0x10010]==w3
        return mm[0x10000],mm[0x10008]
    count=0
    for c in (0,32767,32768,(1<<48)-1):
        for member in (0,1023):
            for k in (1,1<<36,(1<<36)+1,27904000000000,1<<45):
                for u in (0,(1<<21)-1):
                    w2=(member<<3)|((c&32767)<<13)|(((k-1)&((1<<36)-1))<<28)
                    w3=((k-1)>>36)|((c>>15)<<9)|(u<<42)
                    assert forward('fwdq1',w2,w3)==(1024*c+member,0)
                    assert ((w2>>28)|((w3&511)<<36))+1==k
                    count+=1
    oldfq=a64.a64_load(forwarding.replace('ubfx x5,x3,#9,#33','ubfx x5,x3,#1,#41'))
    assert forward('fwdq1',0,512,oldfq)==(1<<33,0)
    controls['old_ordinary_decoder_wrong_value']=True
    mining_max=256*23750000*256*32736-1
    for g in (0,(1<<58)-1):
        for idx in (0,1<<54,mining_max):
            assert forward('fwdq',g,idx<<2)==(g,0)
    result['ordinary_forwarding_cases']=count
    result['mining_forwarding_cases']=6
    result['maximum_mining_index']=mining_max
    result['controls']=controls;result['all_checks_pass']=True
    return result

if __name__=='__main__':print(json.dumps(check(),sort_keys=True))
