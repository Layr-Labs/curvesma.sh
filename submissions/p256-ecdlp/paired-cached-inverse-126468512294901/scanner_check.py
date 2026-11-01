import json,random
import a64

def check(root):
    pr,labels=a64.a64_load((root/'scan_ranges.s').read_text())
    cfg={'seed':101};rng=random.Random(cfg['seed']);cases=0;worst=0

    def check_case(groups,cap):
     nonlocal cases,worst
     pa,ends,out,records=0x10000,0x20000,0x30000,0x40000;m={};n=0;expected=[]
     for bucket,keys in enumerate(groups):
      aa=[]
      for key in keys:
       a=records+32*n;m[pa+8*n]=a;m[a]=key;aa.append(a);n+=1
      expected.extend((aa[i-1],aa[i]) for i in range(1,len(aa)) if keys[i-1]==keys[i])
      m[ends+8*bucket]=pa+8*n
     regs={};ct=a64.a64_run(pr,labels,'scan_ranges',dict(x0=len(groups),x1=pa,x2=ends,x3=out,x4=out+16*cap,x30=0x12345),m,limit=100000,out_regs=regs)
     got=[(m[a],m[a+8]) for a in range(out,regs['x3'],16)]
     assert got==expected[:cap],(groups,cap,got,expected[:cap])
     bound=(51*n+15)//16+24*len(groups)+8*(cap+1)+2
     assert 16*ct<=51*n+16*(24*len(groups)+8*(cap+1)+2),(ct,bound,n,len(groups),cap)
     if not expected and len(groups)==1:worst=max(worst,ct-2-51*n/16)
     cases+=1
    for n in range(98):
     for cap in [0,1,2,1000]:check_case([[rng.randrange(3) for _ in range(n)]],cap)
     check_case([list(range(n))],1000)
    for b in [0,1,2,4,8,16]:
     for _ in range(10):
      g=[[rng.randrange(4) for _ in range(rng.randrange(130))] for _ in range(b)]
      for cap in [0,1,2,8,10000]:check_case(g,cap)
    return {'status':'passed','cases':cases,'exact_adjacencies_and_caps':True,'instruction_bound':'51*n/16+24*b+8*(cap+1)+2','max_single_nonmatch_overhead':worst,'scope':'Finite scanner checks only; no complete search score or success claim.'}
