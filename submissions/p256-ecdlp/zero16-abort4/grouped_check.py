import random,a64
from pathlib import Path
from routing_check import Node,run_k,range_scan_model,LinkView,LINK,SENTINEL
ROOT=Path(__file__).resolve().parent
C4C=0x5200000
def run(last):
 prog,lab=a64.a64_load((ROOT/'route.s').read_text());rg=random.Random(97+last)
 side,n3,n4,cap3,ncl=2,4,8,20,64;st=32*cap3;H=[(a,j) for a in range(side) for j in range(side)]
 V={h:Node() for h in H};R={h:dict(x2=V[h].CT,x4=V[h].C4) for h in H};calls=[];bk={h:[[] for _ in range(n3)] for h in H};ok=True;found=-1;ran=0
 for h,nd in V.items():
  run_k(prog,lab,'init3',dict(x0=n3,x1=nd.CT,x2=nd.BK,x3=st),nd,calls,5*n3+1)
  run_k(prog,lab,'zero4',dict(x0=n4,x1=nd.C4),nd,calls,n4//4*4+1)
 for g in range(4):
  ran+=1
  for h,nd in V.items():
   nd.mem[Node.DC]=0;dc=0;pk=[]
   for s in range(side):
    for t in range(6):
     k2=rg.choice([5<<26,9<<26|2,rg.getrandbits(38)<<26|2]);k3=rg.getrandbits(16)<<48|rg.randrange(n4);w=0
     if h==H[0] and (g,s,t) in ((0,0,0),(last,1,0)):k2,k3,w=77<<26|1,44<<48|3,1
     pk.append([0,0,k2,k3,w<<63|g<<40|s<<8|t,s])
    pk.append([s,1,0,0,0,0])
   for m in pk:
    o={};done=m[1]==1;d3=m[2]&(n3-1);acc=not done and len(bk[h][d3])<cap3
    if acc:bk[h][d3].append((nd.BK+d3*st+32*len(bk[h][d3]),tuple(m[2:])))
    c=a64.a64_run(prog,lab,'receive',dict(x0=LINK,x2=R[h]['x2'],x4=R[h]['x4'],x13=nd.DC,x15=dc,x16=side,x17=1,x18=32,x19=(1<<64)-32,x30=SENTINEL,sp=nd.SP),LinkView(nd.mem,m),out_regs=o)
    R[h]=dict(x2=o['x2'],x4=o['x4']);fin=done and dc+1==side;calls.append(('receive',c,(10 if fin else 8) if done else 14 if acc else 11));dc=o['x15'] if done else dc
   ok&=nd.mem[Node.DC]==1 and dc==0
  for h,nd in V.items():
   seq=[e for t in range(n3) for e in bk[h][t]];cnt=[sum(1 for e in seq if e[1][1]&(n4-1)==d) for d in range(n4)]
   kp={a:nd.mem.get(a) for a in [nd.CT+8*i for i in range(2*n3)]+[nd.C4+8*i for i in range(n4)]+[e[0]+8*k for e in seq for k in range(4)]}
   ok&=[nd.mem[nd.C4+8*d] for d in range(n4)]==cnt
   run_k(prog,lab,'prefix4c',dict(x0=n4,x1=nd.C4,x2=nd.PA,x3=C4C),nd,calls,n4//4*10+1)
   run_k(prog,lab,'pass4',dict(x0=n3,x1=nd.CT,x2=nd.BK,x3=st,x4=C4C,x5=n4-1),nd,calls,11*n3+6*len(seq)+2*sum(len(b)//8 for b in bk[h])+1)
   order=sorted(seq,key=lambda e:e[1][1]&(n4-1));ok&=[nd.mem[nd.PA+8*k] for k in range(len(seq))]==[e[0] for e in order]
   mc,cands=range_scan_model([[e[1][0] for e in order if e[1][1]&(n4-1)==d] for d in range(n4)],ncl)
   o=run_k(prog,lab,'scan_ranges',dict(x0=n4,x1=nd.PA,x2=C4C,x3=nd.CL,x4=nd.CL+16*ncl),nd,calls,mc)
   got=[(nd.mem[c],nd.mem[c+8]) for c in range(nd.CL,o['x3'],16)];ok&=got==[(order[a][0],order[b][0]) for a,b in cands]
   ok&=all(nd.mem.get(a)==v for a,v in kp.items())
   tag={e[0]:e[1][2] for e in seq}
   if any(tag[x]>>63 and tag[y]>>63 and tag[x]>>40&3!=tag[y]>>40&3 for x,y in got):found=g
  if found>=0:break
 return dict(groups_run=ran,stop_group=found,records=sum(len(b) for h in H for b in bk[h]),ok=ok,formulas_ok=all(c==f for _,c,f in calls))
def check():
 s,f=run(1),run(3)
 assert s['ok'] and f['ok'] and s['formulas_ok'] and f['formulas_ok']
 assert (s['groups_run'],s['stop_group'],f['groups_run'],f['stop_group'])==(2,1,4,3)
 return dict(stop_case=s,full_case=f,scope='Synchronous finite grouped replay at toy scope.')
