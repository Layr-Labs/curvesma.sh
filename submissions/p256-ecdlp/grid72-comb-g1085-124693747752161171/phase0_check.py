"""Phase 0 (proof-appendices.md PHASE0) on real P-256 points: mine() is the home dictionary (heads on x bits, links, visit
and list caps), psel() the collector (pair cap, merge sort, abort = None). Finite evidence only."""
import random
import scalar_check as sc
P,ORD=sc.P,sc.R
BC=(sc.G[1]**2-sc.G[0]**3+3*sc.G[0])%P
PC,PD,PH,HCAP,LCAP=102,64,150,48,64# cell x>>102, close 0<x'-x<2^64, head bits from 150, caps
def real_x(x):
 if not 0<=x<P:return None
 w=(x*x*x-3*x+BC)%P;y=pow(w,(P+1)//4,P)
 return (x,y)if y*y%P==w else None
def near_x(x0,lo,hi,rg):
 for _ in range(4000):
  Z=real_x(x0+lo+rg.randrange(hi-lo))
  if Z:return Z
def opairs(cands,cw=PC,dw=PD):
 # reference, by sorting
 l=sorted(cands)
 return[(l[i],l[k])for i in range(len(l))for k in range(i+1,len(l))if l[k][0]-l[i][0]<1<<dw and l[k][0]>>cw==l[i][0]>>cw]
def mine(cands,cw=PC,dw=PD,hs=PH,hb=6,hc=HCAP,lc=LCAP):
 homes,flag,most={},0,0
 for x,t in cands:# arrival order
  hd,nx,rec,ls=homes.setdefault(x>>hs+hb,([-1]*(1<<hb),[],[],[]))
  h=x>>hs&(1<<hb)-1;j=hd[h];n=0
  while j>=0:
   if n==hc:flag=1;break# an earlier record past the visit cap
   n+=1;z=rec[j][0]
   if z==x or z>>cw==x>>cw and abs(x-z)<1<<dw and len(ls)==lc:flag=1# equal x, or a pair past the list cap
   elif z>>cw==x>>cw and abs(x-z)<1<<dw:ls.append(tuple(sorted((rec[j],(x,t)))))
   j=nx[j]
  most=max(most,n);nx.append(hd[h]);hd[h]=len(rec);rec.append((x,t))
 return homes,flag,most
CMPS=[0,0]
def msort(v):
 # merge sort by lower x; comparisons against n ceil(log2 n)
 w,n=1,len(v);CMPS[1]+=n*max(n-1,0).bit_length()
 while w<n:
  o=[]
  for i in range(0,n,2*w):
   a,b=v[i:i+w],v[i+w:i+2*w]
   while a and b:
    CMPS[0]+=1;o.append((b if b[0][0][0]<a[0][0][0]else a).pop(0))
   o+=a+b
  v,w=o,2*w
 return v
def psel(cands,need,pc=1<<20,**k):
 homes,flag,_=mine(cands,**k)
 pr=[x for v in homes.values()for x in v[3]]
 used=[t for a,b in pr for t in(a[1],b[1])]
 if flag or len(pr)>pc or len(set(used))<len(used)or len(pr)<need:return None
 return msort(pr)[:need]
def generator(rg,n=80,cw=250,dw=246):
 # real [s]P, toy widths: 4 homes of 4 heads of 4 cells
 c=[(sc.mul(rg.randrange(1,ORD))[0],i)for i in range(n)]
 k=dict(cw=cw,dw=dw,hs=252,hb=2,hc=n,lc=n)
 pr=opairs(c,cw,dw);f=len(pr);o=c[:];rg.shuffle(o)
 homes,flag,most=mine(c,**k);s1=psel(c,f,**k)
 ts=[t for a,b in pr for t in(a[1],b[1])];sh=len(set(ts))<len(ts)
 hh=[x>>252 for x,_ in c];full=max(map(hh.count,hh));lm=max(len(v[3])for v in homes.values())
 res=dict(candidates=n,found=f,shared_candidate=sh,lists_equal_reference=not flag and sorted(x for v in homes.values()for x in v[3])==pr,
  heads_with_several_cells=sum(len({x>>cw for x,_ in c if x>>252==h})>1 for h in set(hh)),fullest_head=full,
  visits_at_most_fullest_less_one=most==full-1,head_cap_aborts=mine(c,**dict(k,hc=full-2))[1]==1 and mine(c,**dict(k,hc=full-1))[1]==0,
  list_cap_aborts=lm==0 or mine(c,**dict(k,lc=lm-1))[1]==1,collector_cap_aborts=sh or f==0 or psel(c,f,pc=f-1,**k)is None,
  sort_is_reference_order=sh or s1==pr,sort_comparisons_within_n_log_n=CMPS[0]<=CMPS[1],
  arrival_order_invariant=s1==psel(o,f,**k),abort_if_one_more_needed=psel(c,f+1,**k)is None,selected_iff_no_shared=(s1 is None)==sh)
 res['ok']=all(v for q,v in res.items()if q not in('found','shared_candidate'))
 return res
def constructed(rg):
 # production widths, 64 heads a home
 cd=[]
 for i in range(6):
  U=near_x(rg.randrange(1,P>>PC)<<PC|rg.randrange(1<<PC-1),0,1<<40,rg)
  cd+=[(U[0],('R',2*i)),(near_x(U[0],1,1<<PD,rg)[0],('R',2*i+1))]
 n=len(opairs(cd));sel=psel(cd,n);sh=cd[:];rg.shuffle(sh)
 lo=cd[0][0];W=near_x(cd[1][0],1,1<<40,rg)
 H0=rg.randrange(1,P>>PH)<<PH
 H=[(near_x(H0+(i<<120),0,1<<40,rg)[0],('H',i))for i in range(HCAP+2)]
 Hc=(near_x(H[0][0],1,1<<62,rg)[0],('H',99))
 res=dict(pairs=n,all_pairs_selected=sel is not None and len(sel)==6,arrival_order_invariant=sel==psel(sh,n),abort_short=psel(cd,n+1)is None,
  abort_shared_candidate=psel(cd+[(W[0],('W',0))],1)is None,abort_equal_x=psel(cd+[(lo,('E',0))],1)is None,
  one_head_50_cells=len({x>>PH for x,_ in H})==1 and len({x>>PC for x,_ in H})==len(H),
  head_of_49_passes=mine(cd+H[:-1])[1:]==(0,HCAP)and psel(cd+H[:-1],n)==sel,head_of_50_aborts=mine(cd+H)[1]==1 and psel(cd+H,n)is None,
  close_pair_found_in_full_head=(H[0],Hc)in(psel(cd+H[:HCAP]+[Hc],n+1)or()),
  selected_real_oriented=sel is not None and all(real_x(a[0])and real_x(b[0])and 0<b[0]-a[0]<1<<PD and a[0]>>PC==b[0]>>PC for a,b in sel))
 res['ok']=all(v for k,v in res.items()if k!='pairs')
 return res
def kernel_xs(rg,pairs):
 # kernel X's: psel's choice from a shuffled pool of close pairs, decoys and a cell-split pair
 out=[]
 while len(out)<2*pairs:
  U=near_x(rg.randrange(1,P>>PC)<<PC|rg.randrange(1<<PC-1),0,1<<40,rg)
  d=(1+rg.randrange(1<<20),(1<<64)-1-rg.randrange(1<<20),rg.randrange(1,1<<64))[len(out)//2%3]
  V=real_x(U[0]+d)
  if V:out+=[U,(V[0],V[1]if rg.random()<0.5 else P-V[1])]
 B=rg.randrange(1,P>>PC)<<PC
 dec=[sc.mul(rg.randrange(1,ORD))for _ in range(8)]+[near_x(B,-(1<<40),0,rg),near_x(B,0,1<<40,rg)]
 pool=[(Z[0],Z)for Z in out+dec];rg.shuffle(pool)
 X=[Z for a,b in psel(pool,pairs)for Z in(a[1],b[1])]
 return X,sorted(X)==sorted(out)and psel(pool,pairs+1)is None
def gf128(a,b):
 z=0
 for i in range(7):
  if b>>i&1:z^=a<<i
 for i in range(13,6,-1):
  if z>>i&1:z^=0b10000011<<(i-7)
 return z
def class_pool(c):
 # class c=(u,w): lines b=u*a+v, v in 4w..4w+3, as ranks 128a+b
 u,w=c%128,c//128
 return[128*a+(gf128(u,a)^v)for v in range(4*w,4*w+4)for a in range(128)]
def design(rg):
 got={}
 for _ in range(12):
  c,d=rg.randrange(4096),rg.randrange(4096)
  got.setdefault('same'if c==d else'slope'if c%128==d%128 else'diff',set()).add(len(set(class_pool(c))&set(class_pool(d))))
 c=rg.randrange(4096);got.setdefault('slope',set()).add(len(set(class_pool(c))&set(class_pool((c+128)%4096))));got.setdefault('same',set()).add(len(set(class_pool(c))))
 return dict(shared_pairs=str({k:sorted(v)for k,v in got.items()}),ok=got.get('diff')=={16}and got['slope']=={0}and got['same']=={512})
def shared_x_family(rg):
 # proof.md 6 (M) d2: shared X's A, B, chain points U=+-V at two nodes: every pair through them collides at once
 A,B,U=(sc.mul(rg.randrange(1,ORD))for _ in range(3));hit=0
 for V in(U,sc.neg(U)):
  for X in(A,B):
   for s in(1,-1):
    L=sc.add(X,U if s==1 else sc.neg(U));hit+=any(sc.add(X,V if t==1 else sc.neg(V))==L for t in(1,-1))
 W=sc.mul(rg.randrange(1,ORD))
 return dict(simultaneous_collisions=hit,generic_point_collisions=sum(sc.add(A,U)[0]==sc.add(X,W)[0]for X in(A,B)),ok=hit==8)
def check(seed):
 rg=random.Random(seed)
 out=dict(constructed=constructed(rg),generator=generator(rg),design=design(rg),shared_x_family=shared_x_family(rg))
 out['all_checks']=all(v['ok']for v in out.values())
 return out
