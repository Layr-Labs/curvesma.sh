import sys,json,random,bisect,math,re,collections
from pathlib import Path
root=Path(__file__).resolve().parent
if(root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root));import a64
import range_witness,frac_witness as fw
witness=range_witness.check(root)
assert not any(re.search(r'\bret\s+\.',x)for x in(root/'p256_fields.s').read_text().splitlines())
import routing_check,cold_accounting_check
cold=cold_accounting_check.check()
routing=routing_check.check()
assert routing['formulas_ok']and routing['cands_ok']
assert routing['msgs_ok']and routing['slots_ok']and routing['tags_ok']and routing['chain_ok']and routing['stop_ok']and routing['all_paths']
assert routing['signaled']and routing['done_markers']==64*routing['groups_run']and routing['completions']==16*routing['groups_run']
isa=0
for shift in[0,1,35,63]:
 for av in[0,1,(1<<64)-1,1<<63]:
  for bv in[0,1,(1<<64)-1,1<<63]:
   pp,ll=a64.a64_load('check:\n orr x0,x1,x2,lsl #'+str(shift)+'\n ret\n');oo={}
   ct=a64.a64_run(pp,ll,'check',dict(x1=av,x2=bv,x30=777),{},out_regs=oo)
   assert ct==2 and oo['x0']==(av|((bv<<shift)&((1<<64)-1)));isa+=1
for lo,width in[(0,64),(35,1),(36,3),(63,1)]:
 for v in[0,1,(1<<64)-1,1<<63]:
  pp,ll=a64.a64_load('check:\n ubfx x0,x1,#'+str(lo)+',#'+str(width)+'\n ret\n');oo={}
  ct=a64.a64_run(pp,ll,'check',dict(x1=v,x30=777),{},out_regs=oo)
  assert ct==2 and oo['x0']==((v>>lo)&((1<<width)-1));isa+=1
for lo,width in[(0,35),(3,35),(36,28)]:
 for dv,sv in[(0,(1<<64)-1),((1<<64)-1,0),(1<<38,(1<<35)-1),(0xa5a5a5a5a5a5a5a5,0x5a5a5a5a5a5a5a5a)]:
  pp,ll=a64.a64_load('check:\n bfi x0,x1,#'+str(lo)+',#'+str(width)+'\n ret\n');oo={}
  ct=a64.a64_run(pp,ll,'check',dict(x0=dv,x1=sv,x30=777),{},out_regs=oo);mk=((1<<width)-1)<<lo
  assert ct==2 and oo['x0']==((dv&~mk)|((sv<<lo)&mk))&((1<<64)-1);isa+=1
M64=(1<<64)-1
for av in[0,1,(1<<64)-1,1<<63,0xa5a5a5a5a5a5a5a5]:
 for bv in[0,1,(1<<64)-1,1<<63,0x5a5a5a5a5a5a5a5a]:
  for prog_,want in[('add x0,x1,x2,lsr #28',(av+(bv>>28))&M64),('add x0,x1,x2,lsl #3',(av+(bv<<3))&M64),('orr x0,x1,x2,lsr #29',av|(bv>>29)),('extr x0,x1,x2,#35',((av<<64|bv)>>35)&M64)]:
   pp,ll=a64.a64_load('check:\n '+prog_+'\n ret\n');oo={}
   ct=a64.a64_run(pp,ll,'check',dict(x1=av,x2=bv,x30=777),{},out_regs=oo)
   assert ct==2 and oo['x0']==want,(prog_,hex(av),hex(bv));isa+=1
  pp,ll=a64.a64_load('check:\n tst x1,#0xffffe0001fffffff\n cset x0,eq\n ret\n');oo={}
  ct=a64.a64_run(pp,ll,'check',dict(x1=av^bv,x30=777),{},out_regs=oo)
  assert ct==3 and oo['x0']==int((av^bv)&0xffffe0001fffffff==0);isa+=1
for zb in[0x100000,0x100100]:
 mm={zb+8*k:5 for k in range(32)};mm[zb-8]=7;mm[zb+256]=9;oo={}
 ct=a64.a64_run(*a64.a64_load('check:\n dc zva, x1\n ret\n'),'check',dict(x1=zb,x30=777),mm,out_regs=oo)
 assert ct==2 and all(mm[zb+8*k]==0 for k in range(32))and mm[zb-8]==7 and mm[zb+256]==9 and oo['x1']==zb;isa+=1
LS=[(0x200000,0x200100)]
def ls(s,x,sr=LS):
 o={}
 try:a64.a64_run(*a64.a64_load('check:\n '+s+'\n ret\n'),'check',dict(x1=0x200000+x,x30=777),{a:a for a in range(0x200000,0x200100,8)},out_regs=o,ls64=(LS,sr));return o
 except AssertionError as e:return str(e)
r=ls('ld64b x4,[x1]',64);assert[r['x%d'%k]for k in range(4,12)]==list(range(0x200040,0x200080,8))and type(ls('st64b x4,[x1]',64))==dict
r=ls('ld64b x10,[x1]',128);assert[r['x%d'%k]for k in range(10,18)]==list(range(0x200080,0x2000c0,8));isa+=1
for s,x,sr in[('ld64b x4,[x1]',32,LS),('ld64b x5,[x1]',64,LS),('ld64b x24,[x1]',64,LS),('ld64b x4,[x1,#64]',0,LS),('ld64b x4,[x1]',256,LS),('st64b x4,[x1]',64,[])]:
 assert 'LS64' in ls(s,x,sr);isa+=1
assert 'DC ZVA' in ls('dc zva, x1',64)
import hashlib
assert hashlib.sha256('\n'.join(x for raw in a64.expand((root/'p256_inv.s').read_text()).splitlines()if(x:=raw.split('//')[0].rstrip()).strip()).encode()).hexdigest()=='d3b9a146237ea972367016b6b66079aa4a45c2708c28836e923f811ac53b093d'
base_prog,base_labels=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.base_generated_assembly(root));fbase_prog,fbase_labels=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.fused_base(root));unrows=[]
cfg=json.load(sys.stdin);rng=random.Random(cfg['seed'])
def regs(**kw):
 z={f'x{i}':rng.getrandbits(64) for i in range(31)};z.update(kw);return z
P=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
RR=1<<256;R=RR%P;RI=pow(R,-1,P);MASK=(1<<64)-1;M35=(1<<35)-1
G=(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)
def add(A,B):
 if A is None:return B
 if B is None:return A
 x,y=A;u,v=B
 if x==u:
  if(y+v)%P==0:return None
  m=(3*x*x-3)*pow(2*y,-1,P)%P
 else:m=(v-y)*pow(u-x,-1,P)%P
 z=(m*m-x-u)%P;return z,(m*(x-z)-y)%P
def norm(A):return(A[0],-A[1]%P)if A[1]*R%P>P//2 else A
def put(mem,a,v):
 for k in range(4):mem[a+8*k]=(v>>(64*k))&MASK
def get(mem,a):return sum(mem[a+8*k]<<(64*k)for k in range(4))
spg,slb=a64.a64_load((root/'p256_fields.s').read_text());sr=random.Random(3901)
edges={0,1,2,P-2,P-1,P,P+1,RR-2,RR-1}
for bit in[32,64,96,128,160,192,224,255]:
 for off in[-1,0,1]:edges.add(((1<<bit)+off)%RR)
pairs=[(a,b)for a in sorted(edges)for b in sorted(edges)]+[(sr.randrange(RR),sr.randrange(RR))for _ in range(512)]
short_counts={};short_alias=[0,0,0];ni=-pow(P,-1,RR)%RR
for case,(a,b)in enumerate(pairs):
 mode=case%3;short_alias[mode]+=1;z=[0x11000,0x22000,0x33000][mode];mem={}
 for addr in[0x11000,0x22000,0x33000]:
  for j in range(-2,6):mem[addr+8*j]=sr.getrandbits(64)
 put(mem,0x11000,a);put(mem,0x22000,b);before=dict(mem)
 rg={f'x{i}':sr.getrandbits(64) for i in range(31)};rg.update(x0=z,x1=0x11000,x2=0x22000,sp=0x44000,x30=0x77777777);o={}
 n=a64.a64_run(spg,slb,'k_mul',rg,mem,out_regs=o,limit=500)
 t=(a*b+(a*b*ni%RR)*P)//RR;assert 0<=t<RR+P
 assert get(mem,z)==(t if t<RR else t-P)and get(mem,z)%P==a*b*RI%P
 assert n==(104 if t<RR else 110)
 assert all(mem[k]==v for k,v in before.items()if not z<=k<z+32)
 assert all(o[f'x{i}']==rg[f'x{i}'] for i in range(18,31)) and o['sp']==rg['sp']
 short_counts[n]=short_counts.get(n,0)+1
assert set(short_counts)=={104,110}
short=dict(cases=len(pairs),counts=short_counts,aliases=short_alias,full_width=True,canaries_and_preserved_registers=True)
# section16 tail model
def tm(a,b):
 t=(a*b+(a*b*ni%RR)*P)>>256;return(t-P,1)if t>=RR else(t,0)
def tfwd(X,its):
 A=[R];d=[];f=[]
 for x,y,_ in its:d.append(abs(x-X[0]));a,c=tm(A[-1],d[-1]);A.append(a);f.append(c)
 return f,A,d
def tbwd(X,its,A,d,B):
 e=[0]*len(its)
 for j in range(len(its)-1,-1,-1):
  x,y,_=its[j];i,c0=tm(B,A[j]);B,c1=tm(B,d[j]);w,c2=tm(i,i);_,c3=tm(y,w);e[j]=(c0,c1,c2,c3)
 return e
tres={};tinc=0
def build(items):
 seq=sorted(items,key=lambda a:(a[1],a[2]));nodes=[]
 def rec(lo,hi):
  node={'lo':lo,'hi':hi,'items':sorted(seq[lo:hi],key=lambda a:(a[0],a[2])),'base':0x100000000+len(nodes)*0x10000};nodes.append(node)
  if hi-lo>1:
   m=(lo+hi)//2;node['children']=(rec(lo,m),rec(m,hi))
  return node
 return rec(0,len(seq)),seq,nodes
def descriptors(tree,seq,xx,yy):
 rank=bisect.bisect_left([a[1]for a in seq],yy);out=[];cover=[]
 def visit(nd):
  if nd['hi']<=rank:cover.append((nd,1));return
  if nd['lo']>=rank:cover.append((nd,0));return
  for ch in nd['children']:visit(ch)
 visit(tree)
 for nd,yl in cover:
  xs=[a[0]for a in nd['items']]
  cuts=sorted(set([0,len(xs),bisect.bisect_left(xs,xx),len(xs)if xx==0 else bisect.bisect_left(xs,RR-xx)]))
  for a,b in zip(cuts,cuts[1:]):
   if a==b:continue
   dxlo=xs[a]<xx;carry=xs[a]+xx>=RR
   items=nd['items'][a:b]
   assert all((v[0]<xx)==dxlo and(v[0]+xx>=RR)==carry and(v[1]<yy)==bool(yl)for v in items)
   out.append((nd,a,b,int(dxlo),yl,int(carry)))
 assert len(cover)<=math.ceil(math.log2(len(seq)))+1
 assert len(out)<=3*len(cover)
 idck=[v[2]for nd,a,b,_,_,_ in out for v in nd['items'][a:b]]
 assert sorted(idck)==sorted(v[2]for v in seq)
 return out
combin=0
for n in[1,2,4,8,16,32,64]:
 items=[(rng.randrange(P),rng.randrange(1,P//2+1),i)for i in range(n)];tr,seq,_=build(items)
 for xx in[0,1,P-1,P//2]+[rng.randrange(P)for _ in range(20)]:
  for yy in[0,1,P//2,seq[n//2][1]]:
   descriptors(tr,seq,xx,yy);combin+=1
bounds=0
for n in[1,2,3,5,7,9,17,31,33,65]:
 for xx in[0,1,RR-P,RR-P+1,RR//2,P//2,P-1]:
  pool=[0,1,xx,P-1]
  if xx and RR-xx<P:pool+=[RR-xx-1,RR-xx,RR-xx+1]
  pool=[x for x in pool if 0<=x<P]
  items=[(pool[i%len(pool)],1+(i%3),i)for i in range(n)]
  tr,seq,_=build(items)
  for yy in[0,1,2,3,4,P//2]:
   ds=descriptors(tr,seq,xx,yy)
   high=sum(b-a for nd,a,b,dx,yl,c in ds if c)
   assert high==sum(x+xx>=RR for x,y,l in items)
   bounds+=1
prog,labels=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.generated_assembly(root)+(root/'p256_inv.s').read_text())
aborts=0
for entry,delta in[('pair_fwd_lo',-1),('pair_fwd_hi',1)]:
 for zero_at in[0,1,3,7]:
  m={};yp=0x11000;bp=0x22000;xp=0x33000;scratch=0x44000;sp=0x55000
  put(m,xp,P//2);put(m,xp+32,2);put(m,bp-32,R)
  for j in range(zero_at+1):
   put(m,yp+72*j,P//2 if j==zero_at else P//2+delta)
   put(m,yp+72*j+32,1);m[yp+72*j+64]=j<<36
  source=dict(m);out={};sentinel=0x123456789abcdef
  count=a64.a64_run(prog,labels,entry,regs(x0=zero_at+1,x1=yp,x2=bp,x3=xp,x4=scratch,sp=sp,x30=sentinel),m,limit=200000,out_regs=out)
  assert out['x0']==1 and out['sp']==sp and out['x30']==sentinel
  assert count<=441*(zero_at+1)//4+46
  assert all(m[k]==v for k,v in source.items()if yp<=k<yp+72*(zero_at+1)or xp<=k<xp+64)
  assert all((bp-32<=k<bp+64*(zero_at+1))or(scratch<=k<scratch+256)or(sp-16<=k<sp)or k in source for k in m)
  assert all(isinstance(v,int)and 0<=v<(1<<64)for v in m.values())
  n=zero_at+1
  ic=a64.a64_run(prog,labels,'k_inv',regs(x0=scratch,x1=bp+64*(n-1)+32,x2=0x80000,sp=sp,x30=sentinel),m,limit=200000)
  assert ic<=39479
  FBA=0x7000000000;XAA=0x50000;fw.install(m,FBA,XAA,P//2,2,[(P//2+delta,1,j)for j in range(zero_at+1)])
  m0=dict(m)
  WBASE=1<<38;RBASE=0x2000000;SEND=0x3000000;SPORT=0x90000;W2=17|(31<<35);out2={}
  ct=a64.a64_run(prog,labels,'pair_bwd_f',regs(x0=n,x1=yp+72*n+64,x2=bp+64*(n-1),x3=XAA,x4=scratch,x5=WBASE,x6=W2,x7=RBASE,x8=SEND+65536,x9=SPORT,x10=FBA,sp=sp,x30=sentinel),collections.defaultdict(int,m0),limit=600000,out_regs=out2,ls64=fw.ls64_regions(FBA,XAA))
  assert ct<=2173*n//4+34 and out2['sp']==sp and out2['x30']==sentinel
  aborts+=1
limbs=0
for yx,xx in[(1,1<<128),(1<<128,1),(0,P-1),(P-1,0),(1,(1<<128)+1),((1<<192)-1,1<<192),(0,1),(1,0),((1<<128)-1,1<<128),(1<<128,(1<<128)-1)]:
 m={};yp=0x11000;bp=0x22000;xp=0x33000;scratch=0x44000;sp=0x55000
 put(m,yp,yx);put(m,yp+32,1);m[yp+64]=0;put(m,xp,xx);put(m,xp+32,2);put(m,bp-32,R)
 entry='pair_fwd_lo' if yx<xx else 'pair_fwd_hi';out={}
 a64.a64_run(prog,labels,entry,regs(x0=1,x1=yp,x2=bp,x3=xp,x4=scratch,sp=sp,x30=999),m,limit=200000,out_regs=out)
 assert out['x0']==1 and get(m,bp)==abs(xx-yx),(hex(yx),hex(xx),hex(get(m,bp)));limbs+=1
assert limbs==10
import native_tree
nb=native_tree.Builder(root)
points=[];A=G
for i in range(300):points.append(A);A=add(A,G)
rows=[];variants=set();coupl=0;ids=0;lerr=[]
for n in[1,2,4,8,16,32,64,128]:
 for xi in[0,11,250]:
  X0=points[xi];X=norm(X0);original=[Y for Y in points if Y[0]!=X[0]][:n];ys=[norm(Y)for Y in original]
  id_values=([0,(1<<27)-1,1<<27,(1<<28)-1]+list(range(1,n+1)))[:n];rng.shuffle(id_values)
  tr,seq,nodes=nb.build([(Y[0]*R%P,Y[1]*R%P,id_values[l])for l,Y in enumerate(ys)]);ds=descriptors(tr,seq,X[0]*R%P,X[1]*R%P)
  mem={};xp=0x30000;bp=0x40000;sptr=0x60000;sp=0x70000
  put(mem,xp,X[0]*R%P);put(mem,xp+32,X[1]*R%P);put(mem,bp-32,R)
  for nd in nodes:
   for j,(x,y,l)in enumerate(nd['items']):
    a=nd['base']+72*j;put(mem,a,x);put(mem,a+32,y);mem[a+64]=l<<36
  order=[];off=0;segments=[];fc=bc=kc=0
  tits=[v for nd,a,b,*_ in ds for v in nd['items'][a:b]];Xm=(X[0]*R%P,X[1]*R%P);tfl,tA,td=tfwd(Xm,tits)
  for nd,a,b,dxlo,yl,carry in ds:
   cnt=b-a;ptr=nd['base']+72*a;out={};entry='pair_fwd_lo' if dxlo else 'pair_fwd_hi'
   fmem=dict(mem);ct=a64.a64_run(prog,labels,entry,regs(x0=cnt,x1=ptr,x2=bp+64*off,x3=xp,x4=sptr,sp=sp,x30=12345),mem,limit=200000,out_regs=out)
   foct=a64.a64_run(base_prog,base_labels,entry,regs(x0=cnt,x1=ptr,x2=bp+64*off,x3=xp,x4=sptr,sp=sp,x30=12345),fmem,limit=200000)
   assert mem==fmem and ct==foct+6-30*(cnt//8),(cnt,ct,foct)
   unrows.append(dict(n=cnt,new=ct,old=foct,fwd=1))
   assert out['x0']==1 and out['sp']==sp and out['x30']==12345 and ct<=441*cnt//4+46;fc+=ct
   tt=ct-4*sum(tfl[off:off+cnt]);assert 4*tt<=425*cnt+184 and get(mem,bp+64*(off+cnt-1)+32)==tA[off+cnt];tres.setdefault((entry,cnt),set()).add(tt)
   segments.append((off,cnt,ptr,yl,carry));order.extend(nd['items'][a:b]);off+=cnt
  assert off==n
  ic=a64.a64_run(prog,labels,'k_inv',regs(x0=sptr,x1=bp+64*(n-1)+32,x2=0x80000,sp=sp,x30=12345),mem,limit=200000)
  tinc+=n
  expkeys=[]
  for j,(_,_,l)in enumerate(order):
   Y=ys[id_values.index(l)]
   expkeys.append((add(X,Y)[0]*R%P,add(X,(Y[0],-Y[1]%P))[0]*R%P,l))
   assert{expkeys[-1][0]*RI%P,expkeys[-1][1]*RI%P}=={add(X0,original[id_values.index(l)])[0],add(X0,(original[id_values.index(l)][0],-original[id_values.index(l)][1]%P))[0]};coupl+=1
  FBA=0x7000000000;XAA=0x50000;fw.install(mem,FBA,XAA,X[0],X[1],tits)
  W2=17|(31<<35);exp,tl,er=fw.expected(get(mem,sptr),tA,td,tits,X[0],X[1],W2);lerr+=er
  WBASE=1<<38;RBASE=0x2000000;SEND=0x3000000;SPORT=0x90000;win_of={}
  for rec,d1 in exp:
   if d1 not in win_of:win_of[d1]=0x8000000+len(win_of)*8192;mem[WBASE+d1*8]=win_of[d1]+32
  for off,cnt,ptr,yl,carry in reversed(segments):
   call_regs=regs(x0=cnt,x1=ptr+72*cnt+64,x2=bp+64*(off+cnt-1),x3=XAA,x4=sptr,x5=WBASE,x6=W2,x7=RBASE,x8=SEND+65536,x9=SPORT,x10=FBA,sp=sp,x30=12345)
   out={};ct=a64.a64_run(prog,labels,'pair_bwd_f',call_regs,mem,limit=600000,out_regs=out,ls64=fw.ls64_regions(FBA,XAA))
   tt=ct-4*sum(tl[off:off+cnt]);assert tt==fw.tail_free(cnt),(cnt,ct,tt)
   unrows.append(dict(n=cnt,new=ct,tail_free=tt,bwd=1))
   assert ct<=2173*cnt//4+34 and out['sp']==sp and out['x30']==12345;bc+=ct;variants.add((yl,carry));tres.setdefault(('b',cnt),set()).add(tt)
  assert get(mem,sptr)==R
  bywin={}
  for rec,d1 in exp:bywin.setdefault(d1,[]).append(rec)
  for d1,recs in bywin.items():
   win=win_of[d1];assert mem[WBASE+d1*8]==win+32+32*len(recs)
   for t,rec in enumerate(recs):
    assert tuple(mem[win+32+32*t+8*q]for q in range(4))==rec;ids+=1
  high=sum(c*carry for _,c,_,_,carry in segments)
  rows.append(dict(n=n,descriptors=len(ds),forward=fc,backward=bc,keys=kc,high=high,inverse=ic))
assert variants=={(0,0),(0,1),(1,0),(1,1)}and all(len(v)==1 for v in tres.values())and min(lerr)>=-17 and max(lerr)<=9
BC=0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b;cr=random.Random(cfg['seed']*7+1)
def tpt():
 while 1:
  x=cr.randrange(P);w=(x*x*x-3*x+BC)%P;y=pow(w,(P+1)//4,P)
  if y*y%P==w:y=y*R%P;return x*R%P,min(y,P-y),0
cs=[0]*5;cn=0;tmax=0
for t in range(2):
 tys=[tpt()for _ in range(1024)];ts=[0]*5
 for _ in range(8):
  X=tpt();its=sorted(tys,key=lambda v:(v[1]<X[1],v[0]<X[0],v[0]));f,A,d=tfwd(X,its);e=tbwd(X,its,A,d,pow(A[-1]*RI%P,-1,P)*R%P)
  for j in range(1024):
   for z,v in enumerate((f[j],)+e[j]):ts[z]+=v>0;cs[z]+=v>0
  cn+=1024
 tmax=max(tmax,max(ts)/8192)
crate=[c/cn for c in cs];tmean=4*sum(crate)
assert tmax<0.5 and cn==16384
tails=dict(kernel_call_shapes=len(tres),kernel_incidences=tinc,kernel_model_residual_single_valued=True,census_incidences=cn,site_rates=[round(v,5)for v in crate],max_table_site_rate=round(tmax,5),extra_upper_mean=round(tmean,4),charged=10,worst=20)
import merge_check
merge_audit=merge_check.check(root)
import search_check
search_audit=search_check.check(root)
import raw_differential
raw_audit=raw_differential.check(root,build,descriptors)
fw_audit=raw_differential.full_window(root,build,descriptors)
import runpy,contextlib,io
buf=io.StringIO()
with contextlib.redirect_stdout(buf):runpy.run_path(str(root/'range_ledger.py'),run_name='__main__')
ledger=json.loads(buf.getvalue())
import fixed_layout
layout=fixed_layout.check(ledger)
saved_stdin=sys.stdin;sys.stdin=io.StringIO(json.dumps(dict(seed=193)));resolution_buf=io.StringIO()
with contextlib.redirect_stdout(resolution_buf):runpy.run_path(str(root/'resolution_check.py'),run_name='__main__')
sys.stdin=saved_stdin;resolution=json.loads(resolution_buf.getvalue());assert resolution['all_checks']
print(json.dumps(dict(unroll=unrows,short_mul=short,research_only=True,cold_accounting=cold,routing=routing,fixed_layout=layout,isa_checks=isa,resolution=resolution,native_tree=dict(merge_calls=nb.calls,instructions=nb.instructions,records=nb.records,scope='native merges integrated into actual-point arithmetic; caller control not counted here'),merge_helper=merge_audit,search_helper=search_audit,raw_diff=raw_audit,full_window=fw_audit,witness=witness,sparse_ID_checks=True,ledger=ledger,boundary=bounds,zero_denom_garbage=aborts,limb_spec_cases=limbs,random_GPR_init=True,combin_cases=combin,point_couplings=coupl,orig_id_emit=ids,descriptor_kinds=sorted(variants),label_error=[min(lerr),max(lerr)],rows=rows,tails=tails,all_checks=True,scope='Finite evidence only; the analytical argument is in proof.md.')))
