import sys,json,random,bisect,math,re
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root));import a64
import range_witness
witness=range_witness.check(root)
assert not any(re.search(r'\bret\s+\.', x) for x in (root/'p256_fields.s').read_text().splitlines())
import routing_check,scanner_check,cold_accounting_check
cold=cold_accounting_check.check()
routing=routing_check.check();scanner=scanner_check.check(root)
assert routing['formulas_ok'] and routing['candidates_equal_reference_set']
assert routing['messages_equal_emit_reference'] and routing['slots_equal_reference'] and routing['source_tags_ok'] and routing['buckets_pointer_order_and_scan_ok']
isa_cases=0
for shift in [0,1,35,63]:
 for av in [0,1,(1<<64)-1,1<<63]:
  for bv in [0,1,(1<<64)-1,1<<63]:
   pp,ll=a64.a64_load('check:\n orr x0,x1,x2,lsl #'+str(shift)+'\n ret\n');oo={}
   ct=a64.a64_run(pp,ll,'check',dict(x1=av,x2=bv,x30=777),{},out_regs=oo)
   assert ct==2 and oo['x0']==(av|((bv<<shift)&((1<<64)-1)));isa_cases+=1
for lo,width in [(0,64),(35,1),(36,3),(63,1)]:
 for v in [0,1,(1<<64)-1,1<<63]:
  pp,ll=a64.a64_load('check:\n ubfx x0,x1,#'+str(lo)+',#'+str(width)+'\n ret\n');oo={}
  ct=a64.a64_run(pp,ll,'check',dict(x1=v,x30=777),{},out_regs=oo)
  assert ct==2 and oo['x0']==((v>>lo)&((1<<width)-1));isa_cases+=1
import hashlib
assert hashlib.sha256('\n'.join(x for raw in a64.expand((root/'p256_inv.s').read_text()).splitlines() if (x:=raw.split('//')[0].rstrip()).strip()).encode()).hexdigest()=='d3b9a146237ea972367016b6b66079aa4a45c2708c28836e923f811ac53b093d'
assert hashlib.sha256('\n'.join(x for raw in a64.expand((root/'pair_inline_reference.s').read_text()).splitlines() if (x:=raw.split('//')[0].rstrip()).strip()).encode()).hexdigest()=='4f44679819e9d1e50306d4f31ce86604c999f0bb599871599d8f8cbdb8f80787'
base_prog,base_labels=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.base_generated_assembly(root));fbase_prog,fbase_labels=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.fused_base(root));unroll_rows=[]
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
  if (y+v)%P==0:return None
  m=(3*x*x-3)*pow(2*y,-1,P)%P
 else:m=(v-y)*pow(u-x,-1,P)%P
 z=(m*m-x-u)%P;return z,(m*(x-z)-y)%P
def norm(A):return (A[0],-A[1]%P) if A[1]*R%P>P//2 else A
def put(mem,a,v):
 for k in range(4):mem[a+8*k]=(v>>(64*k))&MASK
def get(mem,a):return sum(mem[a+8*k]<<(64*k) for k in range(4))
spg,slb=a64.a64_load((root/'p256_fields.s').read_text());sr=random.Random(3901)
edges={0,1,2,P-2,P-1,P,P+1,RR-2,RR-1}
for bit in [32,64,96,128,160,192,224,255]:
 for off in [-1,0,1]:edges.add(((1<<bit)+off)%RR)
pairs=[(a,b) for a in sorted(edges) for b in sorted(edges)]+[(sr.randrange(RR),sr.randrange(RR)) for _ in range(512)]
short_counts={};short_alias=[0,0,0];ni=-pow(P,-1,RR)%RR
for case,(a,b) in enumerate(pairs):
 mode=case%3;short_alias[mode]+=1;z=[0x11000,0x22000,0x33000][mode];mem={}
 for addr in [0x11000,0x22000,0x33000]:
  for j in range(-2,6):mem[addr+8*j]=sr.getrandbits(64)
 put(mem,0x11000,a);put(mem,0x22000,b);before=dict(mem)
 rg={f'x{i}':sr.getrandbits(64) for i in range(31)};rg.update(x0=z,x1=0x11000,x2=0x22000,sp=0x44000,x30=0x77777777);o={}
 n=a64.a64_run(spg,slb,'k_mul',rg,mem,out_regs=o,limit=500)
 t=(a*b+(a*b*ni%RR)*P)//RR;assert 0<=t<RR+P
 assert get(mem,z)==(t if t<RR else t-P) and get(mem,z)%P==a*b*RI%P
 assert n==(104 if t<RR else 110)
 assert all(mem[k]==v for k,v in before.items() if not z<=k<z+32)
 assert all(o[f'x{i}']==rg[f'x{i}'] for i in range(18,31)) and o['sp']==rg['sp']
 short_counts[n]=short_counts.get(n,0)+1
assert set(short_counts)=={104,110}
short_check=dict(cases=len(pairs),counts=short_counts,aliases=short_alias,full_width=True,canaries_and_preserved_registers=True)
def build(items):
 seq=sorted(items,key=lambda a:(a[1],a[2]));nodes=[]
 def rec(lo,hi):
  node={'lo':lo,'hi':hi,'items':sorted(seq[lo:hi],key=lambda a:(a[0],a[2])),'base':0x100000000+len(nodes)*0x10000};nodes.append(node)
  if hi-lo>1:
   m=(lo+hi)//2;node['children']=(rec(lo,m),rec(m,hi))
  return node
 return rec(0,len(seq)),seq,nodes
def descriptors(tree,seq,xx,yy):
 rank=bisect.bisect_left([a[1] for a in seq],yy);out=[];cover=[]
 def visit(nd):
  if nd['hi']<=rank:cover.append((nd,1));return
  if nd['lo']>=rank:cover.append((nd,0));return
  for ch in nd['children']:visit(ch)
 visit(tree)
 for nd,yl in cover:
  xs=[a[0] for a in nd['items']]
  cuts=sorted(set([0,len(xs),bisect.bisect_left(xs,xx),len(xs) if xx==0 else bisect.bisect_left(xs,RR-xx)]))
  for a,b in zip(cuts,cuts[1:]):
   if a==b:continue
   dxlo=xs[a]<xx;carry=xs[a]+xx>=RR
   items=nd['items'][a:b]
   assert all((v[0]<xx)==dxlo and (v[0]+xx>=RR)==carry and (v[1]<yy)==bool(yl) for v in items)
   out.append((nd,a,b,int(dxlo),yl,int(carry)))
 assert len(cover)<=math.ceil(math.log2(len(seq)))+1
 assert len(out)<=3*len(cover)
 ids=[v[2] for nd,a,b,_,_,_ in out for v in nd['items'][a:b]]
 assert sorted(ids)==sorted(v[2] for v in seq)
 return out
combinatorial=0
for n in [1,2,4,8,16,32,64]:
 items=[(rng.randrange(P),rng.randrange(1,P//2+1),i) for i in range(n)];tr,seq,_=build(items)
 for xx in [0,1,P-1,P//2]+[rng.randrange(P) for _ in range(20)]:
  for yy in [0,1,P//2,seq[n//2][1]]:
   descriptors(tr,seq,xx,yy);combinatorial+=1
boundary_cases=0
for n in [1,2,3,5,7,9,17,31,33,65]:
 for xx in [0,1,RR-P,RR-P+1,RR//2,P//2,P-1]:
  pool=[0,1,xx,P-1]
  if xx and RR-xx<P:pool += [RR-xx-1,RR-xx,RR-xx+1]
  pool=[x for x in pool if 0<=x<P]
  items=[(pool[i%len(pool)],1+(i%3),i) for i in range(n)]
  tr,seq,_=build(items)
  for yy in [0,1,2,3,4,P//2]:
   ds=descriptors(tr,seq,xx,yy)
   high=sum(b-a for nd,a,b,dx,yl,c in ds if c)
   assert high==sum(x+xx>=RR for x,y,l in items)
   boundary_cases+=1
prog,labels=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.generated_assembly(root)+(root/'p256_inv.s').read_text())
abort_cases=0
for entry,delta in [('pair_fwd_lo',-1),('pair_fwd_hi',1)]:
 for zero_at in [0,1,3,7]:
  m={};yp=0x11000;bp=0x22000;xp=0x33000;scratch=0x44000;sp=0x55000
  put(m,xp,P//2);put(m,xp+32,2);put(m,bp-32,R)
  for j in range(zero_at+1):
   put(m,yp+72*j,P//2 if j==zero_at else P//2+delta)
   put(m,yp+72*j+32,1);m[yp+72*j+64]=j<<36
  source=dict(m);out={};sentinel=0x123456789abcdef
  count=a64.a64_run(prog,labels,entry,regs(x0=zero_at+1,x1=yp,x2=bp,x3=xp,x4=scratch,sp=sp,x30=sentinel),m,limit=200000,out_regs=out)
  assert out['x0']==1 and out['sp']==sp and out['x30']==sentinel
  assert count<=449*(zero_at+1)//4+30
  assert all(m[k]==v for k,v in source.items() if yp<=k<yp+72*(zero_at+1) or xp<=k<xp+64)
  assert all((bp-32<=k<bp+64*(zero_at+1)) or (scratch<=k<scratch+256) or (sp-16<=k<sp) or k in source for k in m)
  assert all(isinstance(v,int) and 0<=v<(1<<64) for v in m.values())
  n=zero_at+1
  ic=a64.a64_run(prog,labels,'k_inv',regs(x0=scratch,x1=bp+64*(n-1)+32,x2=0x80000,sp=sp,x30=sentinel),m,limit=200000)
  assert ic<=39479
  m0=dict(m)
  ct0=a64.a64_run(base_prog,base_labels,'pair_bwd_y0_c0',regs(x0=n,x1=yp+72*(n-1),x2=bp+64*(n-1),x3=xp,x4=scratch,sp=sp,x30=sentinel),m,limit=300000)
  assert ct0<=628*n+11
  WBASE=0x1000000;RBASE=0x2000000;SEND=0x3000000;SPORT=0x90000;W2=17|(31<<35);win_of={}
  for j in range(2*n):
   d1=get(m,bp+32*j)&M35
   if d1 not in win_of:win_of[d1]=0x8000000+len(win_of)*4096;m0[WBASE+d1*8]=win_of[d1]+32
  out2={}
  ct=a64.a64_run(prog,labels,'pair_bwd_y0_c0',regs(x0=n,x1=yp+72*(n-1),x2=bp+64*(n-1),x3=xp,x4=scratch,x5=WBASE,x6=W2,x7=RBASE,x8=SEND+65536,x9=SPORT,sp=sp,x30=sentinel),m0,limit=600000,out_regs=out2)
  assert ct<=2533*n//4+33 and out2['sp']==sp and out2['x30']==sentinel
  assert all((bp-32<=k<bp+64*n) or (scratch<=k<scratch+256) or (sp-64<=k<sp+64) or k in source or (0x80000<=k<0x80000+4096) or (WBASE<=k<WBASE+8*(1<<35)) or (0x8000000<=k<0x8000000+2*n*4096) for k in m0)
  abort_cases+=1
limb_cases=0
for yx,xx in [(1,1<<128),(1<<128,1),(0,P-1),(P-1,0),(1,(1<<128)+1),((1<<192)-1,1<<192),(0,1),(1,0),((1<<128)-1,1<<128),(1<<128,(1<<128)-1)]:
 m={};yp=0x11000;bp=0x22000;xp=0x33000;scratch=0x44000;sp=0x55000
 put(m,yp,yx);put(m,yp+32,1);m[yp+64]=0;put(m,xp,xx);put(m,xp+32,2);put(m,bp-32,R)
 entry='pair_fwd_lo' if yx<xx else 'pair_fwd_hi';out={}
 a64.a64_run(prog,labels,entry,regs(x0=1,x1=yp,x2=bp,x3=xp,x4=scratch,sp=sp,x30=999),m,limit=200000,out_regs=out)
 assert out['x0']==1 and get(m,bp)==abs(xx-yx),(hex(yx),hex(xx),hex(get(m,bp)));limb_cases+=1
assert limb_cases==10
import native_tree
native_builder=native_tree.Builder(root)
points=[];A=G
for i in range(300):points.append(A);A=add(A,G)
rows=[];variants=set();couplings=0;idchecks=0
for n in [1,2,4,8,16,32,64,128]:
 for xi in [0,11,250]:
  X0=points[xi];X=norm(X0);original=[Y for Y in points if Y[0]!=X[0]][:n];ys=[norm(Y) for Y in original]
  id_values=([0,(1<<27)-1,1<<27,(1<<28)-1]+list(range(1,n+1)))[:n];rng.shuffle(id_values)
  tr,seq,nodes=native_builder.build([(Y[0]*R%P,Y[1]*R%P,id_values[l]) for l,Y in enumerate(ys)]);ds=descriptors(tr,seq,X[0]*R%P,X[1]*R%P)
  mem={};xp=0x30000;bp=0x40000;sptr=0x60000;sp=0x70000
  put(mem,xp,X[0]*R%P);put(mem,xp+32,X[1]*R%P);put(mem,bp-32,R)
  for nd in nodes:
   for j,(x,y,l) in enumerate(nd['items']):
    a=nd['base']+72*j;put(mem,a,x);put(mem,a+32,y);mem[a+64]=l<<36
  order=[];off=0;segments=[];fc=bc=kc=0
  for nd,a,b,dxlo,yl,carry in ds:
   cnt=b-a;ptr=nd['base']+72*a;out={};entry='pair_fwd_lo' if dxlo else 'pair_fwd_hi'
   fmem=dict(mem);ct=a64.a64_run(prog,labels,entry,regs(x0=cnt,x1=ptr,x2=bp+64*off,x3=xp,x4=sptr,sp=sp,x30=12345),mem,limit=200000,out_regs=out)
   foct=a64.a64_run(base_prog,base_labels,entry,regs(x0=cnt,x1=ptr,x2=bp+64*off,x3=xp,x4=sptr,sp=sp,x30=12345),fmem,limit=200000)
   assert mem==fmem and ct==foct+4-14*(cnt//8),(cnt,ct,foct)
   unroll_rows.append(dict(n=cnt,new=ct,old=foct,fwd=1))
   assert out['x0']==1 and out['sp']==sp and out['x30']==12345 and ct<=449*cnt//4+30;fc+=ct
   segments.append((off,cnt,ptr,yl,carry));order.extend(nd['items'][a:b]);off+=cnt
  assert off==n
  ic=a64.a64_run(prog,labels,'k_inv',regs(x0=sptr,x1=bp+64*(n-1)+32,x2=0x80000,sp=sp,x30=12345),mem,limit=200000)
  expkeys=[]
  for j,(_,_,l) in enumerate(order):
   Y=ys[id_values.index(l)]
   expkeys.append((add(X,Y)[0]*R%P,add(X,(Y[0],-Y[1]%P))[0]*R%P,l))
   assert {expkeys[-1][0]*RI%P,expkeys[-1][1]*RI%P}=={add(X0,original[id_values.index(l)])[0],add(X0,(original[id_values.index(l)][0],-original[id_values.index(l)][1]%P))[0]};couplings+=1
  WBASE=0x1000000;RBASE=0x2000000;SEND=0x3000000;SPORT=0x90000;W2=17|(31<<35);win_of={}
  for (kp,km,l) in expkeys:
   for key in (kp,km):
    d1=key&M35
    if d1 not in win_of:win_of[d1]=0x8000000+len(win_of)*4096;mem[WBASE+d1*8]=win_of[d1]+32
  for off,cnt,ptr,yl,carry in reversed(segments):
   call_regs=regs(x0=cnt,x1=ptr+72*(cnt-1),x2=bp+64*(off+cnt-1),x3=xp,x4=sptr,x5=WBASE,x6=W2,x7=RBASE,x8=SEND+65536,x9=SPORT,sp=sp,x30=12345)
   fmem2=dict(mem);out={};ct=a64.a64_run(prog,labels,f'pair_bwd_y{yl}_c{carry}',call_regs,mem,limit=600000,out_regs=out)
   foct=a64.a64_run(fbase_prog,fbase_labels,f'pair_bwd_y{yl}_c{carry}',call_regs,fmem2,limit=600000)
   assert mem==fmem2 and ct==foct+4-14*(cnt//8),(cnt,ct,foct)
   unroll_rows.append(dict(n=cnt,new=ct,old=foct,fused=1))
   assert ct<=2533*cnt//4+33 and out['sp']==sp and out['x30']==12345;bc+=ct;variants.add((yl,carry))
  assert get(mem,sptr)==R
  emitted=[]
  for off,cnt,ptr,yl,carry in reversed(segments):
   for j in range(off+cnt-1,off-1,-1):
    kp,km,l=expkeys[j];emitted.append((kp,0,l));emitted.append((km,1,l))
  assert len(emitted)==2*n;bywin={}
  for (key,sign,l) in emitted:bywin.setdefault(key&M35,[]).append((key,sign,l))
  for d1,recs in bywin.items():
   win=win_of[d1]
   assert mem[WBASE+d1*8]==win+32+32*len(recs)
   for t,(key,sign,l) in enumerate(recs):
    base=win+32+32*t
    assert mem[base]==((key>>128)&MASK) and mem[base+8]==(key>>192)
    assert mem[base+16]==W2 and mem[base+24]==((key>>64)&M35)|(sign<<35)|(l<<36);idchecks+=1
  high=sum(c*carry for _,c,_,_,carry in segments)
  rows.append(dict(n=n,descriptors=len(ds),forward=fc,backward=bc,keys=kc,high=high,inverse=ic))
assert variants=={(0,0),(0,1),(1,0),(1,1)}
import merge_check
merge_audit=merge_check.check(root)
import search_check
search_audit=search_check.check(root)
import raw_differential
raw_audit=raw_differential.check(root,build,descriptors)
import fused_check
fused_audit=fused_check.check(root)
import runpy,contextlib,io
buf=io.StringIO()
with contextlib.redirect_stdout(buf):runpy.run_path(str(root/'range_ledger.py'),run_name='__main__')
ledger=json.loads(buf.getvalue())
import fixed_layout
layout=fixed_layout.check(ledger)
saved_stdin=sys.stdin;sys.stdin=io.StringIO(json.dumps(dict(seed=193)));resolution_buf=io.StringIO()
with contextlib.redirect_stdout(resolution_buf):runpy.run_path(str(root/'resolution_check.py'),run_name='__main__')
sys.stdin=saved_stdin;resolution=json.loads(resolution_buf.getvalue());assert resolution['all_checks']
print(json.dumps(dict(unroll_counts=unroll_rows,standalone_short_mul=short_check,research_only=True,cold_accounting=cold,routing=routing,scanner=scanner,fixed_layout=layout,isa_extension_checks=isa_cases,resolution=resolution,native_tree=dict(merge_calls=native_builder.calls,instructions=native_builder.instructions,records=native_builder.records,scope='native merges integrated into actual-point arithmetic; caller control not counted here'),merge_helper=merge_audit,search_helper=search_audit,raw_differential=raw_audit,fused_check=fused_audit,witness=witness,sparse_ID_boundary_checks=True,ledger=ledger,boundary_cases=boundary_cases,zero_denominator_garbage_cases=abort_cases,limb_specialization_cases=limb_cases,random_GPR_initialization=True,combinatorial_cases=combinatorial,point_couplings=couplings,original_id_emissions=idchecks,backward_variants=sorted(variants),rows=rows,all_checks=True,scope='Finite range, arithmetic, raw-label, routing and recovery evidence; complete analytical argument supplied separately. Not a production run or formal certification.')))
