import random,bisect,a64
MASK=(1<<64)-1
def check(root):
 p,l=a64.a64_load((root/'search.s').read_text());rng=random.Random(81202);cases=maxcount=0
 for n in [0,1,2,3,7,17,128,257]:
  vals=sorted([rng.getrandbits(256) for _ in range(n)])
  if n>3:vals[1]=vals[0]
  base=0x100000+8*rng.randrange(50);threshold=0x200000;mem={}
  for j,v in enumerate(vals):
   for k in range(9):mem[base+72*j+8*k]=(v>>(64*k))&MASK if k<4 else rng.getrandbits(64)
  for q in [0,1,(1<<256)-1]+vals[:3]+[rng.getrandbits(256) for _ in range(17)]:
   for k in range(4):mem[threshold+8*k]=(q>>(64*k))&MASK
   saved=dict(mem);r={f'x{i}':rng.getrandbits(64) for i in range(31)};r.update(x0=base,x1=n,x2=threshold,sp=0x300000,x30=12345);out={}
   ct=a64.a64_run(p,l,'range_lower_bound',r,mem,limit=10000,out_regs=out)
   assert out['x0']==bisect.bisect_left(vals,q)
   assert ct<=16*(n.bit_length())+10,(n,ct)
   assert mem==saved and out['sp']==r['sp'] and out['x30']==r['x30']
   assert all(out['x'+str(i)]==r['x'+str(i)] for i in range(18,30))
   cases+=1;maxcount=max(maxcount,ct)
 return dict(cases=cases,max_observed_instructions=maxcount,per_iteration_bound=16,fixed_bound=10,production_n_2p28_bound=16*29+10,read_only=True,callee_saved_preserved=True,scope='native binary-search helper only; not full driver or score')
