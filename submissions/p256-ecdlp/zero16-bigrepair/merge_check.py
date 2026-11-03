import a64,random
MASK=(1<<64)-1
def check(root):
 p,l=a64.a64_load((root/'merge.s').read_text());rng=random.Random(91202);cases=0;maxper=0
 for off in [0,32]:
  for na in [0,1,2,7,16]:
   for nb in [0,1,3,9,17]:
    for mode in range(3):
     def rec(j):
      return (rng.getrandbits(256) if mode==0 else j%3,rng.getrandbits(256) if mode==0 else j%2,j)
     key=lambda z:(z[off//32],z[2]);a=sorted([rec(j) for j in range(na)],key=key);b=sorted([rec(j+na) for j in range(nb)],key=key)
     mem={};left=0x10000;right=0x30000;dst=0x50000
     def words(z):return [(z[0]>>(64*i))&MASK for i in range(4)]+[(z[1]>>(64*i))&MASK for i in range(4)]+[z[2]]
     for base,seq in [(left,a),(right,b)]:
      for j,z in enumerate(seq):
       for k,w in enumerate(words(z)):mem[base+72*j+8*k]=w
     saved=dict(mem);r={f'x{i}':rng.getrandbits(64) for i in range(31)};r.update(x0=left,x1=na,x2=right,x3=nb,x4=dst,x5=off,sp=0x70000,x30=12345);out={}
     ct=a64.a64_run(p,l,'range_merge',r,mem,limit=30000,out_regs=out)
     expected=sorted(a+b,key=key)
     for j,z in enumerate(expected):assert [mem[dst+72*j+8*k] for k in range(9)]==words(z)
     assert all(mem[k]==v for k,v in saved.items())
     assert out['x0']==dst+72*(na+nb) and out['sp']==r['sp'] and out['x30']==r['x30']
     assert all(out['x'+str(i)]==r['x'+str(i)] for i in range(18,30))
     assert ct<=64*(na+nb)+8,(ct,na,nb)
     cases+=1;maxper=max(maxper,ct)
 return dict(cases=cases,record_bound=64,fixed_bound=8,source_preserved=True,callee_saved_ok=True,both_coord_keys=True,max_observed=maxper,scope='native merge helper only; caller/tree construction still separate')
