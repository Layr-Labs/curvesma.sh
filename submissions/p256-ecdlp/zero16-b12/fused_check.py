import random
from pathlib import Path
import a64,range_witness
def check(root):
 root=Path(root)
 prog,labels=a64.a64_load(range_witness.fused_base(root))
 rg=random.Random(97531)
 W=0x1000000;RB=0x2000000;SN=0x3000000;PO=0x90000
 n=0
 for t in range(32):
  dr=t>=24;d=rg.randrange(1<<35);w=0x8000000+(rg.randrange(1<<20)<<12);m={};m[RB+d*8]=SN+65536 if dr else SN;sp=0x70000
  m[sp]=W;m[sp+8]=5|(7<<35);m[sp+16]=RB;m[sp+24]=SN+65536;m[sp+32]=PO;pm=[];o={}
  c=a64.a64_run(prog,labels,'flush_publish_b',dict(x14=w+4096,x13=W+d*8,x7=d,sp=sp,x30=12345),m,limit=5000,out_regs=o,dev=(PO,32,lambda k,v:pm.append((k,v))))
  assert o['x30']==12345 and o['sp']==sp and m[W+d*8]==w+32,(t,c)
  assert m[RB+d*8]==(SN+65536+4096 if dr else SN+4096) and c==(12 if dr else 15),(t,c)
  g={k:v for k,v in pm}
  if dr:assert pm==[],(t,pm)
  else:assert g.get(0)==d and g.get(8)==SN and g.get(16)==w+32 and g.get(24)==128,(t,g)
  n+=1
 return dict(flush_cases=n,flush_taken=15,flush_dropped=12,frame_layout_ok=True,structural_liveness=structural_liveness(root),scope='finite block-flush-subroutine unit checks; integration via probe fused runs')
def _regs(s):
 import re
 o=[int(m.group(1)) for m in re.finditer(r'\bx(\d+)\b',s)]+[int(m.group(1)) for m in re.finditer(r'\bw(\d+)\b',s)]
 if re.search(r'\bsp\b',s):o.append(31)
 return o
def _rw(s):
 import re
 mn=s.split()[0];rs=_regs(s)
 if not rs:return set(),set()
 if mn in ('stp','str','cbz','cmp','b','b.ne','b.cs','b.cc','b.hs','b.lo','tst','ret'):
  w=set()
  if mn in ('stp','str') and (s.rstrip().endswith('!') or re.search(r'\],\s*#',s)):w.add(rs[-1])
  return set(rs)-w,w
 if mn in ('ldp','ldr'):
  k=2 if mn=='ldp' else 1;r=set(rs[k:]);w=set(rs[:k])
  if re.search(r'\],\s*#',s):w.add(rs[k])
  return r,w
 return set(rs[1:]),{rs[0]}
def structural_liveness(root):
 import re
 from pathlib import Path
 root=Path(root)
 asm=range_witness.base_generated_assembly(root).splitlines()
 pk={}
 for yl in [0,1]:
  for cy in [0,1]:
   suf=f'_y{yl}_c{cy}'
   i0=next(i for i,l in enumerate(asm) if l.startswith(f'pair_bwd{suf}:'))
   nx=[i for i,l in enumerate(asm) if i>i0 and re.match(r'pair_bwd_y\d_c\d:',l)]
   i1=min(nx) if nx else len(asm)
   bd=[(i0+j,t) for j,l in enumerate(asm[i0:i1]) for t in [l.strip()] if t and not t.endswith(':') and not t.startswith('.')]
   lo=set();peak=0
   for _ in range(8):
    lv=set(lo);peak=0
    for (ln,t) in reversed(bd):
     r,w=_rw(t);lv=(lv-w)|r;peak=max(peak,len(lv))
    lo=set(lv);lo.add(30)
   assert peak<=30,(suf,peak)
   pk[suf]=peak
   ou=[ln for (ln,t) in bd if re.match(r'stp x\d+,x\d+,\[x21',t)]
   assert len(ou)==4,(suf,ou)
   a1,b1,a2,b2=ou
   for ln in (a1,a2):
    lv=set(lo)
    for (m,t) in reversed(bd):
     if m<ln:break
     r,w=_rw(t);lv=(lv-w)|r
    assert set(range(9,17))<=set(range(31))-lv,(suf,ln)
   af=[t for (ln,t) in bd if b1<ln<a2]
   f7=next(t for t in af if re.search(r'\bx7\b',t));f8=next(t for t in af if re.search(r'\bx8\b',t))
   assert re.match(r'(mul|umulh|mov|ldr|ldp) x[78]\b',f7) and re.match(r'(mul|umulh|mov|ldr|ldp) x[78]\b',f8),(suf,f7,f8)
   for (ln,t) in bd:
    if ln>a1:assert not re.match(r'ldp .*\[x21',t),(suf,ln,t)
   wr=[(ln,t) for (ln,t) in bd if re.match(r'(mov|sub|add|ldp|ldr) x20\b',t) or re.search(r'\[x20\],\s*#',t)]
   assert len(wr)==2,(suf,wr)
   co=None
   for (ln,t) in bd:
    if ln>b1 and t.split()[0] in ('adc','adcs','sbc','sbcs','b.cs','b.cc','b.hs','b.lo'):
     co=(ln,t);break
   assert co and any(t.split()[0] in ('adds','adcs','subs','sbcs','cmp','tst') for (ln,t) in bd if b1<ln<co[0]),(suf,co)
 return dict(peak_live_regs=pk,output_scratch_free=True,key_halves_dead_after_output1=True,no_batch_reads_after_outputs=True,x20_pure_walk=True,flag_boundary_ok=True)
