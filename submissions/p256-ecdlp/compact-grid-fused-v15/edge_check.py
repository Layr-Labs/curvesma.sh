import random
import a64

def check(root):
 P=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
 R=1<<256;M=(1<<64)-1;D=R-P;L=R-(1<<224);NP=-pow(P,-1,R)%R
 pr,lab=a64.a64_load((root/'fused.s').read_text());rng=random.Random(20261101)
 vals=[0,1,2,P-2,P-1,P,P+1,R-2,R-1,R//2]
 pairs=[(a,c) for a in vals for c in vals]+[(rng.randrange(R),rng.randrange(R)) for _ in range(1000)]
 bad=kept=0;combos=set()
 for a,c in pairs:
  mem={}
  for ptr,v in [(4096,a),(8192,c)]:
   for i in range(4):mem[ptr+8*i]=(v>>(64*i))&M
  a64.a64_run(pr,lab,'k_fused',dict(x0=12288,x1=4096,x2=8192,x30=12345),mem,limit=1000)
  t=(a*a+(a*a*NP%R)*P)//R;combos.add((t//R,int(t%R<c)));z=t-c;h=z//R;u=(z-h*P)%R;k=z%P
  got=sum(mem[12288+8*i]<<(64*i) for i in range(4));assert got==u
  if u<L:assert u==k and u<P;kept+=1
  else:assert k<D or k>=P-D;bad+=1
 assert bad>0 and combos=={(0,0),(0,1),(1,0),(1,1)}
 return dict(cases=len(pairs),carry_borrow_combinations=len(combos),retained=kept,filtered=bad,retained_exact=True,all_rejected_true_keys_in_proved_set=True)
