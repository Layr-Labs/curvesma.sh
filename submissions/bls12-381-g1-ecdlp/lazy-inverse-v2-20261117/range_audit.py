"""Finite emitter boundary checks and exact production-parameter inequalities for lazy inverse."""
import random
from fractions import Fraction
import a64,rows,table_trade_audit
P,R,B=rows.P,rows.RA,1<<64
D=1<<376

def run():
 h=Fraction((P+D)**2,R)
 assert 9*P+h<R and 7*P<R and 14*P<2*R
 A=['x%d'%i for i in range(19,25)];O=['x17']+['x%d'%i for i in range(11,17)]
 W=['x0','x1','x3','x4','x25','x26']
 body=rows.rowprod(A,'x29',rows.yoff,O,W,'x27','x28')
 pg,lb=a64.a64_load('entry:\n'+''.join(' '+l+'\n' for l in body)+' ret\n')
 rng=random.Random(733);cases=0;maxr=0
 vals=[0,1,P-1,P,P+D-1,2*P,7*P,8*P,R-1]
 tests=[(a,y) for a in vals for y in [0,1,P-1,P]]
 tests += [(rng.randrange(R),rng.randrange(P+1)) for _ in range(80)]
 for a,y in tests:
  yr=rows.rows_of(y);assert all(0<=v<=P for v in yr)
  mem={}
  for j,v in enumerate(yr):
   for k in range(6):mem[0x10000+rows.yoff(j)+8*k]=(v>>(64*k))&(B-1)
  regs=table_trade_audit.words(A,a)|table_trade_audit.words(rows.PM,P)|{'x18':rows.N0,'x29':0x10000}
  out={};count=a64.a64_run(pg,lb,'entry',regs,mem,limit=1000,out_regs=out)
  rr=sum(out[O[j+1]]<<(64*j) for j in range(6));want=rows.rp_raw(a,yr)
  assert count==len(body)+1 and rr==want and rr<7*P and (rr-a*y*rows.RINV)%P==0
  assert all(out[r]==regs[r] for r in A)
  maxr=max(maxr,rr);cases+=1
 # Raw inverses substantially beyond the old range, including near-bound products.
 raw=[rows.tmr_raw(a*b) for a in [0,1,P-1,P+D-1] for b in [0,1,P-1,P+D-1]]
 raw += [rows.tmr_raw(rng.randrange(P+D)*rng.randrange(P+D)) for _ in range(200)]
 assert all(P<=v<8*P+h<R for v in raw) and sum(v>2*P for v in raw)>100
 # Prefix sums cannot overflow; signed e's represented addition wraps exactly once iff e<0.
 for _ in range(200):
  mu=rng.randrange(P+D);e=rng.randrange(-P+1,P);prod=mu*mu;acc=prod>>384;low=prod&(R-1)
  for i,(sh,w) in enumerate(zip(rows.TMR_SHIFTS,rows.TMR_WIDTHS)):
   acc+=rows.tmr_T(i,(low>>sh)&((1<<w)-1));assert 0<acc<R
  total=acc+e;rep=acc+(e%R)
  assert 0<total<R and rep//R==int(e<0) and rep%R==total
 return {'status':'passed','rp6_emitter_cases':cases,'rp6_instructions':len(body),'raw_inverse_cases':len(raw),'raw_inverse_gt_2p':sum(v>2*P for v in raw),'raw_inverse_max_over_p':max(raw)/P,'signed_prefix_cases':200,'exact_bounds':{'9p_plus_h_below_R':True,'7p_below_R':True,'14p_below_2R':True},'scope':'production-parameter algebraic inequalities plus finite full-width emitter and raw-inverse checks'}
