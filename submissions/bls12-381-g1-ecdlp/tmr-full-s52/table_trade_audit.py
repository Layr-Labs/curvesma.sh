"""Exact table-trade construction/address/MADD audits; finite evidence only."""
import random
import a64,rows
P=rows.P
MASK=(1<<64)-1
A=['x%d'%i for i in range(6)]
D=['x%d'%i for i in range(6,12)]
Q=['x%d'%i for i in range(12,18)]
V=['x%d'%i for i in range(18,24)]
def chain(op,tail,d,a,b):
 return ['%s %s,%s,%s'%((op if j==0 else (tail if j==5 else op+'c')),d[j],a[j],b[j]) for j in range(6)]
def build_body(biased=False):
 add=lambda dst,a,b:['%s %s,%s,%s'%('adds' if j==0 else ('adc' if j==5 else 'adcs'),dst[j],a[j],b[j]) for j in range(6)]
 sub=lambda dst,a,b:['%s %s,%s,%s'%('subs' if j==0 else 'sbcs',dst[j],a[j],b[j]) for j in range(6)]
 out=add(V,A,Q) if biased else []
 src=V if biased else A
 out+=['stp %s,%s,[x24,#%d]'%(src[2*j],src[2*j+1],16*j) for j in range(3)]+['add x24,x24,#384']
 out+=add(A,A,D)+sub(V,A,Q)+['csel %s,%s,%s,lo'%(A[j],A[j],V[j]) for j in range(6)]+['subs x25,x25,#1','b.ne loop']
 assert len(out)==(30 if biased else 24)
 return out
def words(rs,v):return {r:(v>>(64*j))&MASK for j,r in enumerate(rs)}
def read(mem,addr):return sum(mem[addr+8*j]<<(64*j) for j in range(6))
def run():
 rng=random.Random(419);cases=0
 for i in range(7):
  body=build_body(i==0);pg,lb=a64.a64_load('entry:\nloop:\n'+''.join(' '+x+'\n' for x in body)+' ret\n')
  step=pow(2,rows.TMR_SHIFTS[i],P)*rows.RINV%P
  for start,count in [(0,19),((1<<rows.TMR_WIDTHS[i])-1,1),(rng.randrange(1<<rows.TMR_WIDTHS[i]),7)]:
   initial=step*start%P; addr=256+384*start+48*i-256
   regs=words(A,initial)|words(D,step)|words(Q,P)|{'x24':addr,'x25':count}
   mem={};n=a64.a64_run(pg,lb,'entry',regs,mem,limit=10000)
   assert n==len(body)*count+1
   assert all(read(mem,addr+384*j)==rows.tmr_T(i,start+j) for j in range(count))
   assert set(mem)=={addr+384*j+8*k for j in range(count) for k in range(6)}
   cases+=1
 # Packed bounds, last-row boundaries and disjoint correction slot.
 alloc=512;B=alloc+256;end=alloc+rows.TMR_BYTES
 for i,w in enumerate(rows.TMR_WIDTHS):
  for c in [0,1,(1<<w)-1]:
   lo=B+384*c+48*i-256;hi=lo+48
   assert alloc<=lo<hi<=end and lo%16==0
 for c in [0,1,511]:
  lo=B+384*c+80;assert lo==alloc+384*c+7*48 and lo+48<=end
 assert B+384*((1<<55)-1)+80+48==end
 # Madd interpreter: modulo, input/output alias, and zero register.
 pg,lb=a64.a64_load('entry:\n madd x0,x0,x1,x2\n madd x1,x0,x1,xzr\n ret\n')
 for _ in range(50):
  aa,bb,cc=[rng.getrandbits(64) for _ in range(3)];out={}
  assert a64.a64_run(pg,lb,'entry',{'x0':aa,'x1':bb,'x2':cc},{},out_regs=out)==3
  x=(aa*bb+cc)&MASK;assert out['x0']==x and out['x1']==(x*bb)&MASK
 # Wrong byte stride and omitted +p give distinct first/last records: negative controls.
 assert rows.tmr_T(0,0)==P and rows.tmr_T(0,0)!=0
 assert B+384*((1<<55)-1)+48*6-256 != B+336*((1<<55)-1)+48*6-256
 return {'status':'passed','construction_cases':cases,'tmr_first_entry_instructions':30,'tmr_other_entry_instructions':24,'madd_alias_cases':50,'address_boundaries':'passed','negative_controls':2,'tmr_bytes':rows.TMR_BYTES,'tmr_entries_charged':rows.TMR_ENTRIES,'scope':'finite exact construction emulator plus algebraic layout assertions; not production table construction'}
