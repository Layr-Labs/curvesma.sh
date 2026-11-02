"""Finite integration: native merge helpers construct all sorted lists."""
import a64
MASK=(1<<64)-1
class Builder:
 def __init__(self,root):
  self.prog,self.labels=a64.a64_load((root/'merge.s').read_text());self.calls=0;self.instructions=0;self.records=0
 def merge(self,a,b,off):
  mem={};left=0x10000;right=0x30000;dst=0x50000
  def words(z):return [(z[0]>>(64*i))&MASK for i in range(4)]+[(z[1]>>(64*i))&MASK for i in range(4)]+[z[2]]
  for base,seq in [(left,a),(right,b)]:
   for j,z in enumerate(seq):
    for k,w in enumerate(words(z)):mem[base+72*j+8*k]=w
  ct=a64.a64_run(self.prog,self.labels,'range_merge',dict(x0=left,x1=len(a),x2=right,x3=len(b),x4=dst,x5=off,x30=12345,sp=0x70000),mem,limit=100000)
  assert ct<=64*(len(a)+len(b))+8
  ans=[]
  for j in range(len(a)+len(b)):
   ws=[mem[dst+72*j+8*k] for k in range(9)];ans.append((sum(ws[k]<<(64*k) for k in range(4)),sum(ws[k+4]<<(64*k) for k in range(4)),ws[8]))
  self.calls+=1;self.instructions+=ct;self.records+=len(ans);return ans
 def build(self,items):
  assert items
  runs=[[x] for x in items]
  while len(runs)>1:runs=[self.merge(runs[i],runs[i+1] if i+1<len(runs) else [],32) for i in range(0,len(runs),2)]
  seq=runs[0];assert seq==sorted(items,key=lambda z:(z[1],z[2]));nodes=[]
  def rec(lo,hi):
   nd={'lo':lo,'hi':hi,'base':0x100000000+len(nodes)*0x10000};nodes.append(nd)
   if hi-lo==1:nd['items']=seq[lo:hi]
   else:
    mid=(lo+hi)//2;a=rec(lo,mid);b=rec(mid,hi);nd['children']=(a,b);nd['items']=self.merge(a['items'],b['items'],0)
   assert nd['items']==sorted(seq[lo:hi],key=lambda z:(z[0],z[2]));return nd
  return rec(0,len(seq)),seq,nodes
