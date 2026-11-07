"""Focused exact arithmetic and certificate checks; no fwd/bwd or routing execution."""
import contextlib
import io
import json
import os
import random
import sys
import time

started=time.monotonic()
cfg=json.load(sys.stdin)
here=os.path.dirname(os.path.abspath(__file__))
here=next(b for b in (here,os.path.join(here,'files'),'/snapshot/files')
 if os.path.isfile(os.path.join(b,'kernels','a64.py'))and os.path.isfile(os.path.join(b,'certificate.py')))
sys.path.insert(0,here)
sys.path.insert(0,os.path.join(here,'kernels'))
import a64
capture=io.StringIO()
with contextlib.redirect_stdout(capture):
 import certificate as cert
lines=capture.getvalue().splitlines()
p,C0,M=cert.p,cert.C0,(1<<64)-1
rng=random.Random(cfg['seed'])
vals=[0,1,p-1,p,2**256-1]+[2**k+d for k in range(37,256,37)for d in (-1,0,1)]
vals += [rng.getrandbits(256)for _ in range(128)]
radix_ok=True
cases=0
for v in vals:
 digits=[(v>>(37*k))&((1<<37)-1)for k in range(7)]
 radix_ok &= sum(d<<(37*k)for k,d in enumerate(digits))==v and digits[-1]<2**34
 for c,add in ((p-1,p-1),(rng.randrange(p),rng.randrange(p))):
  base=rng.randrange(p)
  z=base+sum((c*d*pow(2,37*k,p)+(add if k==0 else 0))%p for k,d in enumerate(digits))
  t=(z%(1<<256))+(z>>256)*C0
  folded=t if t<1<<256 else t-(1<<256)+C0
  radix_ok &= z<8*p and 0<=folded<1<<256 and folded%p==(base+c*v+add)%p
  cases+=1
programs={}
counts={}
for name in ('leaf','route','secp_fe'):
 with open(os.path.join(here,'kernels',name+'.s'))as f:
  programs[name]=a64.a64_load(f.read())
 counts[name]=len(programs[name][0])
expanded_ok=counts==dict(leaf=25196,route=16084,secp_fe=1524)and sum(counts.values())==cert.KERNEL_INSTR

def put(mem,a,v):
 for k in range(4):
  mem[a+8*k]=(v>>(64*k))&M

def get(mem,a):
 return sum(mem[a+8*k]<<(64*k)for k in range(4))

constructor_ok=True
runs=[]
for c,add,n in ((rng.randrange(p),rng.randrange(p),17),(p-1,p-1,9),(1,p-2,5),(0,rng.randrange(p),3)):
 mem={}
 for k in range(8):
  put(mem,0x1000+32*k,c*pow(2,37*k,p)%p)
 put(mem,0x1200,add)
 prog,labels=programs['route']
 steps=a64.a64_run(prog,labels,'tbld8',dict(x0=0x100000,x1=0x1000,x2=0x1200,x3=n,sp=0x7000000,x30=0x5E5E5E5E5E5E5E5E),mem)
 words={a for a in mem if a>=0x100000}
 ok=steps==92+96*n and words=={0x100000+8*w for w in range(32*n)}and all(
  get(mem,0x100000+256*t+32*k)==(c*t*pow(2,37*k,p)+(add if k==0 else 0))%p
  for t in range(n)for k in range(8))
 constructor_ok &= ok
 runs.append(dict(rows=n,instructions=steps,passed=ok))
certificate_ok=not cert.FAILED and lines[-1]=='PASS'
print(json.dumps(dict(checks_certificate=dict(passed=certificate_ok,failed=cert.FAILED,checks=sum(' ok' in line for line in lines)),
 radix37=dict(cases=cases,passed=radix_ok),expanded_instructions=counts,expanded_counts_passed=expanded_ok,
 tbld8=dict(runs=runs,passed=constructor_ok),scope='Exact certificate, radix identities, assembly expansion and short table construction only; fwd/bwd and routing are not executed.',
 passed=certificate_ok and radix_ok and expanded_ok and constructor_ok,seconds=round(time.monotonic()-started,3)),indent=1))
