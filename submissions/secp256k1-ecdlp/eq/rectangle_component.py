"""Native ygen, strict-memory emulator, formula and source-codec controls.
Research evidence only; does not execute the whole attack or price CPU time.
"""
import ctypes
import hashlib
import json
import subprocess
from pathlib import Path
OUT=Path('/out')
P=json.loads(Path('/profile.json').read_text())
p=int(P['field_modulus'],16)
r=int(P['subgroup_order'],16)
G=tuple(int(P['generator_P'][k],16) for k in ('x','y'))
assert p==2**256-2**32-977
beta=int('7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee',16)
mask64=2**64-1
def add(A,B):
 if A is None:return B
 if B is None:return A
 x,y=A;u,v=B
 if x==u:
  if (y+v)%p==0:return None
  z=3*x*x*pow(2*y,-1,p)%p
 else:z=(v-y)*pow(u-x,-1,p)%p
 w=(z*z-x-u)%p
 return w,(z*(x-w)-y)%p
def mul(k,A=G):
 B=None
 while k:
  if k&1:B=add(B,A)
  A=add(A,A);k>>=1
 return B
def record(Y):
 x,y=Y
 vals=[0,0]
 def words(v,n=4):return [(v>>(64*j))&mask64 for j in range(n)]
 vals+=words(y*y%p)+words(pow(y,4,p))+words(pow(x,4,p))
 vals+=sum((words(y*pow(2,64*k,p)%p) for k in range(4)),[])
 vals+=words(-((x<<128)//p)%2**128,2)+words(-(((beta*x%p)<<128)//p)%2**128,2)
 for n,v in enumerate((x,pow(x,3,p),3*x*x%p)):
  vals += [((v>>(32*k))&(2**32-1))*256+64*(k>>1)+(n<<40) for k in range(8)]
 assert len(vals)==58
 return vals
def formula(m,n):return -23285+5*m+22598*n+1606*m*n+45680*(m-1).bit_length()
src=Path('/source/kernels/secp_fe.s').read_text()
wrapper='''
 .text
 .balign 4
 .global yg_test
 .type yg_test, %function
yg_test:
 sub sp,sp,#96
 stp x19,x20,[sp]
 stp x21,x22,[sp,#16]
 stp x23,x24,[sp,#32]
 stp x25,x26,[sp,#48]
 stp x27,x28,[sp,#64]
 stp x29,x30,[sp,#80]
 cmp x7,#1
 bl ygen
 ldp x19,x20,[sp]
 ldp x21,x22,[sp,#16]
 ldp x23,x24,[sp,#32]
 ldp x25,x26,[sp,#48]
 ldp x27,x28,[sp,#64]
 ldp x29,x30,[sp,#80]
 add sp,sp,#96
 ret
 .section .note.GNU-stack,"",%progbits
'''
(OUT/'native-ygen.s').write_text(src+wrapper)
cp=subprocess.run(['gcc','-shared','-fPIC','-march=armv8-a','-o',str(OUT/'native-ygen.so'),str(OUT/'native-ygen.s')],capture_output=True,text=True,timeout=12)
(OUT/'native-compile.log').write_text(cp.stdout+cp.stderr)
assert cp.returncode==0
lib=ctypes.CDLL(str(OUT/'native-ygen.so'))
f=lib.yg_test
f.argtypes=[ctypes.c_void_p,ctypes.c_uint64,ctypes.c_uint64,ctypes.c_void_p,ctypes.c_void_p,ctypes.c_void_p,ctypes.c_void_p,ctypes.c_uint64]
f.restype=ctypes.c_uint64
C,T,S=mul(7),mul(11),mul(13)
def point(A):return (ctypes.c_uint64*8)(*[((v>>(64*j))&mask64) for v in A for j in range(4)])
def native(m,n,c=C,carry=0,verify=True):
 outsize=512*m*n;scratchsize=96*m+1056
 output=ctypes.create_string_buffer(b'\xa5'*(outsize+128))
 scratch=ctypes.create_string_buffer(b'\xa5'*(scratchsize+128))
 cbuf,tbuf,sbuf=point(c),point(T),point(S)
 op,sp=ctypes.addressof(output)+64,ctypes.addressof(scratch)+64
 flag=f(op,m,n,cbuf,tbuf,sbuf,sp,carry)
 assert output.raw[:64]==output.raw[64+outsize:128+outsize]==b'\xa5'*64
 assert scratch.raw[:64]==scratch.raw[64+scratchsize:128+scratchsize]==b'\xa5'*64
 checked=0
 if verify:
  row=c
  for nn in range(n):
   Y=row
   for mm in range(m):
    v=(ctypes.c_uint64*58).from_address(op+512*(nn*m+mm))
    assert list(v)==record(Y),(m,n,mm,nn)
    checked+=58
    Y=add(Y,T)
   row=add(row,S)
  assert flag==0
 return dict(m=m,n=n,carry=carry,flag=flag,record_words_checked=checked,scratch_bytes=scratchsize,guards_ok=True,
  payload_sha256=hashlib.sha256(ctypes.string_at(op,outsize)).hexdigest())
native_results=[native(8192,4,carry=0),native(8192,4,carry=1),native(64,512),native(4,2)]
assert native_results[0]['payload_sha256']==native_results[1]['payload_sha256']
bad=native(1,2,c=S,verify=False)
assert bad['flag']!=0
codec=0
for u in range(32768):
 m,n=u%8192,u//8192
 assert 0<=m<8192 and 0<=n<4 and u==8192*n+m
 for j in [0,1,511,512,1023]:
  l=1024*u+j
  assert l<2**25
  for sign in [1,-1]:
   code=2 if sign==1 else 4
   lp=l if sign==1 else (l-1)%2**25
   for path in range(3):
    high=code|(lp<<3)
    transported=high^((1<<28)-1) if path==1 else high
    restored=transported^((1<<28)-1) if path==1 else transported
    cc=restored&7;ll=restored>>3
    assert cc==code
    recovered=ll if cc==2 else (ll+1)%2**25
    assert recovered==l and divmod(recovered,1024)==(u,j)
    codec+=1
result=dict(scope='ygen and codec component only, not a full attack or protected submission gate',
 authorship='GPT 6.1 Sol / MAX / Codex via T3 Code; inherited kernel/a64/base record formulas credited',
 source_sha256=hashlib.sha256(src.encode()).hexdigest(),native_compile_ok=True,native=native_results,
 native_negative=bad,production_codec_cases=codec,
 production_ygen_count=formula(8192,4),source_record_count=32768,
 scratch_formula='1024 + 64*M + 32*(M+1) = 96*M+1056',
 instruction_count_scope='Formula/count evidence is the declared experiment and analytical loop derivation; native execution is correctness evidence, never timing-derived work.',
 completed=True)
(OUT/'component-result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['native','emulator']},indent=2))
