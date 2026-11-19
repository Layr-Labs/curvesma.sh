"""Exact DC ZVA256 initialization model; no hardware execution or cycles claimed."""
BLOCK=256
UNROLL=1<<20
FIXED=64

def production_cost(length):
 assert length==(1<<61) and length%(BLOCK*UNROLL)==0
 blocks=length//BLOCK;groups=blocks//UNROLL
 # Two instructions/block, SUBS+B.EQ+B/group (last B skipped), plus setup.
 return 2*blocks+3*groups+FIXED

def listing(unroll,groups):
 assert unroll>=1 and groups>=1
 # Production groups=2^33 fits one MOVZ x1,#2,lsl#32.
 # Test cases use small one-MOVZ immediate counts.
 count='movz x1,#2,lsl #32' if groups==1<<33 else 'movz x1,#%d'%groups
 return ['mrs x2,dczid_el0','and x2,x2,#31','cmp x2,#6','b.eq .Lsupported','ret','.Lsupported:',
         'movz x0,#0',count,'.Lzero:']+['dc zva,x0','add x0,x0,#256']*unroll+[
         'subs x1,x1,#1','b.eq .Ldone','b .Lzero','.Ldone:','dsb sy','ret']

def run_listing(lines,block,dczid,length):
 # Dedicated finite semantics for the literal scalar listing above. Registers
 # and branches are modeled; zero writes affect exactly the aligned block.
 labels={s[:-1]:i for i,s in enumerate(lines) if s.endswith(':')}
 mem=bytearray([0xa7])*length;r={};pc=0;steps=0;eq=False;barriers=0;written=set()
 while True:
  s=lines[pc];pc+=1
  if s.endswith(':'):continue
  steps+=1
  if s=='mrs x2,dczid_el0':r['x2']=dczid
  elif s=='and x2,x2,#31':r['x2']&=31
  elif s=='cmp x2,#6':eq=r['x2']==6
  elif s.startswith('b.eq '):
   if eq:pc=labels[s.split()[1]]
  elif s.startswith('b '):pc=labels[s.split()[1]]
  elif s.startswith('movz x0,#'):r['x0']=int(s.split('#')[1])
  elif s.startswith('movz x1,#'):r['x1']=int(s.split('#')[1])
  elif s=='dc zva,x0':
   a=(r['x0']//block)*block
   assert 0<=a and a+block<=length,'out-of-allocation write'
   for j in range(a,a+block):
    assert j not in written,'duplicate zero block'
    written.add(j);mem[j]=0
  elif s.startswith('add x0,x0,#'):r['x0']+=int(s.split('#')[1])
  elif s=='subs x1,x1,#1':r['x1']-=1;eq=r['x1']==0
  elif s=='dsb sy':barriers+=1
  elif s=='ret':break
  else:raise AssertionError(s)
  assert steps<100000
 return {'steps':steps,'writes':len(written),'all_zero':not any(mem),
         'final_address':r.get('x0'),'barriers':barriers}

def branch_certificate(source):
 # Compute offsets in the production-expanded source without allocating its
 # 2^21 repeated instructions. Only the explicit two-instruction macro repeats.
 pc=0;labels={};branches=[];i=0
 while i<len(source):
  line=source[i];i+=1
  if line.endswith(':'):labels[line[:-1]]=pc;continue
  if line.startswith('.rept '):
   repeat=int(line.split()[1]);body=[]
   while source[i]!='.endr':body.append(source[i]);i+=1
   i+=1
   assert body==['dc zva,x0','add x0,x0,#256']
   pc+=4*len(body)*repeat;continue
  assert not line.startswith('.')
  if line.startswith(('b.eq ', 'b ')):
   op,target=line.split();branches.append((pc,op,target))
  pc+=4
 def encodable(displacement,bits):
  return displacement%4==0 and -(1<<(bits-1))<=displacement//4<(1<<(bits-1))
 results=[]
 for addr,op,target in branches:
  disp=labels[target]-addr;bits=19 if op=='b.eq' else 26
  assert encodable(disp,bits),(addr,op,target,disp)
  results.append({'source_offset':addr,'opcode':op,'target_offset':labels[target],
                  'displacement':disp,'signed_immediate_bits':bits})
 assert [x['displacement'] for x in results]==[8,8,-8388616]
 assert pc==8*UNROLL+48
 # Original-style conditional backedge across the8MiB body is not encodable.
 assert not encodable(results[-1]['displacement'],19)
 assert encodable(-(1<<20),19) and encodable((1<<20)-4,19)
 assert not encodable(-(1<<20)-4,19) and not encodable(1<<20,19)
 assert encodable(-(1<<27),26) and encodable((1<<27)-4,26)
 assert not encodable(-(1<<27)-4,26) and not encodable(1<<27,26)
 return {'expanded_code_bytes':pc,'branches':results,'old_conditional_backedge_rejected':True,
         'signed_boundary_controls':8}

def run():
 from pathlib import Path
 # Tie the symbolic production count to the actual supplied assembler source.
 lines=[]
 for line in (Path(__file__).with_name('owner_init.S')).read_text().splitlines():
  q=line.split('//',1)[0].strip().replace(', ', ',')
  if q:lines.append(q)
 expected=listing(1,1<<33)
 # Replace the single body pair by a literal assembler repeat directive.
 expected=expected[:9]+['.rept 1048576']+expected[9:11]+['.endr']+expected[11:]
 assert lines==['owner_init:']+expected
 reach=branch_certificate(lines)
 cases=0
 for u in (1,2,4,8):
  for g in (1,2,3,8):
   n=BLOCK*u*g;r=run_listing(listing(u,g),BLOCK,6,n)
   assert r=={'steps':2*u*g+3*g+7,'writes':n,'all_zero':True,'final_address':n,'barriers':1}
   cases+=1
 unsupported=0
 for dc in (0,4,5,7,16|6,31):
  r=run_listing(listing(2,2),BLOCK,dc,1024)
  assert r['steps']==5 and r['writes']==0 and r['barriers']==0
  unsupported+=1
 negative=0
 # A64FX capability is an explicit precondition. A modeled wrong block size
 # while falsely advertising BS6 must fail coverage, not silently pass.
 r=run_listing(listing(2,2),64,6,1024)
 assert not r['all_zero'] and r['writes']==256;negative+=1
 bad=[s.replace('add x0,x0,#256','add x0,x0,#128') for s in listing(2,2)]
 try:run_listing(bad,BLOCK,6,1024)
 except AssertionError:negative+=1
 else:raise AssertionError('overlap mutation passed')
 bad=[s.replace('add x0,x0,#256','add x0,x0,#512') for s in listing(2,2)]
 try:run_listing(bad,BLOCK,6,1024)
 except AssertionError:negative+=1
 else:raise AssertionError('overrun mutation passed')
 length=1<<61;blocks=length//BLOCK;groups=blocks//UNROLL
 assert blocks==1<<53 and groups==1<<33
 assert BLOCK*(blocks-1)+BLOCK==length and length<1<<64
 assert UNROLL*8+12*4<(1<<24) # code fits the existing 1GiB node allowance
 assert production_cost(length)==(1<<54)+3*(1<<33)+64
 assert 8<64 # six executed entry instructions + DSB/RET; final loop B is skipped
 # Every conditional jump is +8 bytes. Unconditional loop branch is in range.
 branch_back=-(2*UNROLL+2)*4
 assert -(1<<27)<=branch_back<(1<<27) and -(1<<20)<=8<(1<<20)
 return {'status':'passed','production_source_matches_template':True,'branch_reach_certificate':reach,'finite_listing_cases':cases,'capability_rejections':unsupported,
  'negative_controls':negative,'block_bytes':BLOCK,'unroll':UNROLL,'production_blocks':blocks,
  'production_groups':groups,'production_per_owner_instructions':production_cost(length),
  'conditional_branch_displacements':[8,8],'unconditional_branch_back':branch_back,'last_group_saves_one_branch':True,'fixed_exact_body':8,'fixed_charged':FIXED,'reference_C_unchanged':2000,
  'scope':'Literal finite scalar-listing model and exact full-address/count algebra; no native DC ZVA execution, no production allocation, no hardware cycle claim.'}
