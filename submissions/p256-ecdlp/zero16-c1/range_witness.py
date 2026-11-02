import re
import chain_witness
def base_generated_assembly(root):
 source,meta=chain_witness.generate(root)
 source=chain_witness.normalize_transform(source)
 t=source;f,b=t.split('.globl pair_bwd',1)
 f=f.replace('[x20], #64','[x20], #72').replace('[x20, #-48]','[x20, #-56]')
 assert '[x20], #64' not in f
 parts=[f]
 for yl in [0,1]:
  for carry in [0,1]:
   z='.globl pair_bwd'+b
   old=' b.hs abs_store_3\n subs x11,xzr,x11\n sbcs x12,xzr,x12\n sbcs x13,xzr,x13\n sbc x14,xzr,x14\n';assert old in z;z=z.replace(old,'')
   if yl:
    # Four-limb operand swap, one limb per line (backward y-partition variant).
    z=z.replace('subs x11, x3, x7','subs x11, x7, x3')
    z=z.replace('sbcs x12, x4, x8','sbcs x12, x8, x4')
    z=z.replace('sbcs x13, x5, x9','sbcs x13, x9, x5')
    z=z.replace('sbcs x14, x6, x10','sbcs x14, x10, x6')
   old=' adc x15, xzr, xzr\n cbz x15,k_add_store_2\n';assert old in z;z=z.replace(old,'')
   if not carry:
    old=' adds x11,x11,#1\n sbcs x12,x12,x0\n sbcs x13,x13,xzr\n sbc x14,x14,x2\nk_add_store_2:';assert old in z;z=z.replace(old,'k_add_store_2:')
   z=z.replace('sub x20,x20,#64','sub x20,x20,#72')
   suffix=f'_y{yl}_c{carry}'
   for lab in re.findall(r'^([A-Za-z_]\w*):',z,re.M):z=re.sub(r'\b'+lab+r'\b',lab+suffix,z)
   parts.append(z)
 return '\n'.join(parts)
def unroll8(code,loop,done):
 start=code.index(loop+':\n');stop=code.index(done+':\n',start)
 block=code[start+len(loop)+2:stop]
 end=f' subs x19,x19,#1\n b.ne {loop}\n'
 assert block.endswith(end)
 body=block[:-len(end)]
 assert not re.search(r'\bx19\b',body)
 t=loop+'_tail';u=loop+'_single';g=loop+'_eight'
 new=f' subs x19,x19,#8\n b.lo {t}\n{g}:\n'
 for k in range(8):
  b=body
  for lab in re.findall(r'^([A-Za-z_]\w*):',body,re.M):b=re.sub(r'\b'+lab+r'\b',lab+f'_u{k}',b)
  new+=b
 new+=f' subs x19,x19,#8\n b.hs {g}\n{t}:\n adds x19,x19,#8\n cbz x19,{done}\n{u}:\n'+body+f' subs x19,x19,#1\n b.ne {u}\n'
 return code[:start]+new+code[stop:]
def pingpong_templates(body):
 km=body.index(' // inline k_mul')
 tail=body.index(' add x23,x21,#32')
 assert body.count(' // inline k_mul')==1 and body.count(' add x23,x21,#32')==1
 head,region,rest=body[:km],body[km:tail],body[tail:]
 assert region.count(' ldp x3, x4, [x23]\n')==1 and region.count(' ldp x5, x6, [x23, #16]\n')==1
 assert re.search(r'\bx1\b',region) is None
 region=region.replace(' ldp x3, x4, [x23]\n','').replace(' ldp x5, x6, [x23, #16]\n','')
 rest=rest.replace(' add x23,x21,#32\n','',1)
 assert 'add x23,x21,#32' not in rest
 outs=[]
 for bank_left,bank_out in [([28,29,30,18],[24,25,26,27]),([24,25,26,27],[28,29,30,18])]:
  mapping=dict(zip([3,4,5,6],bank_left))
  mapping.update(zip([9,10],bank_out[:2]));mapping.update(zip([24,25],bank_out[2:]));mapping.update({26:15,27:16})
  new_region=re.sub(r'\b([xw])(\d+)\b',lambda m:m.group(1)+str(mapping.get(int(m.group(2)),int(m.group(2)))),region)
  new_body=head+new_region+rest
  assert len(new_body.splitlines())==len(body.splitlines())-3
  assert 'x23' not in new_body
  outs.append(new_body)
 return outs[0],outs[1]
def unroll8_fwd(code,loop,done):
 start=code.index(loop+':\n');stop=code.index(done+':\n',start)
 block=code[start+len(loop)+2:stop]
 end=f' subs x19,x19,#1\n b.ne {loop}\n'
 assert block.endswith(end)
 body=block[:-len(end)]
 assert not re.search(r'\bx19\b',body)
 even,odd=pingpong_templates(body)
 t=loop+'_tail';u=loop+'_single';g=loop+'_eight';e=loop+'_epilogue'
 new=f' str x30,[sp,#-16]!\n subs x19,x19,#8\n b.lo {t}\n ldp x28,x29,[x21,#-32]\n ldp x30,x18,[x21,#-16]\n{g}:\n'
 for k in range(8):
  b=even if k%2==0 else odd
  for lab in re.findall(r'^([A-Za-z_]\w*):',b,re.M):b=re.sub(r'\b'+lab+r'\b',lab+f'_u{k}',b)
  new+=b
 new+=f' subs x19,x19,#8\n b.hs {g}\n{t}:\n adds x19,x19,#8\n cbz x19,{e}\n sub x23,x21,#32\n{u}:\n'+body+f' subs x19,x19,#1\n b.ne {u}\n{e}:\n ldr x30,[sp],#16\n'
 return code[:start]+new+code[stop:]
def generated_assembly(root):
 code=fused_base(root)
 for loop,done in [('pf_loop_lo','pf_success_lo'),('pf_loop_hi','pf_success_hi')]:
  code=unroll8_fwd(code,loop,done)
 for yl in [0,1]:
  for carry in [0,1]:
   code=unroll8(code,'pb_loop'+f'_y{yl}_c{carry}','pb_done'+f'_y{yl}_c{carry}')
 return code
def block_production(root):
 return unroll8((root/'block.s').read_text(),'pk_loop','pk_done')
FLUSH_SUB_B=""".globl flush_publish_b
flush_publish_b:
 ldr x9,[sp,#16]
 ldr x10,[sp,#24]
 ldr x11,[sp,#32]
 mov x12,#4064
 sub x14,x14,x12
 str x14,[x13]
 add x9,x9,x7,lsl #3
 mov x12,#4096
 ldadd x12,x13,[x9]
 cmp x13,x10
 b.hs fpb_skip
 mov x12,#128
 stp x7,x13,[x11]
 stp x14,x12,[x11,#16]
fpb_skip:
 ret
"""
def fused_base(root):
 code=base_generated_assembly(root)
 parts=re.split(r'(?=\.globl pair_bwd_y\d_c\d)',code)
 assert len(parts)==5
 out=[parts[0]]
 for z in parts[1:]:
  m=re.match(r'\.globl pair_bwd(_y\d_c\d)',z);assert m;suf=m.group(1)
  pro=f' cbz x19,pb_done{suf}\n mov w0,#-1\n mov x2,#0xffffffff00000001\n'
  assert z.count(pro)==1
  z=z.replace(pro,pro+' sub sp,sp,#48\n stp x5,x6,[sp]\n stp x7,x8,[sp,#16]\n str x9,[sp,#32]\n')
  epi=f'pb_done{suf}:\n ldr x30,[sp],#16'
  assert z.count(epi)==1
  z=z.replace(epi,f'pb_done{suf}:\n add sp,sp,#48\n ldr x30,[sp],#16')
  for o,sn,rg,k in [(' stp x7,x8,[x21]\n stp x5,x6,[x21, #16]\n','','x5,x6',1),(' stp x7,x8,[x21, #32]\n stp x3,x4,[x21, #48]\n',' orr x11,x11,#0x800000000\n','x3,x4',2)]:
   assert z.count(o)==1
   sk='fskip%d%s'%(k,suf)
   z=z.replace(o,f' ldp x9,x10,[sp]\n ldr x11,[x20,#64]\n'+sn+f' and x7,x7,#0x7ffffffff\n and x8,x8,#0x7ffffffff\n mov x12,#32\n add x13,x9,x7,lsl #3\n ldadd x12,x14,[x13]\n orr x15,x8,x11\n stp {rg},[x14]\n stp x10,x15,[x14,#16]\n tst x14,#0xfe0\n b.ne {sk}\n mov x16,x30\n bl flush_publish_b\n mov x30,x16\n{sk}:\n')
  out.append(z)
 out.append(FLUSH_SUB_B)
 fused='\n'.join(out)
 assert fused.count('\n bl flush_publish_b')==8
 assert fused.count('stp x7,x8,[x21')==0
 assert fused.count('sub sp,sp,#48')==4 and fused.count('add sp,sp,#48')==4
 return fused
FWD_LIMBS=['subs x11, x3, x7','sbcs x12, x4, x8','sbcs x13, x5, x9','sbcs x14, x6, x10']
REV_LIMBS=['subs x11, x7, x3','sbcs x12, x8, x4','sbcs x13, x9, x5','sbcs x14, x10, x6']
def assert_limb_specialization(base):
 fwd=base.split('.globl pair_bwd',1)[0]
 lo,hi=fwd.split('pair_fwd_hi',1)
 for line in REV_LIMBS:assert lo.count(line)==1,line
 for line in FWD_LIMBS:assert lo.count(line)==0,line
 for line in FWD_LIMBS:assert hi.count(line)==1,line
 for line in REV_LIMBS:assert hi.count(line)==0,line
 variants=re.split(r'(?=\.globl pair_bwd_y\d_c\d)',base)[1:]
 assert len(variants)==4
 for z in variants:
  want=REV_LIMBS if re.match(r'\.globl pair_bwd_y1',z) else FWD_LIMBS
  other=FWD_LIMBS if want is REV_LIMBS else REV_LIMBS
  for line in want:assert z.count(line)==1,(line,z[:40])
  for line in other:assert z.count(line)==0,(line,z[:40])
 return True
def check(root):
 import hashlib
 code=generated_assembly(root)
 generated_hash=hashlib.sha256(code.encode()).hexdigest()
 base=base_generated_assembly(root)
 assert hashlib.sha256(base.encode()).hexdigest()=='a4bf5e12e181237e0e1bc107b75a478f1adf8234737810d2264c44346618812e'
 assert assert_limb_specialization(base)
 assert len(code.encode())<1<<20
 import a64
 encoded_bytes=4*len(a64.a64_load(code)[0])
 assert encoded_bytes<(1<<18)
 helper_encoded_bytes=sum(4*len(a64.a64_load(f.read_text())[0]) for f in root.glob("*.s") if f.name!='block.s')
 blk=block_production(root);blk_enc=4*len(a64.a64_load(blk)[0])
 assert encoded_bytes+helper_encoded_bytes+blk_enc+65536<=1<<20
 return dict(exact_assembly_derivation=True,generated_sha256=generated_hash,forward_bound="453/4",forward_fixed_bound=45,fused_rolled_per_item=637,fused_rolled_fixed=16,fused_backward_per_item="2541/4",backward_per_carry=4,fused_backward_fixed=33,flush_taken=18,flush_dropped=15,generated_source_bytes=len(code.encode()),generated_encoded_bytes=encoded_bytes,helper_encoded_bytes=helper_encoded_bytes,count_delta="fwd 6-38*floor(n/8)+2*(n>=8)+(n%8>=1), bwd/emit 4-14*floor(n/8) for nonempty call")
