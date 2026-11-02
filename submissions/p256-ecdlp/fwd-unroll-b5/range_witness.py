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
    z=z.replace('subs x11, x3, x7','subs x11, x7, x3').replace('sbcs x12, x4, x8','sbcs x12, x8, x4').replace('sbcs x13, x5, x9','sbcs x13, x9, x5').replace('sbcs x14, x6, x10','sbcs x14, x10, x6')
   old=' adc x15, xzr, xzr\n cbz x15,k_add_store_2\n';assert old in z;z=z.replace(old,'')
   if not carry:
    old=' adds x11,x11,#1\n sbcs x12,x12,x0\n sbcs x13,x13,xzr\n sbc x14,x14,x2\nk_add_store_2:';assert old in z;z=z.replace(old,'k_add_store_2:')
   z=z.replace('sub x20,x20,#64','sub x20,x20,#72')
   suffix=f'_y{yl}_c{carry}'
   for lab in re.findall(r'^([A-Za-z_]\w*):',z,re.M):z=re.sub(r'\b'+lab+r'\b',lab+suffix,z)
   parts.append(z)
 return '\n'.join(parts)
def generated_assembly(root):
 code=base_generated_assembly(root)
 pairs=[('pf_loop_lo','pf_success_lo'),('pf_loop_hi','pf_success_hi')]
 for yl in [0,1]:
  for carry in [0,1]:
   pairs.append(('pb_loop'+f'_y{yl}_c{carry}','pb_done'+f'_y{yl}_c{carry}'))
 for loop,done in pairs:
   start=code.index(loop+':\n');stop=code.index(done+':\n',start)
   block=code[start+len(loop)+2:stop]
   end=f' subs x19,x19,#1\n b.ne {loop}\n'
   assert block.endswith(end)
   body=block[:-len(end)]
   assert not re.search(r'\bx19\b',body),'Loop body unexpectedly uses count'
   tail=loop+'_tail';single=loop+'_single';group=loop+'_eight'
   new=f' subs x19,x19,#8\n b.lo {tail}\n{group}:\n'
   for k in range(8):
    b=body
    for lab in re.findall(r'^([A-Za-z_]\w*):',body,re.M):b=re.sub(r'\b'+lab+r'\b',lab+f'_u{k}',b)
    new+=b
   new+=f' subs x19,x19,#8\n b.hs {group}\n{tail}:\n adds x19,x19,#8\n cbz x19,{done}\n{single}:\n'+body+f' subs x19,x19,#1\n b.ne {single}\n'
   code=code[:start]+new+code[stop:]
 return code
def check(root):
 import hashlib
 code=generated_assembly(root)
 generated_hash=hashlib.sha256(code.encode()).hexdigest()
 assert hashlib.sha256(base_generated_assembly(root).encode()).hexdigest()=='1e5508db4b8c2b8742a43337e0ab3f47fb654a4dc7c94151e764761b438e192a'
 assert len(code.encode())<1<<20
 import a64
 encoded_bytes=4*len(a64.a64_load(code)[0])
 assert encoded_bytes<(1<<18)
 helper_encoded_bytes=sum(4*len(a64.a64_load(f.read_text())[0]) for f in root.glob("*.s"))
 assert encoded_bytes+helper_encoded_bytes+65536<=1<<20
 t=(root/'keys_reference.s').read_text().replace(' mov x24,#0',' mov x28,x5').replace(' mov x27,#0x800000000\n','').replace('pk_loop:\n','pk_loop:\n ldr x24,[x28],#72\n').replace(' add x24,x24,x27\n','')
 assert t==(root/'range_keys.s').read_text()
 return dict(exact_assembly_derivation=True,generated_sha256=generated_hash,Y_stride=72,batch_stride=64,forward_bound="481/4",forward_fixed_bound=30,backward_base_per_item="2457/4",backward_per_carry=4,backward_fixed_bound=28,emitter_bound=29,variant_count=4,generated_source_bytes=len(code.encode()),generated_encoded_bytes=encoded_bytes,helper_encoded_bytes=helper_encoded_bytes,unroll_factor=8,count_delta="4-14*floor(n/8) for nonempty fwd/bwd call")
