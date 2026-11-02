import re
import chain_witness
def generated_assembly(root):
 source,meta=chain_witness.generate(root)
 source=chain_witness.normalize_transform(source)
 t=source;f,b=t.split('.globl pair_bwd',1)
 f=f.replace('[x20], #64','[x20], #72').replace('[x20, #-48]','[x20, #-56]')
 # Whitespace exact forms from frozen source.
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
    old=' adds x11,x11,#1\n sbcs x12,x12,x0\n sbcs x13,x13,xzr\n adc x14,x14,x2\nk_add_store_2:';assert old in z;z=z.replace(old,'k_add_store_2:')
   z=z.replace('sub x20,x20,#64','sub x20,x20,#72')
   suffix=f'_y{yl}_c{carry}'
   for lab in re.findall(r'^([A-Za-z_]\w*):',z,re.M):z=re.sub(r'\b'+lab+r'\b',lab+suffix,z)
   parts.append(z)
 return '\n'.join(parts)
def check(root):
 import hashlib
 code=generated_assembly(root)
 assert hashlib.sha256(code.encode()).hexdigest()=='f1924d38228e8367840244af3d58810ea61aead4f066f99de62f599ee2e4dacd'
 assert len(code.encode())<1<<20
 t=(root/'keys_reference.s').read_text().replace(' mov x24,#0',' mov x28,x5').replace(' mov x27,#0x800000000\n','').replace('pk_loop:\n','pk_loop:\n ldr x24,[x28],#72\n').replace(' add x24,x24,x27\n','')
 assert t==(root/'range_keys.s').read_text()
 return dict(exact_assembly_derivation=True,Y_stride=72,batch_stride=64,forward_bound=121,backward_bound=632,emitter_bound=29,variant_count=4)
