import re
import chain_witness
def base_generated_assembly(root):
 source,meta,fwd_variants=chain_witness.generate(root)
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
    z=z.replace('subs x11, x3, x7','subs x11, x7, x3')
    z=z.replace('sbcs x12, x4, x8','sbcs x12, x8, x4')
    z=z.replace('sbcs x13, x5, x9','sbcs x13, x9, x5')
    z=z.replace('sbcs x14, x6, x10','sbcs x14, x10, x6')
   old=' adc x15, xzr, xzr\n cbz x15,k_add_store_2\n';assert old in z;z=z.replace(old,'')
   if not carry:
    old=' adds x11,x11,#1\n sbcs x12,x12,x0\n sbcs x13,x13,xzr\n sbc x14,x14,x2\nk_add_store_2:';assert old in z;z=z.replace(old,'k_add_store_2:')
   for old,new in [(' ldp x3, x4, [x20]\n',' ldp x3, x4, [x20], #-72\n'),
                   (' ldp x5, x6, [x20, #16]\n',' ldp x5, x6, [x20, #88]\n'),
                   (' ldp x3, x4, [x20, #32]\n',' ldp x3, x4, [x20, #104]\n'),
                   (' ldp x5, x6, [x20, #48]\n',' ldp x5, x6, [x20, #120]\n'),
                   (' ldp x7, x8, [x21]\n',' ldp x7, x8, [x21], #-64\n'),
                   (' ldp x13, x14, [x21, #16]\n',' ldp x13, x14, [x21, #80]\n'),
                   (' stp x7,x8,[x21]\n',' stp x7,x8,[x21, #64]\n'),
                   (' stp x5,x6,[x21, #16]\n',' stp x5,x6,[x21, #80]\n'),
                   (' stp x7,x8,[x21, #32]\n',' stp x7,x8,[x21, #96]\n'),
                   (' stp x3,x4,[x21, #48]\n',' stp x3,x4,[x21, #112]\n')]:
    assert z.count(old)==1,(old,yl,carry)
    z=z.replace(old,new)
   for dead in [' sub x20,x20,#64\n',' sub x21,x21,#64\n']:
    assert z.count(dead)==1,dead
    z=z.replace(dead,'')
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
def fwd_alt_bodies(root):
 _,_,variants=chain_witness.generate(root)
 out={}
 for parity,tag in ((0,'e'),(1,'o')):
  lo,hi=chain_witness.split_fwd_directions(variants[parity],tag)
  for direction,text in (('lo',lo),('hi',hi)):
   text=text.replace('[x20], #64','[x20], #72').replace('[x20, #-48]','[x20, #-56]')
   assert '[x20], #64' not in text and '[x20, #-48]' not in text
   lp='pf_loop_'+direction;dp='pf_success_'+direction
   s=text.index(lp+':\n');e=text.index(dp+':\n',s)
   block=text[s+len(lp)+2:e]
   end=f' subs x19,x19,#1\n b.ne {lp}\n'
   assert block.endswith(end)
   body=block[:-len(end)]
   assert not re.search(r'\bx19\b',body)
   out[(direction,parity)]=body
 return out
def unroll8_alt(code,loop,done,bodies,glob):
 start=code.index(loop+':\n');stop=code.index(done+':\n',start)
 est=code.index('.globl '+glob)
 assert est<start
 block=code[start+len(loop)+2:stop]
 end=f' subs x19,x19,#1\n b.ne {loop}\n'
 assert block.endswith(end)
 pinned=block[:-len(end)]
 assert not re.search(r'\bx19\b',pinned)
 for p in (0,1):
  assert not re.search(r'\bx19\b',bodies[p])
  assert pinned.count('\n')-bodies[p].count('\n')==2,(p,)
 head=code[est:start]
 old=' ldp x9, x10, [x3, #16]\n'
 assert head.count(old)==1
 head=head.replace(old,old+' ldp x28, x29, [x21, #-32]\n ldp x17, x18, [x21, #-16]\n')
 t=loop+'_tail';u=loop+'_single';g=loop+'_eight'
 new=f' subs x19,x19,#8\n b.lo {t}\n{g}:\n'
 for k in range(8):
  b=bodies[k%2]
  for lab in re.findall(r'^([A-Za-z_]\w*):',b,re.M):b=re.sub(r'\b'+lab+r'\b',lab+f'_u{k}',b)
  new+=b
 new+=f' subs x19,x19,#8\n b.hs {g}\n{t}:\n adds x19,x19,#8\n cbz x19,{done}\n{u}:\n'+pinned+f' subs x19,x19,#1\n b.ne {u}\n'
 return code[:est]+head+new+code[stop:]
def _keep(lines):
 o=[];r={3:0,4:0,5:0}
 for x in lines:
  w=x.strip().split();op=w[0] if w else ''
  d=w[1].split(',')[0] if len(w)>1 and op not in ('stp','str','b','cbz','cbnz','cmp','tst','tbz','tbnz','ret') else ''
  assert d!='x6' and not re.search(r'\bx28\b',x)
  if d in ('x3','x4','x5') and op!='ldp':r[int(d[1])]=1
  if r[3]:x=re.sub(r'\bx3\b','x1',x)
  if r[4]:x=re.sub(r'\bx4\b','x17',x)
  if r[5]:x=re.sub(r'\bx5\b','x28',x)
  o.append(x)
 assert r[3] and r[4] and r[5];return o
def optimize_backward(code):
 line=' ldp x3, x4, [x23]\n';high=' ldp x5, x6, [x23, #16]\n';p=re.split(r'(?=\.globl pair_bwd_y\d_c\d)',code);o=[p[0]]
 for z in p[1:]:
  a=z.index(line);b=z.index(line,a+1);rest=z[b+len(line):];c=rest.index(high)
  assert not re.search(r'\bx(5|6|28)\b',rest[:c])
  o.append(z[:a]+'\n'.join(_keep(z[a:b].splitlines()))+'\n'+rest[:c]+rest[c+len(high):])
 return ''.join(o)
def generated_assembly(root):
 code=optimize_backward(fused_base(root))
 alt=fwd_alt_bodies(root)
 for direction in ['lo','hi']:
  code=unroll8_alt(code,'pf_loop_'+direction,'pf_success_'+direction,{0:alt[(direction,0)],1:alt[(direction,1)]},'pair_fwd_'+direction)
 for yl in [0,1]:
  for carry in [0,1]:
   code=unroll8(code,'pb_loop'+f'_y{yl}_c{carry}','pb_done'+f'_y{yl}_c{carry}')
 return code
def block_production(root):
 return unroll8((root/'block.s').read_text(),'pk_loop','pk_done')
FLUSH_SUB_B=""".globl flush_publish_b
flush_publish_b:
 and x7,x7,#0x7ffffffff
 ldr x9,[sp,#16]
 ldr x10,[sp,#24]
 ldr x11,[sp,#32]
 mov x12,#16352
 sub x14,x14,x12
 str x14,[x13]
 add x9,x9,x7,lsl #3
 mov x12,#16384
 ldadd x12,x13,[x9]
 cmp x13,x10
 b.hs fpb_skip
 mov x12,#512
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
  o1=' stp x7,x8,[x21, #64]\n stp x5,x6,[x21, #80]\n';o2=' stp x7,x8,[x21, #96]\n stp x3,x4,[x21, #112]\n'
  assert z.count(o1)==1 and z.count(o2)==1
  e1=f' ldp x13,x10,[sp]\n ldr x1,[x20,#136]\n bfi x13,x7,#3,#35\n mov x17,#32\n ldadd x17,x14,[x13]\n bfi x1,x8,#0,#35\n stp x5,x6,[x14]\n stp x10,x1,[x14,#16]\n tbz x14,#14,fskip1{suf}\n mov x16,x30\n bl flush_publish_b\n mov x30,x16\nfskip1{suf}:\n'
  e2=f' ldp x13,x10,[sp]\n orr x15,x1,#0x800000000\n bfi x13,x7,#3,#35\n ldadd x17,x14,[x13]\n bfi x15,x8,#0,#35\n stp x3,x4,[x14]\n stp x10,x15,[x14,#16]\n tbz x14,#14,fskip2{suf}\n mov x16,x30\n bl flush_publish_b\n mov x30,x16\nfskip2{suf}:\n'
  z=z.replace(o1,e1).replace(o2,e2)
  # b19: the ID word (x1) and the ldadd addend 32 (x17) are carried from the first emission to the
  # second; the intervening slope product (index 8) never touches x1/x17 and the fused square
  # (index 9) has its x1/x17 scratch renamed to the dead inverse-bank registers x24/x25.
  a=z.index(f'fskip1{suf}:\n')+len(f'fskip1{suf}:\n');m=z.index(f'k_mul_store_8{suf}:\n',a);b=z.index(e2,m)
  assert not re.search(r'\b[xw](1|17)\b',z[a:m]) and not re.search(r'\b[xw](24|25)\b',z[m:b])
  mid=re.sub(r'\bx1\b','x24',z[m:b]);mid=re.sub(r'\bx17\b','x25',mid)
  assert not re.search(r'\b[xw](1|17)\b',mid)
  z=z[:m]+mid+z[b:]
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
 assert hashlib.sha256(base.encode()).hexdigest()=='07f98f57a417d080a5ccdf49fc8d1381b501b174251f89606e5ea69a48b95a29'
 assert assert_limb_specialization(base)
 assert len(code.encode())<1<<20
 import a64
 encoded_bytes=4*len(a64.a64_load(code)[0])
 assert encoded_bytes<(1<<18)
 helper_encoded_bytes=sum(4*len(a64.a64_load(f.read_text())[0]) for f in root.glob("*.s") if f.name!='block.s')
 blk=block_production(root);blk_enc=4*len(a64.a64_load(blk)[0])
 assert encoded_bytes+helper_encoded_bytes+blk_enc+65536<=1<<20
 return dict(exact_assembly_derivation=True,generated_sha256=generated_hash,forward_bound="441/4",forward_fixed_bound=46,fused_rolled_per_item=626,fused_rolled_fixed=16,fused_backward_per_item="2493/4",backward_per_carry=4,fused_backward_fixed=33,flush_taken=19,flush_dropped=16,generated_source_bytes=len(code.encode()),generated_encoded_bytes=encoded_bytes,helper_encoded_bytes=helper_encoded_bytes,count_delta="6-30*floor(n/8) fwd (ping-pong mains); 4-14*floor(n/8)-2n bwd/emit (b16 low-pair and b18 high-pair first-operand reuse; b17 TBZ/BFI emission and b19 carried ID/addend are in both rolled and unrolled fused bodies)")
