import a64

def check(root):
 pr,lab=a64.a64_load((root/'pair_keys.s').read_text());M=(1<<64)-1;mask=(1<<35)-1;cases=0
 for row in [0,mask]:
  for index in [0,(1<<29)-1,377461701]:
   for leaf in [0,(1<<28)-1]:
    mem={};base=0x10000;port=0x20000;dst=0x40000;packets=[];packet=[0]*6
    keys=[mask|(mask<<64)|(7<<128)|(9<<192),3|(2<<64)|(11<<128)|(13<<192)]
    for j,key in enumerate(keys):
     for k in range(4):mem[base+32*j+8*k]=(key>>(64*k))&M
     mem[port+64+8*(key&mask)]=dst
    def device(off,val):
     packet[off//8]=val
     if off==40:packets.append(tuple(packet))
    out={};w2=row|(index<<35)
    a64.a64_run(pr,lab,'pk_loop',dict(x19=1,x20=base,x21=port,x22=port+64,x23=w2,x24=leaf<<36,x25=dst+64,x26=32,x27=1<<35,x28=0xffffffff00000000,x30=12345),mem,limit=1000,out_regs=out,dev=(port,48,device))
    assert len(packets)==2
    for sign,v in enumerate(packets):
     assert v[4]==w2 and v[4]&mask==row and v[4]>>35==index
     assert v[5]&mask==(keys[sign]>>64)&mask
     assert (v[5]>>35)&1==sign and v[5]>>36==leaf
    assert out['x24']==((leaf+1)<<36)&M
    cases+=1
 return dict(cases=cases,production_width_tags_roundtrip=True,last_index_wrap_is_post_emission=True)
