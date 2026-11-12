"""AArch64 interpreter/counter (one step per executed instruction): scalar GPR ops, W regs, BL/RET, device stores, reg-offset LDR/STR, LDADD, SWP, DC ZVA, LD64B/ST64B.
Checks only: LS64 operand form, alignment, permitted regions and base outside Xt..Xt+7; DC ZVA alignment; unique named labels; unknown mnemonics raise;
TST with a shifted register and W-register forms other than MOV/MOVK/CBZ/CBNZ/TBZ/TBNZ raise. It does NOT validate operand shapes in general and does not check
A64 encodability (immediate ranges, register classes): probe.py restricts the shipped listings to flip_check.FORMS and clang assembles every listing."""
M64=(1<<64)-1
def expand(text):
 macros,out,cur={},[],None
 for raw in text.split("\n"):
  s=raw.split("//")[0].strip()
  if s.startswith(".macro"):
   parts=s[len(".macro"):].replace(","," ").split()
   cur=macros[parts[0]]=(parts[1:],[])
   continue
  if s.startswith(".endm"):
   cur=None
   continue
  (cur[1]if cur is not None else out).append(raw)
 def emit(lines,depth):
  assert depth<16
  res=[]
  for raw in lines:
   s=raw.split("//")[0].strip()
   name=s.split()[0]if s else ""
   if raw[:1].isspace()and name in macros:
    params,body=macros[name]
    args=[a.strip()for a in s[len(name):].split(",")]if s[len(name):].strip()else[]
    sub=[]
    for b in body:
     for p_,a in zip(params,args):
      b=b.replace("\\"+p_,a)
     sub.append(b)
    res+=emit(sub,depth+1)
   else:
    res.append(raw)
  return res
 return "\n".join(emit(out,0))
def _a64_ops(s):
 ops,cur,depth=[],"",0
 for ch in s:
  if ch=="[":
   depth+=1
  if ch=="]":
   depth-=1
  if ch=="," and depth==0:
   ops.append(cur.strip())
   cur=""
  else:
   cur+=ch
 if cur.strip():
  ops.append(cur.strip())
 return ops
def a64_load(text):
 prog,labels=[],{}
 for raw in expand(text).split("\n"):
  line=raw.split("//")[0].rstrip()
  if not line.strip():
   continue
  if not line[0].isspace():
   labels.setdefault(line.rstrip(":"),[]).append(len(prog))
   continue
  ins=line.strip()
  if ins.startswith("."):
   continue
  parts=ins.split(None,1)
  mn,rest=parts[0],parts[1]if len(parts)>1 else ""
  ops=_a64_ops(rest)
  assert mn!="tst" or len(ops)==2,"TST shifted-register form not implemented"
  assert mn in("mov","movk","cbz","cbnz","tbz","tbnz")or not any(x=="wzr" or x[:1]=="w" and x[1:].isdigit()for x in ops),"W-register form not implemented: "+ins
  prog.append((mn,ops))
 assert all(len(v)==1 for k,v in labels.items()if not k.isdigit()),"duplicate label"
 return prog,labels
def a64_run(prog,labels,entry,regs,mem,limit=10**7,out_regs=None,visits=None,dev=None,ls64=([],[])):
 # ls64: (LD64B, ST64B) permitted (lo,hi) regions (proof.md section4)
 M=M64
 R={("x%d"%i):0 for i in range(31)}
 R["sp"]=0
 R.update(regs)
 lr0=R["x30"]
 flags={"C":0,"Z":0}
 pc=labels[entry][0]
 steps=0
 depth=0
 def rd(o):
  if o in("xzr","wzr"):
   return 0
  if o[0]=="w":
   return R["x"+o[1:]]&0xFFFFFFFF
  return R[o]
 def wr(o,v):
  if o in("xzr","wzr"):
   return
  if o[0]=="w":
   R["x"+o[1:]]=v&0xFFFFFFFF
  else:
   R[o]=v&M
 def imm(o):
  return int(o.lstrip("#"),0)
 def ea(op,post,regoff=False):
  inner=op.strip()
  pre=inner.endswith("!")
  inner=inner.rstrip("!")[1:-1]
  parts=[p.strip()for p in inner.split(",")]
  base=parts[0]
  if len(parts)>1 and parts[1][:1]=="x":
   assert regoff and post is None and not pre,"register offset only for plain LDR/STR"
   sh=parts[2].split()if len(parts)>2 else["lsl","#0"]
   assert len(parts)<=3 and sh[0]=="lsl" and imm(sh[1])in(0,3)
   return(rd(base)+(rd(parts[1])<<imm(sh[1])))&M,None,base
  off=imm(parts[1])if len(parts)>1 else 0
  if post is not None:
   return rd(base),(rd(base)+imm(post))&M,base
  if pre:
   return(rd(base)+off)&M,(rd(base)+off)&M,base
  return(rd(base)+off)&M,None,base
 def cond(c):
  return{"eq":flags["Z"]==1,"ne":flags["Z"]==0,"hs":flags["C"]==1,"cs":flags["C"]==1,
    "lo":flags["C"]==0,"cc":flags["C"]==0}[c]
 def store_dev(addr,v):
  if dev[0]<=addr<dev[0]+dev[1]:
   dev[2](addr-dev[0],v)
 def tgt(t):
  if t[:-1].isdigit():
   c=labels[t[:-1]]
   return[p for p in c if p<pc][-1]if t[-1]=="b" else[p for p in c if p>=pc][0]
  return labels[t][0]
 while True:
  mn,o=prog[pc]
  steps+=1
  assert steps<limit
  if visits is not None and pc in visits:
   visits[pc]+=1
  pc+=1
  if mn in("mul","umulh"):
   a,b=rd(o[1]),rd(o[2])
   wr(o[0],(a*b)&M if mn=="mul" else(a*b)>>64)
  elif mn in("adds","adc","adcs","add")and not(mn=="add" and(o[2].startswith("#")or len(o)==4)):
   a=rd(o[1])
   b=imm(o[2])if o[2].startswith("#")else rd(o[2])
   if len(o)==4:
    sh=o[3].split()
    assert mn=="adds" and sh[0]in("lsl","lsr")
    b=((b<<imm(sh[1]))&M)if sh[0]=="lsl" else(b>>imm(sh[1]))
   cin=flags["C"]if mn in("adc","adcs")else 0
   t=a+b+cin
   wr(o[0],t)
   if mn in("adds","adcs"):
    flags["C"]=t>>64
    flags["Z"]=1 if(t&M)==0 else 0
  elif mn=="add":
   if o[2].startswith("#"):
    wr(o[0],rd(o[1])+imm(o[2]))
   else:
    sh=o[3].split();n=imm(sh[1]);assert sh[0]in("lsl","lsr")and 0<=n<64
    wr(o[0],rd(o[1])+(((rd(o[2])<<n)&M)if sh[0]=="lsl" else rd(o[2])>>n))
  elif mn in("sub","subs","sbc","sbcs"):
   a=rd(o[1])
   b=imm(o[2])if o[2].startswith("#")else rd(o[2])
   if len(o)==4:
    sh=o[3].split()
    assert mn in("sub","subs")and sh[0]in("lsl","lsr")
    b=((b<<imm(sh[1]))&M)if sh[0]=="lsl" else(b>>imm(sh[1]))
   cin=flags["C"]if mn in("sbc","sbcs")else 1
   t=a+((~b)&M)+cin
   wr(o[0],t)
   if mn in("subs","sbcs"):
    flags["C"]=t>>64
    flags["Z"]=1 if(t&M)==0 else 0
  elif mn=="csel":
   wr(o[0],rd(o[1])if cond(o[3])else rd(o[2]))
  elif mn=="madd":
   wr(o[0],rd(o[1])*rd(o[2])+rd(o[3]))
  elif mn in("cmp","cmn","ccmp","ccmn"):
   if mn[1]=="c" and not cond(o[3]):
    flags["C"],flags["Z"]=imm(o[2])>>1&1,imm(o[2])>>2&1
   else:
    a,b=rd(o[0]),imm(o[1])if o[1].startswith("#")else rd(o[1])
    if len(o)==3:
     assert mn=="cmp" and o[2][:4]=="asr "
     b=((b-(b>>63<<64))>>imm(o[2][4:]))&M
    t=a+b if mn[-1]=="n" else a+((~b)&M)+1
    flags["C"],flags["Z"]=t>>64,1 if(t&M)==0 else 0
  elif mn=="tst":
   t=rd(o[0])&(imm(o[1])if o[1].startswith("#")else rd(o[1]))
   flags["C"],flags["Z"]=0,1 if t==0 else 0
  elif mn in("orr","and","eor"):
   b=imm(o[2])if o[2].startswith("#")else rd(o[2])
   assert len(o)in(3,4)
   if len(o)==4:
    assert not o[2].startswith("#")
    sh=o[3].split();n=imm(sh[1]);assert sh[0]in("lsl","lsr")and 0<=n<64
    b=((b<<n)&M)if sh[0]=="lsl" else b>>n
   wr(o[0],rd(o[1])|b if mn=="orr" else rd(o[1])&b if mn=="and" else rd(o[1])^b)
  elif mn=="ubfx":
   lo,width=imm(o[2]),imm(o[3]);assert 0<=lo<64 and 1<=width<=64-lo
   wr(o[0],(rd(o[1])>>lo)&((1<<width)-1))
  elif mn=="lsl":
   wr(o[0],rd(o[1])<<imm(o[2]))
  elif mn=="cinc":
   wr(o[0],rd(o[1])+(1 if cond(o[2])else 0))
  elif mn=="cset":
   wr(o[0],1 if cond(o[1])else 0)
  elif mn=="bfi":
   mask=((1<<imm(o[3]))-1)<<imm(o[2])
   wr(o[0],(rd(o[0])&~mask)|((rd(o[1])<<imm(o[2]))&mask))
  elif mn=="extr":
   wr(o[0],((rd(o[1])<<64|rd(o[2]))>>imm(o[3])))
  elif mn=="lsr":
   wr(o[0],rd(o[1])>>imm(o[2]))
  elif mn=="mov":
   wr(o[0],imm(o[1])if o[1].startswith("#")else rd(o[1]))
  elif mn=="movk":
   sh=imm(o[2].split()[1])if len(o)>2 else 0
   wr(o[0],(rd(o[0])&~(0xFFFF<<sh))|((imm(o[1])&0xFFFF)<<sh))
  elif mn in("ldp","stp"):
   post=o[3]if len(o)>3 else None
   addr,newb,base=ea(o[2],post)
   if mn=="ldp":
    wr(o[0],mem[addr])
    wr(o[1],mem[addr+8])
   else:
    mem[addr]=rd(o[0])
    mem[addr+8]=rd(o[1])
    if dev is not None:
     store_dev(addr,mem[addr])
     store_dev(addr+8,mem[addr+8])
   if newb is not None:
    R[base]=newb
  elif mn in("ldr","str","stlr"):
   if mn=="stlr":
    assert len(o)==2 and o[1].startswith("[")and o[1].endswith("]")and "," not in o[1]and "!" not in o[1]
   post=o[2]if len(o)>2 else None
   addr,newb,base=ea(o[1],post,regoff=True)
   if mn=="ldr":
    wr(o[0],mem[addr])
   else:
    mem[addr]=rd(o[0])
    if dev is not None:
     store_dev(addr,mem[addr])
   if newb is not None:
    R[base]=newb
  elif mn=="stadd":
   addr,newb,base=ea(o[1],None)
   mem[addr]=(mem[addr]+rd(o[0]))&M
  elif mn=="swp":
   assert o[2][0]+o[2][-1]=="[]"and","not in o[2];a=rd(o[2][1:-1]);v=mem[a];mem[a]=rd(o[0]);wr(o[1],v)
  elif mn=="ldadd":
   addr,newb,base=ea(o[2],None)
   old=mem[addr]
   mem[addr]=(old+rd(o[0]))&M
   wr(o[1],old)
  elif mn in("ld64b","st64b"):
   n,b=int(o[0][1:]),o[1]
   assert o[0][0]=="x" and n%2==0 and n<=22 and len(o)==2 and b[0]+b[-1]=="[]" and "," not in b,"LS64 Xt/form"
   assert b[1:-1].strip()not in["x%d"%(n+k)for k in range(8)],"LS64 base in Xt..Xt+7"
   addr=rd(b[1:-1].strip());assert addr%64==0 and any(lo<=addr<=hi-64 for lo,hi in ls64[mn=="st64b"]),"LS64 alignment/region"
   for k in range(8):
    if mn=="ld64b":wr("x%d"%(n+k),mem[addr+8*k])
    else:
     mem[addr+8*k]=rd("x%d"%(n+k))
     if dev is not None:store_dev(addr+8*k,mem[addr+8*k])
  elif mn=="dc":
   assert o[0]=="zva" and len(o)==2 and rd(o[1])%256==0,"DC ZVA: aligned 256-byte block (DCZID_EL0.BS=6)"
   addr=rd(o[1])
   for k in range(32):
    mem[addr+8*k]=0
  elif mn in("b","b.eq","b.ne","b.lo","b.hs","b.cs","b.cc"):
   if mn=="b" or cond(mn[2:]):
    pc=tgt(o[0])
  elif mn in("cbz","cbnz","tbz","tbnz"):
   v=rd(o[0])if mn[0]=="c" else rd(o[0])>>imm(o[1])&1
   if(v==0)==(mn in("cbz","tbz")):
    pc=tgt(o[-1])
  elif mn=="bl":
   R["x30"]=pc
   depth+=1
   pc=tgt(o[0])
  elif mn=="ret":
   if depth:
    depth-=1
    pc=R["x30"]
    continue
   assert R["x30"]==lr0,"outermost RET with a clobbered x30"
   if out_regs is not None:
    out_regs.update(R)
   return steps
  else:
   raise ValueError("unsupported instruction: "+mn)
