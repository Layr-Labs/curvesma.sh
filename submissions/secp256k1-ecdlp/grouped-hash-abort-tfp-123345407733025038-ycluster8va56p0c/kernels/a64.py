import re as _re
M64=(1<<64)-1
def va56(addr,n=8):
 assert 0<=addr and addr+n<=1<<64 and(addr+n<=1<<56 or addr>=(1<<64)-(1<<56)),"data address 0x%x outside the 56-bit VA ranges"%addr
 return addr
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
 count=[0]
 def emit(lines,depth):
  assert depth<16
  res=[]
  for raw in lines:
   s=raw.split("//")[0].strip()
   name=s.split()[0]if s else""
   if not s.endswith(":")and name in macros:
    params,body=macros[name]
    args=[a.strip()for a in s[len(name):].split(",")]if s[len(name):].strip()else[]
    count[0]+=1
    sub=[]
    for b in body:
     for p_,a in zip(params,args):
      b=b.replace("\\"+p_,a)
     sub.append(b.replace("\\@",str(count[0])))
    res+=emit(sub,depth+1)
   else:
    res.append(raw)
  return res
 return"\n".join(emit(out,0))
def _a64_ops(s):
 ops,cur,depth=[],"",0
 for ch in s:
  if ch=="[":depth+=1
  if ch=="]":depth-=1
  if ch==","and depth==0:
   ops.append(cur.strip());cur=""
  else:
   cur+=ch
 if cur.strip():ops.append(cur.strip())
 return ops
def _int(o):
 m=_re.fullmatch(r"#?(-?(0x[0-9a-fA-F]+|0b[01]+|[0-9]+))",o.strip())
 if not m:raise ValueError("bad immediate: "+o)
 return int(m.group(1),0)
def _xreg(o):
 return _re.fullmatch(r"x([12]?[0-9]|30)|sp|xzr",o.strip(),_re.ASCII)is not None and not _re.fullmatch(r"x0[0-9]",o.strip())
def _bitmask64(v):
 if not-(1<<63)<=v<1<<64:return False
 v&=M64
 if v in(0,M64):return False
 for e in(2,4,8,16,32,64):
  m=(1<<e)-1;el=v&m
  if any((v>>k)&m!=el for k in range(0,64,e)):continue
  if el in(0,m):continue
  if sum(1 for i in range(e)if el>>i&1 and not el>>((i-1)%e)&1)==1:return True
 return False
def _arith_imm(v):
 v=abs(v)
 return 0<=v<4096 or(v%4096==0 and v>>12<4096)
def _mov_imm(v):
 if not-(1<<63)<=v<1<<64:return False
 v&=M64
 return any(v& ~(0xFFFF<<s)&M64==0 for s in(0,16,32,48))or any(~v&M64& ~(0xFFFF<<s)&M64==0 for s in(0,16,32,48))or _bitmask64(v)
def _mem(op):
 op=op.strip();pre=op.endswith("!");op=op[: -1]if pre else op
 if not(op[:1]=="["and op[-1:]=="]"):raise AssertionError("not a memory operand: "+op)
 parts=[q.strip()for q in op[1: -1].split(",")]
 if not _xreg(parts[0])or parts[0]=="xzr":raise AssertionError("bad base register: "+op)
 if pre and len(parts)==1:raise AssertionError("not encodable: "+op)
 if len(parts)>1 and not parts[1].startswith("#"):return parts[0],0,pre,parts[1:]
 if len(parts)>2:raise AssertionError("not encodable: "+op)
 try:return parts[0],_int(parts[1])if len(parts)>1 else 0,pre,None
 except ValueError:raise AssertionError("not encodable: "+op)
def _regs_ok(mn,ops):
 ops=["xzr"if o.strip()=="wzr"else o for o in ops]
 regs=[o.strip()for o in ops if _xreg(o.strip())]
 if mn in("add","sub","adds","subs","cmp","cmn"):
  k=1 if mn in("cmp","cmn")else 2
  dst,srcs=([]if k==1 else[ops[0]]),ops[k-1:k]
  imm=len(ops)>k and ops[k].startswith("#")
  ext=imm or(len(ops)<=k+2 and"sp"in ops[:k]and(len(ops)==k+1 or ops[k+1].split()[0]=="lsl"))
  if not ext:return"sp"not in ops
  if"sp"in ops[k:]:return False
  if mn in("adds","subs","cmp","cmn"):return"xzr"not in srcs and"sp"not in dst
  return"xzr"not in ops[:k]
 if mn in("and","orr","eor")and len(ops)==3 and ops[2].startswith("#"):
  return ops[1]!="sp"and ops[0]!="xzr"
 if mn=="mov":
  if ops[1].startswith("#"):
   v=_int(ops[1])&M64
   wide=any(v& ~(0xFFFF<<q)&M64==0 or~v&M64& ~(0xFFFF<<q)&M64==0 for q in(0,16,32,48))
   return wide if ops[0]=="xzr"else _bitmask64(v)if ops[0]=="sp"else True
  return not("sp"in ops and"xzr"in ops)
 return"sp"not in regs or mn in("ldp","stp","ldr","str","stlur","stlr","ld64b","st64b","stadd","staddl","ldadd","swp","ldaxr")
_COND={"eq","ne","hs","cs","lo","cc","mi","pl","vs","vc","hi","ls","ge","lt","gt","le"}
_DP3=["rrr"]
_ARITH=["rri","rrr","rrrs"]
_SHAPES=dict(mul=_DP3,umulh=_DP3,adc=_DP3,adcs=_DP3,sbc=_DP3,sbcs=_DP3,madd=["rrrr"],csel=["rrrc"],cset=["rc"],cinc=["rrc"],mvn=["rr","rrs"],lsl=["rri","rrr"],lsr=["rri","rrr"],ubfx=["rrii"],bfi=["rrii"],bfxil=["rrii"],extr=["rrri"],ccmp=["rric","riic"],ccmn=["rric","riic"],add=_ARITH,adds=_ARITH,sub=_ARITH,subs=_ARITH,cmp=["ri","rr","rrs"],cmn=["ri","rr","rrs"],tst=["ri","rr","rrs"],mov=["rr","ri"],movk=["ri","ris"],movn=["ri","ris"],movz=["ri","ris"],ldp=["rrm","rrmi"],stp=["rrm","rrmi"],ldr=["rm","rmi"],str=["rm","rmi"],stlur=["rm"],stlr=["rm"],ld64b=["rm"],st64b=["rm"],stadd=["rm"],staddl=["rm"],ldadd=["rrm"],swp=["rrm"],ldaxr=["rm"],blr=["r"],wfe=[""],ret=["","r"],b=["l"],bl=["l"],cbz=["rl"],cbnz=["rl"],tbz=["ril"],tbnz=["ril"],**{k:_ARITH for k in("and","orr","eor","ands","bic","eon")},**{"b."+c:["l"]for c in("eq","ne","lo","hs","cs","cc","lt","ge")})
_W_OK={"mul","adc","adcs","sbc","sbcs","madd","csel","cset","cinc","mvn","lsl","lsr","ubfx","bfi","bfxil","extr","add","adds","sub","subs","cmp","cmn","tst","mov","and","orr","eor","ands","bic","eon"}
def _kind(o):
 o=o.strip()
 if _xreg(o):return"r"
 if _re.fullmatch(r"w([12]?[0-9]|30)|wzr",o):return"w"
 if o.startswith("#"):
  try:_int(o);return"i"
  except ValueError:return"?"
 if o.startswith("["):return"m"
 if _re.fullmatch(r"(lsl|lsr|asr|ror) #[0-9]+",o):return"s"
 if o in _COND:return"c"
 if _re.fullmatch(r"[xwqdshbp][0-9]+|[vz][0-9]+(\.[0-9]*[bhsdq])?|zt[0-9]+|ffr|lr|fp|wsp|ip[01]|nzcv|al|nv",o):return"?"
 if _re.fullmatch(r"[A-Za-z_.][A-Za-z0-9_.]*|[0-9]+[fb]",o):return"l"
 return"?"
def _shape_ok(mn,ops):
 kinds=[_kind(o)for o in ops]
 if"w"in kinds:
  if mn not in _W_OK or"r"in kinds:return False
  kinds=["r"if k=="w"else k for k in kinds]
  if any(k=="i"for k in kinds)and mn not in("add","adds","sub","subs","cmp","cmn","lsl","lsr","ubfx","bfi","bfxil","extr"):
   return False
  if any(k=="s"and int(o.split("#")[1])>31 for k,o in zip(kinds,ops)):return False
 return"".join(kinds)in _SHAPES[mn]
_OK,_IM,_EA,_OP=set(),{},{},{}
def encodable(mn,ops):
 k=(mn,*ops)
 if k in _OK:return
 try:_encodable(mn,ops)
 except ValueError as ex:raise AssertionError("not encodable (bad immediate): %s %s (%s)"%(mn,",".join(ops),ex))
 _OK.add(k)
def _encodable(mn,ops):
 def bad(w=""):raise AssertionError("not encodable%s: %s %s"%(w,mn,",".join(ops)))
 if mn=="dmb":
  if ops!=["ish"]:bad()
  return
 if mn=="dc":
  if not(len(ops)==2 and ops[0]=="zva"and _xreg(ops[1])and ops[1]not in("sp","xzr")):bad()
  return
 if mn not in _SHAPES or not _shape_ok(mn,ops):bad()
 if any(_kind(o)=="w"for o in ops)and any(_kind(o)=="i"for o in ops):
  imms=[_int(o)for o in ops if _kind(o)=="i"]
  if mn in("lsl","lsr","extr")and not all(0<=v<=31 for v in imms):bad()
  if mn in("ubfx","bfi","bfxil")and not(0<=imms[0]<=31 and 1<=imms[1]<=32-imms[0]):bad()
 if not _regs_ok(mn,ops):bad()
 if mn in("ldr","str","ldp","stp","stlur","stlr","ld64b","st64b","ldaxr")and ops[0]=="sp":bad()
 if mn in("stadd","staddl","ldadd","swp")and"sp"in ops[: -1]:bad()
 if mn in("tbz","tbnz")and not 0<=_int(ops[1])<=63:bad()
 for o in ops:
  if _kind(o)=="s":
   kind,amt=o.split()[0],_int(o.split()[1])
   if not 0<=amt<=63 or(mn in("movk","movn","movz")and kind!="lsl")or(mn in("add","adds","sub","subs","cmp","cmn")and kind=="ror"):bad()
 if mn in("cbz","cbnz","tbz","tbnz")and ops[0]=="sp":bad()
 if mn in("ldp","stp"):
  if not(3<=len(ops)<=4 and _xreg(ops[0])and _xreg(ops[1])and"sp"not in ops[:2]):bad()
  base,off,pre,roff=_mem(ops[2])
  if roff is not None:bad()
  if len(ops)==4:
   if pre or off:bad()
   off=_int(ops[3])
  if off%8 or not-512<=off<=504 or((pre or len(ops)==4)and base!="sp"and base in ops[:2])or(mn=="ldp"and ops[0]==ops[1]):bad()
 elif mn in("ldr","str"):
  if not(2<=len(ops)<=3 and _xreg(ops[0])and ops[0]!="sp"):bad()
  base,off,pre,roff=_mem(ops[1])
  if roff is not None:
   if pre or len(ops)==3 or not _xreg(roff[0])or roff[0]=="sp"or len(roff)>2 or(len(roff)==2 and roff[1].replace(" ","")not in("lsl#0","lsl#3")):bad()
   return
  if len(ops)==3:
   if pre or off:bad()
   off=_int(ops[2])
  if pre or len(ops)==3:
   if not-256<=off<=255 or(base==ops[0]and base!="xzr"):bad()
  elif not((off%8==0 and 0<=off<=32760)or-256<=off<=255):bad()
 elif mn=="stlur":
  base,off,pre,roff=_mem(ops[1])
  if pre or roff is not None or not-256<=off<=255:bad()
 elif mn in("stlr","ldaxr","ld64b","st64b","stadd","staddl","ldadd","swp"):
  base,off,pre,roff=_mem(ops[-1])
  if pre or roff is not None or","in ops[-1]:bad()
  if mn in("ld64b","st64b")and not(ops[0][:1]=="x"and ops[0][1:].isdigit()and int(ops[0][1:])%2==0 and int(ops[0][1:])<=22):bad(" LS64")
 elif mn in("add","adds","sub","subs","cmp","cmn"):
  k=1 if mn in("cmp","cmn")else 2
  if len(ops)==k+1 and ops[k].startswith("#"):
   if not _arith_imm(_int(ops[k])):bad()
  elif len(ops)==k+2:
   if ops[k].startswith("#"):bad()
   kind,amt=ops[k+1].split();lim=4 if"sp"in ops[:k+1]else 63
   if kind not in("lsl","lsr","asr")or not 0<=_int(amt)<=lim or("sp"in ops[:k+1]and kind!="lsl"):bad()
 elif mn in("and","orr","eor","tst","ands","bic","eon"):
  k=1 if mn=="tst"else 2
  if ops[k].startswith("#"):
   if mn in("bic","eon")or len(ops)!=k+1 or not _bitmask64(_int(ops[k])):bad()
  elif len(ops)==k+2:
   kind,amt=ops[k+1].split()
   if kind not in("lsl","lsr","asr","ror")or not 0<=_int(amt)<=63:bad()
 elif mn=="mov":
  if ops[1].startswith("#")and not _mov_imm(_int(ops[1])):bad()
 elif mn in("movk","movn","movz"):
  sh=_int(ops[2].split()[1])if len(ops)>2 else 0
  if not(0<=_int(ops[1])<=0xFFFF and sh in(0,16,32,48)):bad()
 elif mn in("lsl","lsr"):
  if ops[2].startswith("#")and not 0<=_int(ops[2])<=63:bad()
 elif mn in("ubfx","bfi","bfxil"):
  lsb,w=_int(ops[2]),_int(ops[3])
  if not(0<=lsb<=63 and 1<=w<=64-lsb):bad()
 elif mn=="extr":
  if not 0<=_int(ops[3])<=63:bad()
 elif mn in("ccmp","ccmn"):
  if(ops[1].startswith("#")and not 0<=_int(ops[1])<=31)or not 0<=_int(ops[2])<=15:bad()
_LC={}
_LS8={("x%d"%t):tuple("x%d"%(t+k)for k in range(8))for t in range(0,24,2)}
def a64_load(text):
 if text in _LC:
  p,lb=_LC[text];return list(p),{k:list(v)for k,v in lb.items()}
 prog,labels=[],{}
 main,cold,cur=[],[],None
 for raw in expand(text).split("\n"):
  s=raw.split("//")[0].strip()
  if s.startswith(".pushsection"):
   assert cur is None;cur=cold
  elif s.startswith(".popsection"):
   assert cur is cold;cur=None
  else:
   (main if cur is None else cold).append(raw)
 assert cur is None
 for raw in main+cold:
  line=raw.split("//")[0].rstrip()
  if not line.strip():continue
  if line.strip().endswith(":"):
   labels.setdefault(line.strip().rstrip(":"),[]).append(len(prog));continue
  ins=line.strip()
  if ins.startswith("."):continue
  parts=ins.split(None,1)
  mn,rest=parts[0],parts[1]if len(parts)>1 else""
  if rest not in _OP:_OP[rest]=_a64_ops(rest)
  ops=list(_OP[rest])
  encodable(mn,ops)
  prog.append((mn,ops))
 for i,(mn,ops)in enumerate(prog):
  lim=1<<15 if mn in("tbz","tbnz")else 1<<20 if mn.startswith("b.")or mn in("cbz","cbnz")else 1<<27 if mn in("b","bl")else 0
  if lim:
   t=ops[-1]
   num=t[: -1].isdigit()and t[-1:]in("f","b")
   pos=labels[t[: -1]]if num else labels[t]
   tgt=([q for q in pos if q>i][0]if t[-1]=="f"else[q for q in pos if q<=i][-1])if num else pos[0]
   if not((num or len(pos)==1)and-lim<=4*(tgt-i)<lim):
    raise AssertionError(("branch reach",mn,t))
 _LC[text]=(tuple(prog),{k:tuple(v)for k,v in labels.items()})
 return prog,labels
def a64_run(prog,labels,entry,regs,mem,limit=10**7,out_regs=None,visits=None,dev=None,ls64=None,wfe=None,ev0=0):
 M=M64
 R={("x%d"%i):0 for i in range(31)};R["sp"]=0;R["xzr"]=0
 R.update(regs);R["xzr"]=0
 lr0=R["x30"]
 fC=fZ=fN=fV=0
 excl,ev=[None],[ev0]
 def mon_wr(a):
  if excl[0]is not None and(a>>11)==excl[0]:excl[0]=None;ev[0]=1
 pc,steps,depth=labels[entry][0],0,0
 def rd(o):
  return R["x"+o[1:]]&0xFFFFFFFF if o[0]=="w"else R[o]
 def wr(o,v):
  if o[0]=="w":
   if o!="wzr":R["x"+o[1:]]=v&0xFFFFFFFF
  elif o!="xzr":R[o]=v&M
 def imm(o):
  if o not in _IM:_IM[o]=int(o.lstrip("#"),0)
  return _IM[o]
 def ea(op,post,regoff=False):
  if op not in _EA:inner=op.strip();_EA[op]=inner.endswith("!"),[p.strip()for p in inner.rstrip("!")[1: -1].split(",")]
  pre,parts=_EA[op]
  base=parts[0]
  if len(parts)>1 and parts[1][:1]=="x":
   assert regoff and post is None and not pre
   sh=parts[2].split()if len(parts)>2 else["lsl","#0"]
   assert len(parts)<=3 and sh[0]=="lsl"and imm(sh[1])in(0,3)
   return(rd(base)+(rd(parts[1])<<imm(sh[1])))&M,None,base
  off=imm(parts[1])if len(parts)>1 else 0
  if post is not None:return rd(base),(rd(base)+imm(post))&M,base
  if pre:return(rd(base)+off)&M,(rd(base)+off)&M,base
  return(rd(base)+off)&M,None,base
 def shifted(v,sh):
  kind,k=sh.split()
  assert kind in("lsl","lsr","asr")and 0<=imm(k)<64
  if kind=="asr":return((v-(v>>63<<64))>>imm(k))&M
  return(v<<imm(k))&M if kind=="lsl"else v>>imm(k)
 def cond(c):
  assert c in("eq","ne","hs","cs","lo","cc","lt","ge")
  return fC==0 if c in("lo","cc")else fC==1 if c in("hs","cs")else fZ==1 if c=="eq"else fZ==0 if c=="ne"else(fN!=fV if c=="lt"else fN==fV)
 def store_dev(addr,v):
  if dev[0]<=addr<dev[0]+dev[1]:dev[2](addr-dev[0],v)
 TG={}
 def tgt(t):
  k=t,pc
  if k not in TG:
   c=labels[t[: -1]]if t[: -1].isdigit()else None
   TG[k]=labels[t][0]if c is None else[p for p in c if p<pc][-1]if t[-1]=="b"else[p for p in c if p>=pc][0]
  return TG[k]
 while True:
  mn,o=prog[pc];steps+=1;assert steps<limit
  if visits is not None and pc in visits:visits[pc]+=1
  pc+=1
  if mn=="madd":wr(o[0],rd(o[1])*rd(o[2])+rd(o[3]))
  elif mn in("mul","umulh"):
   a,b=rd(o[1]),rd(o[2])
   wr(o[0],(a*b)&M if mn=="mul"else(a*b)>>64)
  elif mn in("adds","adc","adcs","add")and not(mn=="add"and(o[2].startswith("#")or len(o)==4)):
   a,b=rd(o[1]),rd(o[2])
   if len(o)==4:
    sh=o[3].split();assert mn=="adds"and sh[0]=="lsl"
    b=(b<<imm(sh[1]))&M
   t=a+b+(fC if mn in("adc","adcs")else 0)
   wr(o[0],t)
   if mn in("adds","adcs"):st=(t>>63)&1;fC,fZ,fN,fV=t>>64,int((t&M)==0),st,int((a>>63)==(b>>63)and st!=(a>>63))
  elif mn=="add":
   wr(o[0],rd(o[1])+(imm(o[2])if o[2].startswith("#")else shifted(rd(o[2]),o[3])))
  elif mn in("sub","subs","sbc","sbcs"):
   a,b=rd(o[1]),imm(o[2])if o[2].startswith("#")else rd(o[2])
   if len(o)==4:
    sh=o[3].split();assert mn in("sub","subs")and sh[0]=="lsl"and not o[2].startswith("#")
    b=(b<<imm(sh[1]))&M
   assert len(o)<=4
   t=a+((~b)&M)+(fC if mn in("sbc","sbcs")else 1)
   wr(o[0],t)
   if mn in("subs","sbcs"):st=(t>>63)&1;fC,fZ,fN,fV=t>>64,int((t&M)==0),st,int((a>>63)!=(b>>63)and st!=(a>>63))
  elif mn=="csel":wr(o[0],rd(o[1])if cond(o[3])else rd(o[2]))
  elif mn in("cmp","cmn","ccmp","ccmn"):
   if mn[1]=="c"and not cond(o[3]):
    nzcv=imm(o[2]);fV,fC,fZ,fN=nzcv&1,nzcv>>1&1,nzcv>>2&1,nzcv>>3&1
   else:
    a,b=rd(o[0]),imm(o[1])if o[1].startswith("#")else rd(o[1])
    if len(o)==3:
     assert mn=="cmp"and o[2][:4]=="asr "
     b=((b-(b>>63<<64))>>imm(o[2][4:]))&M
    t=a+b if mn[-1]=="n"else a+((~b)&M)+1
    st=(t>>63)&1;fC,fZ,fN,fV=t>>64,int((t&M)==0),st,int(((a>>63)==(b>>63)if mn[-1]=="n"else(a>>63)!=(b>>63))and st!=(a>>63))
  elif mn=="tst":
   t=rd(o[0])&(imm(o[1])if o[1].startswith("#")else rd(o[1]))
   fC,fZ=0,int(t==0)
  elif mn in("orr","and","eor","bic","eon"):
   assert mn not in("bic","eon")or not o[2].startswith("#")
   b=imm(o[2])if o[2].startswith("#")else rd(o[2])
   if len(o)==4:assert not o[2].startswith("#");b=shifted(b,o[3])
   assert len(o)<=4
   if mn in("bic","eon"):b= ~b&M
   wr(o[0],rd(o[1])|b if mn=="orr"else rd(o[1])^b if mn in("eor","eon")else rd(o[1])&b)
  elif mn=="lsl":wr(o[0],rd(o[1])<<imm(o[2]))
  elif mn=="cinc":wr(o[0],rd(o[1])+int(cond(o[2])))
  elif mn=="cset":wr(o[0],int(cond(o[1])))
  elif mn=="bfi":
   mask=((1<<imm(o[3]))-1)<<imm(o[2])
   wr(o[0],(rd(o[0])& ~mask)|((rd(o[1])<<imm(o[2]))&mask))
  elif mn=="bfxil":
   mask=(1<<imm(o[3]))-1
   wr(o[0],(rd(o[0])& ~mask)|((rd(o[1])>>imm(o[2]))&mask))
  elif mn=="extr":wr(o[0],(rd(o[1])<<64|rd(o[2]))>>imm(o[3]))
  elif mn=="lsr":wr(o[0],rd(o[1])>>imm(o[2]))
  elif mn=="ubfx":
   assert imm(o[2])+imm(o[3])<=64
   wr(o[0],(rd(o[1])>>imm(o[2]))&((1<<imm(o[3]))-1))
  elif mn=="mvn":assert len(o)==2;wr(o[0],~rd(o[1])&M)
  elif mn=="mov":wr(o[0],imm(o[1])if o[1].startswith("#")else rd(o[1]))
  elif mn=="movn":wr(o[0],~(imm(o[1])<<(imm(o[2].split()[1])if len(o)>2 else 0))&M)
  elif mn=="movk":
   sh=imm(o[2].split()[1])if len(o)>2 else 0
   wr(o[0],(rd(o[0])& ~(0xFFFF<<sh))|((imm(o[1])&0xFFFF)<<sh))
  elif mn in("ldp","stp"):
   addr,newb,base=ea(o[2],o[3]if len(o)>3 else None);va56(addr,16)
   if mn=="ldp":wr(o[0],mem[addr]);wr(o[1],mem[addr+8])
   else:
    mem[addr],mem[addr+8]=rd(o[0]),rd(o[1])
    if dev is not None:store_dev(addr,mem[addr]);store_dev(addr+8,mem[addr+8])
   if newb is not None:R[base]=newb
  elif mn in("ldr","str"):
   addr,newb,base=ea(o[1],o[2]if len(o)>2 else None,regoff=True);va56(addr)
   if mn=="ldr":wr(o[0],mem[addr])
   else:
    mem[addr]=rd(o[0])
    if dev is not None:store_dev(addr,mem[addr])
   if newb is not None:R[base]=newb
  elif mn in("stlr","stlur"):
   inner=o[1].strip()[1: -1].split(",")
   assert len(o)==2 and(len(inner)==1 if mn=="stlr"else-256<=imm(inner[1].strip())<256)
   addr,newb,base=ea(o[1],None);va56(addr);mem[addr]=rd(o[0]);mon_wr(addr)
   if dev is not None:store_dev(addr,mem[addr])
  elif mn=="dmb":
   assert o==["ish"]
   if dev is not None and len(dev)>3:dev[3]()
  elif mn in("ld64b","st64b"):
   t,inner=o[0],o[1].strip()
   assert t in _LS8,"LS64 Xt must be even and <= 22"
   xn=inner[1: -1].replace(",#0","").strip()
   assert len(o)==2 and inner[0]=="["and inner[-1]=="]"and xn[0]=="x"
   assert mn!="ld64b"or xn not in _LS8[t],"LS64 LD64B base overlaps Xt..Xt+7"
   addr=rd(xn)
   assert addr%64==0,"LS64 address not 64-byte aligned"
   va56(addr,64)
   assert ls64 is not None and any(lo<=addr and addr+64<=hi for lo,hi in ls64[mn=="st64b"]),"LS64 outside a permitted region"
   if mn=="ld64b":
    for k,rn in enumerate(_LS8[t]):R[rn]=mem[addr+8*k]&M
   else:
    for k,rn in enumerate(_LS8[t]):mem[addr+8*k]=R[rn]
  elif mn=="dc":
   assert o[0]=="zva"and rd(o[1])%2048==0 and dev is None,"DC ZVA on a 2048-B aligned Normal-memory block only"
   va56(rd(o[1]),2048)
   for k in range(0,2048,8):mem[rd(o[1])+k]=0
  elif mn in("stadd","staddl"):
   addr,newb,base=ea(o[1],None);va56(addr);mem[addr]=(mem[addr]+rd(o[0]))&M
  elif mn in("ldadd","swp"):
   addr,newb,base=ea(o[2],None);va56(addr);old=mem[addr]
   mem[addr]=(old+rd(o[0]))&M if mn=="ldadd"else rd(o[0]);wr(o[1],old)
  elif mn in("b","b.eq","b.ne","b.lo","b.hs","b.cs","b.cc","b.lt","b.ge"):
   if mn=="b"or cond(mn[2:]):pc=tgt(o[0])
  elif mn in("cbz","cbnz","tbz","tbnz"):
   v=rd(o[0])if mn[0]=="c"else rd(o[0])>>imm(o[1])&1
   if(v==0)==(mn in("cbz","tbz")):pc=tgt(o[-1])
  elif mn=="ldaxr":
   addr,_,_=ea(o[1],None);va56(addr);excl[0]=addr>>11;wr(o[0],mem[addr])
  elif mn=="wfe":
   if ev[0]:ev[0]=0
   else:assert wfe is not None,"WFE sleep without wake source";wfe(mem,excl,ev);assert ev[0],"WFE wake left event unset";ev[0]=0
  elif mn=="blr":
   t_a=rd(o[0]);R["x30"]=pc;depth+=1;pc=t_a
  elif mn=="bl":R["x30"]=pc;depth+=1;pc=tgt(o[0])
  elif mn=="ret":
   if depth:depth-=1;pc=R["x30"];continue
   assert R["x30"]==lr0,"outermost RET with a clobbered x30"
   assert dev is None or len(dev)<4 or dev[4](),"a completed message without a DMB ISH before RET"
   if out_regs is not None:out_regs.update(R)
   return steps
  else:raise ValueError("unsupported instruction: "+mn)
