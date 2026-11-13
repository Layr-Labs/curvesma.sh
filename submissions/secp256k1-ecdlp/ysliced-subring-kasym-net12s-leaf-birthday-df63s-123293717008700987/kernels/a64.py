import re as _re;M64=2**64-1
DEV=[]
def nd(o):assert not any(a<=o<b for a,b in DEV),"Device atomic 0x%x"%o;return o
def va55(o,n=8):assert 0<=o and o+n<=1<<64 and(o+n<=1<<55 or o>=2**64-(1<<55)),"data address 0x%x outside the two 55-bit VA ranges"%o;return o
def expand(text):
 o,out,a6={},[],None
 for J in text.split("\n"):
  s=J.split("//")[0].strip()
  if s.startswith(".macro"):a4=s[len(".macro"):].replace(","," ").split();a6=o[a4[0]]=(a4[1:],[]);continue
  if s.startswith(".endm"):a6=None;continue
  (a6[1]if a6 is not None else out).append(J)
 a7=[0]
 def emit(aa,a5):
  assert a5<16;a8=[]
  for J in aa:
   s=J.split("//")[0].strip();I=s.split()[0]if s else""
   if not s.endswith(":")and I in o:
    a9,ab=o[I];ac=[a.strip()for a in s[len(I):].split(",")]if s[len(I):].strip()else[];a7[0]+=1;sub=[]
    for b in ab:
     for p_,a in zip(a9,ac):b=b.replace("\\"+p_,a)
     sub+=[b.replace("\\@",str(a7[0]))]
    a8+=emit(sub,a5+1)
   else:
    a8+=[J]
  return a8
 return"\n".join(emit(out,0))
def _ops(s):
 J,I,o=[],"",0
 for ch in s:
  if ch=="[":o+=1
  if ch=="]":o-=1
  if ch==","and o==0:J+=[I.strip()];I=""
  else:
   I+=ch
 if I.strip():J+=[I.strip()]
 return J
def _i(o):
 m=_re.fullmatch(r"#?(-?(0x[0-9a-fA-F]+|0b[01]+|[0-9]+))",o.strip())
 if not m:raise ValueError("bad immediate: "+o)
 return int(m.group(1),0)
def _xr(o):return _re.fullmatch(r"x([12]?[0-9]|30)|sp|xzr",o.strip(),_re.ASCII)is not None and not _re.fullmatch(r"x0[0-9]",o.strip())
def _bm(v):
 if not-(1<<63)<=v<1<<64:return False
 v&=M64
 if v in(0,M64):return False
 for e in(2,4,8,16,32,64):
  m=(1<<e)-1;o=v&m
  if any((v>>k)&m!=o for k in range(0,64,e)):continue
  if o in(0,m):continue
  if sum(1 for i in range(e)if o>>i&1 and not o>>((i-1)%e)&1)==1:return True
 return False
def _aim(v):v=abs(v);return 0<=v<4096 or(v%4096==0 and v>>12<4096)
def _mvi(v):
 if not-(1<<63)<=v<1<<64:return False
 v&=M64;return any(v&~(0xFFFF<<s)&M64==0 for s in(0,16,32,48))or any(~v&M64&~(0xFFFF<<s)&M64==0 for s in(0,16,32,48))or _bm(v)
def _mem(I):
 I=I.strip();pre=I.endswith("!");I=I[:-1]if pre else I
 if not(I[:1]=="["and I[-1:]=="]"):raise AssertionError("not a memory operand: "+I)
 o=[q.strip()for q in I[1:-1].split(",")]
 if not _xr(o[0])or o[0]=="xzr":raise AssertionError("bad base register: "+I)
 if pre and len(o)==1:raise AssertionError("not encodable: "+I)
 if len(o)>1 and not o[1][:1]=="#":return o[0],0,pre,o[1:]
 if len(o)>2:raise AssertionError("not encodable: "+I)
 try:return o[0],_i(o[1])if len(o)>1 else 0,pre,None
 except ValueError:raise AssertionError("not encodable: "+I)
def _rgk(J,I):
 regs=[o.strip()for o in I if _xr(o.strip())]
 if J in("add","sub","adds","subs","cmp","cmn"):
  k=1 if J in("cmp","cmn")else 2;a8,a5=([]if k==1 else[I[0]]),I[k-1:k];a6=len(I)>k and I[k][:1]=="#";a7=a6 or(len(I)<=k+2 and"sp"in I[:k]and(len(I)==k+1 or I[k+1].split()[0]=="lsl"))
  if not a7:return"sp"not in I
  if"sp"in I[k:]:return False
  if J in("adds","subs","cmp","cmn"):return"xzr"not in a5 and"sp"not in a8
  return"xzr"not in I[:k]
 if J in("and","orr","eor")and len(I)==3 and I[2][:1]=="#":return I[1]!="sp"and I[0]!="xzr"
 if J=="mov":
  if I[1][:1]=="#":v=_i(I[1])&M64;a4=any(v&~(0xFFFF<<q)&M64==0 or~v&M64&~(0xFFFF<<q)&M64==0 for q in(0,16,32,48));return a4 if I[0]=="xzr"else _bm(v)if I[0]=="sp"else True
  return not("sp"in I and"xzr"in I)
 return"sp"not in regs or J in("ldp","stp","ldr","str","stlur","stlr","ld64b","st64b","stadd","staddl","ldadd","swp","ldaxr")
_COND={"eq","ne","hs","lo","ge","lt"};_D3=["rrr"];_AR=["rri","rrr","rrrs"];_SHAPES=dict(mul=_D3,umulh=_D3,adc=_D3,adcs=_D3,sbc=_D3,sbcs=_D3,madd=["rrrr"],csel=["rrrc"],cset=["rc"],cinc=["rrc"],mvn=["rr","rrs"],lsl=["rri","rrr"],lsr=["rri","rrr"],ubfx=["rrii"],bfi=["rrii"],bfxil=["rrii"],extr=["rrri"],add=_AR,adds=_AR,sub=_AR,subs=_AR,cmp=["ri","rr","rrs"],cmn=["ri","rr","rrs"],tst=["ri","rr","rrs"],mov=["rr","ri"],movk=["ri","ris"],movn=["ri","ris"],ldp=["rrm","rrmi"],stp=["rrm","rrmi"],ldr=["rm","rmi"],str=["rm","rmi"],stlur=["rm"],stlr=["rm"],ld64b=["rm"],st64b=["rm"],stadd=["rm"],staddl=["rm"],ldadd=["rrm"],swp=["rrm"],ldaxr=["rm"],blr=["r"],wfe=[""],ret=["","r"],b=["l"],bl=["l"],cbz=["rl"],cbnz=["rl"],cbhs=["rrl"],cbhi=["rrl"],tbz=["ril"],tbnz=["ril"],**{k:_AR for k in("and","orr","eor","bic")},**{"b."+c:["l"]for c in _COND})
def _kd(o):
 o=o.strip()
 if _xr(o):return"r"
 if _re.fullmatch(r"w([12]?[0-9]|30)|wzr",o):return"w"
 if o[:1]=="#":
  try:_i(o);return"i"
  except ValueError:return"?"
 if o.startswith("["):return"m"
 if _re.fullmatch(r"(lsl|lsr|asr|ror) #[0-9]+",o):return"s"
 if o in _COND:return"c"
 if _re.fullmatch(r"[xwqdshbp][0-9]+|[vz][0-9]+(\.[0-9]*[bhsdq])?|zt[0-9]+|ffr|lr|fp|wsp|ip[01]|nzcv|al|nv",o):return"?"
 if _re.fullmatch(r"[A-Za-z_.][A-Za-z0-9_.]*|[0-9]+[fb]",o):return"l"
 return"?"
def _shp(mn,J):
 I=[_kd(o)for o in J]
 if"w"in I:
  if mn!="cset"or"r"in I:return False
  I=["r"if k=="w"else k for k in I]
 return"".join(I)in _SHAPES[mn]
_OK,_IM,_EA,_OP=set(),{},{},{}
def encodable(I,o):
 k=(I,*o)
 if k in _OK:return
 try:_enc(I,o)
 except ValueError as ex:raise AssertionError("not encodable (bad immediate): %s %s (%s)"%(I,",".join(o),ex))
 _OK.add(k)
def _enc(J,I):
 def bad(w=""):raise AssertionError("not encodable%s: %s %s"%(w,J,",".join(I)))
 if J=="dmb":
  if I!=["ish"]:bad()
  return
 if J=="dc":
  if not(len(I)==2 and I[0]=="zva"and _xr(I[1])and I[1]not in("sp","xzr")):bad()
  return
 if J not in _SHAPES or not _shp(J,I):bad()
 if not _rgk(J,I):bad()
 if J in("ldr","str","ldp","stp","stlur","stlr","ld64b","st64b","ldaxr")and I[0]=="sp":bad()
 if J in("stadd","staddl","ldadd","swp")and"sp"in I[:-1]:bad()
 if J in("tbz","tbnz")and not 0<=_i(I[1])<=63:bad()
 for o in I:
  if _kd(o)=="s":
   kd,a6=o.split()[0],_i(o.split()[1])
   if not 0<=a6<=63 or(J in("movk","movn")and kd!="lsl")or(J in("add","adds","sub","subs","cmp","cmn")and kd=="ror"):bad()
 if J in("cbz","cbnz","tbz","tbnz")and I[0]=="sp"or J in("cbhs","cbhi")and"sp"in I[:2]:bad()
 if J in("ldp","stp"):
  if not(3<=len(I)<=4 and _xr(I[0])and _xr(I[1])and"sp"not in I[:2]):bad()
  bs,off,pre,a4=_mem(I[2])
  if a4 is not None:bad()
  if len(I)==4:
   if pre or off:bad()
   off=_i(I[3])
  if off%8 or not-512<=off<=504 or((pre or len(I)==4)and bs!="sp"and bs in I[:2])or(J=="ldp"and I[0]==I[1]):bad()
 elif J in("ldr","str"):
  if not(2<=len(I)<=3 and _xr(I[0])and I[0]!="sp"):bad()
  bs,off,pre,a4=_mem(I[1])
  if a4 is not None:
   if pre or len(I)==3 or not _xr(a4[0])or a4[0]=="sp"or len(a4)>2 or(len(a4)==2 and a4[1].replace(" ","")not in("lsl#0","lsl#3")):bad()
   return
  if len(I)==3:
   if pre or off:bad()
   off=_i(I[2])
  if pre or len(I)==3:
   if not-256<=off<=255 or(bs==I[0]and bs!="xzr"):bad()
  elif not((off%8==0 and 0<=off<=32760)or-256<=off<=255):bad()
 elif J=="stlur":
  bs,off,pre,a4=_mem(I[1])
  if pre or a4 is not None or not-256<=off<=255:bad()
 elif J in("stlr","ldaxr","ld64b","st64b","stadd","staddl","ldadd","swp"):
  bs,off,pre,a4=_mem(I[-1])
  if pre or a4 is not None or","in I[-1]:bad()
  if J in("ld64b","st64b")and not(I[0][:1]=="x"and I[0][1:].isdigit()and int(I[0][1:])%2==0 and int(I[0][1:])<=22):bad(" LS64")
 elif J in("add","adds","sub","subs","cmp","cmn"):
  k=1 if J in("cmp","cmn")else 2
  if len(I)==k+1 and I[k][:1]=="#":
   if not _aim(_i(I[k])):bad()
  elif len(I)==k+2:
   if I[k][:1]=="#":bad()
   kd,a6=I[k+1].split();a8=4 if"sp"in I[:k+1]else 63
   if kd not in("lsl","lsr","asr")or not 0<=_i(a6)<=a8 or("sp"in I[:k+1]and kd!="lsl"):bad()
 elif J in("and","orr","eor","tst","bic"):
  k=1 if J=="tst"else 2
  if I[k][:1]=="#":
   if J=="bic"or len(I)!=k+1 or not _bm(_i(I[k])):bad()
  elif len(I)==k+2:
   kd,a6=I[k+1].split()
   if kd not in("lsl","lsr","asr","ror")or not 0<=_i(a6)<=63:bad()
 elif J=="mov":
  if I[1][:1]=="#"and not _mvi(_i(I[1])):bad()
 elif J in("movk","movn"):
  sh=_i(I[2].split()[1])if len(I)>2 else 0
  if not(0<=_i(I[1])<=0xFFFF and sh in(0,16,32,48)):bad()
 elif J in("lsl","lsr"):
  if I[2][:1]=="#"and not 0<=_i(I[2])<=63:bad()
 elif J in("ubfx","bfi","bfxil"):
  a7,w=_i(I[2]),_i(I[3])
  if not(0<=a7<=63 and 1<=w<=64-a7):bad()
 elif J=="extr":
  if not 0<=_i(I[3])<=63:bad()
_LC={};_LS8={("x%d"%t):tuple("x%d"%(t+k)for k in range(8))for t in range(0,24,2)}
def a64_load(text):
 if text in _LC:p,lb=_LC[text];return list(p),{k:list(v)for k,v in lb.items()}
 o,d=[],{};main,cold,I=[],[],None
 for a6 in expand(text).split("\n"):
  s=a6.split("//")[0].strip()
  if s.startswith(".pushsection"):assert I is None;I=cold
  elif s.startswith(".popsection"):assert I is cold;I=None
  else:
   (main if I is None else cold).append(a6)
 assert I is None
 for a6 in main+cold:
  a5=a6.split("//")[0].rstrip()
  if not a5.strip():continue
  if a5.strip().endswith(":"):d.setdefault(a5.strip().rstrip(":"),[]).append(len(o));continue
  ab=a5.strip()
  if ab.startswith("."):continue
  a4=ab.split(None,1);mn,J=a4[0],a4[1]if len(a4)>1 else""
  if J not in _OP:_OP[J]=_ops(J)
  a8=list(_OP[J]);encodable(mn,a8);o+=[(mn,a8)]
 for i,(mn,a8)in enumerate(o):
  aa=1<<15 if mn in("tbz","tbnz")else 1<<20 if mn.startswith("b.")or mn in("cbz","cbnz")else 1<<10 if mn in("cbhs","cbhi")else 1<<27 if mn in("b","bl")else 0
  if aa:
   t=a8[-1];a9=t[:-1].isdigit()and t[-1:]in("f","b");a7=d[t[:-1]]if a9 else d[t];ac=([q for q in a7 if q>i][0]if t[-1]=="f"else[q for q in a7 if q<=i][-1])if a9 else a7[0]
   if not((a9 or len(a7)==1)and-aa<=4*(ac-i)<aa):raise AssertionError(("branch reach",mn,t))
 _LC[text]=(tuple(o),{k:tuple(v)for k,v in d.items()});return o,d
def a64_run(ab,d,entry,regs,mem,limit=10**7,out_regs=None,visits=None,dev=None,ls64=None,wfe=None,ev0=0):
 M=M64;R={("x%d"%i):0 for i in range(31)};R["sp"]=0;R["xzr"]=0;R.update(regs);R["xzr"]=0;ad=R["x30"];fC=fZ=fN=fV=0;a8,ev=[None],[ev0]
 def mon_wr(a):
  if a8[0]is not None and(a>>11)==a8[0]:a8[0]=None;ev[0]=1
 pc,steps,a7=d[entry][0],0,0
 def rd(o):return R[o]
 def wr(o,v):
  if o[0]=="w":o="x"+o[1:]
  if o!="xzr":R[o]=v&M
 def im(o):
  if o not in _IM:_IM[o]=int(o.lstrip("#"),0)
  return _IM[o]
 def ea(op,post,regoff=False):
  if op not in _EA:a5=op.strip();_EA[op]=a5.endswith("!"),[p.strip()for p in a5.rstrip("!")[1:-1].split(",")]
  pre,a4=_EA[op];bs=a4[0]
  if len(a4)>1 and a4[1][:1]=="x":assert regoff and post is None and not pre;sh=a4[2].split()if len(a4)>2 else["lsl","#0"];assert len(a4)<=3 and sh[0]=="lsl"and im(sh[1])in(0,3);return(rd(bs)+(rd(a4[1])<<im(sh[1])))&M,None,bs
  off=im(a4[1])if len(a4)>1 else 0
  if post is not None:return rd(bs),(rd(bs)+im(post))&M,bs
  if pre:return(rd(bs)+off)&M,(rd(bs)+off)&M,bs
  return(rd(bs)+off)&M,None,bs
 def shifted(v,sh):
  kd,k=sh.split();assert kd in("lsl","lsr","asr")and 0<=im(k)<64
  if kd=="asr":return((v-(v>>63<<64))>>im(k))&M
  return(v<<im(k))&M if kd=="lsl"else v>>im(k)
 def cond(c):assert c in _COND;return fC==0 if c=="lo"else fC==1 if c=="hs"else fZ==1 if c=="eq"else fZ==0 if c=="ne"else(fN!=fV if c=="lt"else fN==fV)
 def sd(I,v):
  if dev[0]<=I<dev[0]+dev[1]:dev[2](I-dev[0],v)
 TG={}
 def tgt(t):
  k=t,pc
  if k not in TG:c=d[t[:-1]]if t[:-1].isdigit()else None;TG[k]=d[t][0]if c is None else[p for p in c if p<pc][-1]if t[-1]=="b"else[p for p in c if p>=pc][0]
  return TG[k]
 while True:
  J,o=ab[pc];steps+=1;assert steps<limit
  if visits is not None and pc in visits:visits[pc]+=1
  pc+=1
  if J=="madd":wr(o[0],rd(o[1])*rd(o[2])+rd(o[3]))
  elif J in("mul","umulh"):a,b=rd(o[1]),rd(o[2]);wr(o[0],(a*b)&M if J=="mul"else(a*b)>>64)
  elif J in("adds","adc","adcs","add")and not(J=="add"and(o[2][:1]=="#"or len(o)==4)):
   a,b=rd(o[1]),rd(o[2])
   if len(o)==4:sh=o[3].split();assert J=="adds"and sh[0]=="lsl";b=(b<<im(sh[1]))&M
   t=a+b+(fC if J in("adc","adcs")else 0);wr(o[0],t)
   if J in("adds","adcs"):st=(t>>63)&1;fC,fZ,fN,fV=t>>64,int((t&M)==0),st,int((a>>63)==(b>>63)and st!=(a>>63))
  elif J=="add":wr(o[0],rd(o[1])+(im(o[2])if o[2][:1]=="#"else shifted(rd(o[2]),o[3])))
  elif J in("sub","subs","sbc","sbcs"):
   a,b=rd(o[1]),im(o[2])if o[2][:1]=="#"else rd(o[2])
   if len(o)==4:sh=o[3].split();assert J in("sub","subs")and sh[0]=="lsl"and not o[2][:1]=="#";b=(b<<im(sh[1]))&M
   assert len(o)<=4;t=a+((~b)&M)+(fC if J in("sbc","sbcs")else 1);wr(o[0],t)
   if J in("subs","sbcs"):st=(t>>63)&1;fC,fZ,fN,fV=t>>64,int((t&M)==0),st,int((a>>63)!=(b>>63)and st!=(a>>63))
  elif J=="csel":wr(o[0],rd(o[1])if cond(o[3])else rd(o[2]))
  elif J in("cmp","cmn"):
   a,b=rd(o[0]),im(o[1])if o[1][:1]=="#"else rd(o[1])
   if len(o)==3:assert J=="cmp"and o[2][:4]=="asr ";b=((b-(b>>63<<64))>>im(o[2][4:]))&M
   t=a+b if J[-1]=="n"else a+((~b)&M)+1;st=(t>>63)&1;fC,fZ,fN,fV=t>>64,int((t&M)==0),st,int(((a>>63)==(b>>63)if J[-1]=="n"else(a>>63)!=(b>>63))and st!=(a>>63))
  elif J=="tst":t=rd(o[0])&(im(o[1])if o[1][:1]=="#"else rd(o[1]));fC,fZ=0,int(t==0)
  elif J in("orr","and","eor","bic"):
   assert J!="bic"or not o[2][:1]=="#";b=im(o[2])if o[2][:1]=="#"else rd(o[2])
   if len(o)==4:assert not o[2][:1]=="#";b=shifted(b,o[3])
   assert len(o)<=4
   if J=="bic":b=~b&M
   wr(o[0],rd(o[1])|b if J=="orr"else rd(o[1])^b if J=="eor"else rd(o[1])&b)
  elif J=="lsl":wr(o[0],rd(o[1])<<im(o[2]))
  elif J=="cinc":wr(o[0],rd(o[1])+int(cond(o[2])))
  elif J=="cset":wr(o[0],int(cond(o[1])))
  elif J=="bfi":mask=((1<<im(o[3]))-1)<<im(o[2]);wr(o[0],(rd(o[0])&~mask)|((rd(o[1])<<im(o[2]))&mask))
  elif J=="bfxil":mask=(1<<im(o[3]))-1;wr(o[0],(rd(o[0])&~mask)|((rd(o[1])>>im(o[2]))&mask))
  elif J=="extr":wr(o[0],(rd(o[1])<<64|rd(o[2]))>>im(o[3]))
  elif J=="lsr":wr(o[0],rd(o[1])>>im(o[2]))
  elif J=="ubfx":assert im(o[2])+im(o[3])<=64;wr(o[0],(rd(o[1])>>im(o[2]))&((1<<im(o[3]))-1))
  elif J=="mvn":assert len(o)==2;wr(o[0],~rd(o[1])&M)
  elif J=="mov":wr(o[0],im(o[1])if o[1][:1]=="#"else rd(o[1]))
  elif J=="movn":wr(o[0],~(im(o[1])<<(im(o[2].split()[1])if len(o)>2 else 0))&M)
  elif J=="movk":sh=im(o[2].split()[1])if len(o)>2 else 0;wr(o[0],(rd(o[0])&~(0xFFFF<<sh))|((im(o[1])&0xFFFF)<<sh))
  elif J in("ldp","stp"):
   I,a6,bs=ea(o[2],o[3]if len(o)>3 else None);va55(I,16)
   if J=="ldp":wr(o[0],mem[I]);wr(o[1],mem[I+8])
   else:
    mem[I],mem[I+8]=rd(o[0]),rd(o[1])
    if dev is not None:sd(I,mem[I]);sd(I+8,mem[I+8])
   if a6 is not None:R[bs]=a6
  elif J in("ldr","str"):
   I,a6,bs=ea(o[1],o[2]if len(o)>2 else None,regoff=True);va55(I)
   if J=="ldr":
    if wfe is not None and not mem[I]:wfe(mem,a8,ev)
    wr(o[0],mem[I])
   else:
    mem[I]=rd(o[0])
    if dev is not None:sd(I,mem[I])
   if a6 is not None:R[bs]=a6
  elif J in("stlr","stlur"):
   a5=o[1].strip()[1:-1].split(",");assert len(o)==2 and(len(a5)==1 if J=="stlr"else-256<=im(a5[1].strip())<256);I,a6,bs=ea(o[1],None);va55(I);mem[I]=rd(o[0]);mon_wr(I)
   if dev is not None:sd(I,mem[I])
  elif J=="dmb":
   assert o==["ish"]
   if dev is not None and len(dev)>3:dev[3]()
  elif J in("ld64b","st64b"):
   t,a5=o[0],o[1].strip();assert t in _LS8,"LS64 Xt must be even and <= 22";xn=a5[1:-1].replace(",#0","").strip();assert len(o)==2 and a5[0]=="["and a5[-1]=="]"and xn[0]=="x";assert J!="ld64b"or xn not in _LS8[t],"LS64 LD64B base overlaps Xt..Xt+7";I=rd(xn);assert I%64==0,"LS64 address not 64-byte aligned";va55(I,64);assert ls64 is not None and any(lo<=I and I+64<=hi for lo,hi in ls64[J=="st64b"]),"LS64 outside a permitted region"
   if J=="ld64b":
    for k,rn in enumerate(_LS8[t]):R[rn]=mem[I+8*k]&M
   else:
    for k,rn in enumerate(_LS8[t]):mem[I+8*k]=R[rn]
  elif J=="dc":
   assert o[0]=="zva"and rd(o[1])%2048==0 and dev is None,"DC ZVA on a 2048-B aligned Normal-memory block only";va55(nd(rd(o[1])),2048)
   for k in range(0,2048,8):mem[rd(o[1])+k]=0
  elif J in("stadd","staddl"):I,a6,bs=ea(o[1],None);va55(nd(I));mem[I]=(mem[I]+rd(o[0]))&M
  elif J in("ldadd","swp"):I,a6,bs=ea(o[2],None);va55(nd(I));aa=mem[I];mem[I]=(aa+rd(o[0]))&M if J=="ldadd"else rd(o[0]);wr(o[1],aa)
  elif J=="b"or J[:2]=="b.":
   if J=="b"or cond(J[2:]):pc=tgt(o[0])
  elif J in("cbhs","cbhi"):
   if rd(o[0])>rd(o[1])or J=="cbhs"and rd(o[0])==rd(o[1]):pc=tgt(o[2])
  elif J in("cbz","cbnz","tbz","tbnz"):
   v=rd(o[0])if J[0]=="c"else rd(o[0])>>im(o[1])&1
   if(v==0)==(J in("cbz","tbz")):pc=tgt(o[-1])
  elif J=="ldaxr":I,_,_=ea(o[1],None);va55(nd(I));a8[0]=I>>11;wr(o[0],mem[I])
  elif J=="wfe":
   if ev[0]:ev[0]=0
   else:assert wfe is not None,"WFE sleep without wake source";wfe(mem,a8,ev);assert ev[0],"WFE wake left event unset";ev[0]=0
  elif J=="blr":ac=rd(o[0]);R["x30"]=pc;a7+=1;pc=ac
  elif J=="bl":R["x30"]=pc;a7+=1;pc=tgt(o[0])
  elif J=="ret":
   if a7:a7-=1;pc=R["x30"];continue
   assert R["x30"]==ad,"outermost RET with a clobbered x30";assert dev is None or len(dev)<4 or dev[4](),"a completed message without a DMB ISH before RET"
   if out_regs is not None:out_regs.update(R)
   return steps
  else:raise ValueError("unsupported instruction: "+J)
