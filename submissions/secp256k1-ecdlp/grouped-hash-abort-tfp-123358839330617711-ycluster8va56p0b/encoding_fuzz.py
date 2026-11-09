import os,random,re,shutil,subprocess,sys,tempfile
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(HERE,"kernels"))
from a64 import encodable,_a64_ops,expand,_SHAPES,_kind,a64_load
CLANG=["clang","-target","aarch64-linux-gnu","-march=armv8.7-a+ls64+lse","-c","-ferror-limit=0"]
if subprocess.run(["which","clang"],capture_output=True).returncode:CLANG=["xcrun"]+CLANG
def errors(lines,tmp):
 src=os.path.join(tmp,"f.s")
 open(src,"w").write("\t.text\n"+"".join("\t%s\n"%x for x in lines))
 p=subprocess.run(CLANG+[src,"-o",src+".o"],capture_output=True,text=True)
 err={int(l.split(":")[1])-2 for l in p.stderr.split("\n")if": error:"in l}
 if err:
  ok=[i for i in range(len(lines))if i not in err]
  err|={ok[i]for i in errors([lines[i]for i in ok],tmp)}
 return err
MNEMS=sorted(_SHAPES)+["dmb","dc","mulx","ldur","smulh"]
def gen_shape(rg):
 pool=[lambda:rg.choice(["x1","x9","x30","sp","xzr"]),lambda:rg.choice(["w1","w9","wzr"]),lambda:"#%d"%rg.choice([0,1,7,31,32,63,64,255,4095,4096,-8,65535]),lambda:"#0x%x"%rg.getrandbits(rg.choice([8,32,64])),lambda:rg.choice(["lsl #3","lsr #19","asr #1","ror #7","lsl #40","lsl #12"]),lambda:rg.choice(["eq","hs","lo","ne","zz"]),lambda:rg.choice(["[x2]","[x2,#16]","[sp,#-16]!","[x3,x4]","[x3,#8]"]),lambda:"lbl",lambda:rg.choice(["ish","zva","sy"])]
 mn=rg.choice(MNEMS)
 ops=[rg.choice(pool)()for _ in range(rg.choice([0,1,2,2,3,3,3,4,4]))]
 return(mn+" "+",".join(ops)).strip()
def gen(rg):
 if rg.random()<0.5:return gen_shape(rg)
 R=lambda:rg.choice(["x1","x2","x3","x9","x30","sp","xzr"])
 I=lambda:rg.choice([rg.randrange(-1100,1100),rg.randrange(-1100,1100)//8*8,rg.randrange(0,33000)//8*8,rg.randrange(-300,300),32760,32768,504,512,-512,-520,255,256,-256,-257,0])
 a,b,c,k=R(),R(),R(),rg.randrange(25)
 big=lambda:rg.choice([rg.getrandbits(64),rg.getrandbits(16)<<rg.choice([0,16,32,48]),((1<<rg.randrange(1,64))-1)<<rg.randrange(0,8),rg.randrange(0,5000),rg.randrange(0,4096)<<12])
 return[lambda:"ldp %s,%s,[%s,#%d]"%(a,b,c,I()),lambda:"stp %s,%s,[%s,#%d]!"%(a,b,c,I()),lambda:"ldp %s,%s,[%s],#%d"%(a,b,c,I()),lambda:"ldr %s,[%s,#%d]"%(a,c,I()),lambda:"str %s,[%s],#%d"%(a,c,I()),lambda:"ldr %s,[%s,#%d]!"%(a,c,I()),lambda:"ldr %s,[%s,%s,lsl #%d]"%(a,b,c,rg.choice([0,3,2])),lambda:"%s %s,%s,#%d"%(rg.choice(["add","sub","adds","subs"]),a,b,rg.choice([rg.randrange(-5000,5000),big()])),lambda:"%s %s,%s,%s,%s #%d"%(rg.choice(["add","adds","sub","subs"]),a,b,c,rg.choice(["lsl","lsr"]),rg.randrange(0,66)),lambda:"%s %s,%s,#0x%x"%(rg.choice(["and","orr","eor"]),a,b,big()),lambda:"tst %s,#0x%x"%(a,big()),lambda:"mov %s,#0x%x"%(a,big()),lambda:"mov %s,%s"%(a,b),lambda:"movk %s,#%d,lsl #%d"%(a,rg.randrange(0,70000),rg.choice([0,8,16,32,48])),lambda:"ubfx %s,%s,#%d,#%d"%(a,b,rg.randrange(0,66),rg.randrange(0,66)),lambda:"bfi %s,%s,#%d,#%d"%(a,b,rg.randrange(0,66),rg.randrange(0,66)),lambda:"bfxil %s,%s,#%d,#%d"%(a,b,rg.randrange(0,66),rg.randrange(0,66)),lambda:"cmp %s,#%d"%(a,rg.randrange(-5000,5000)),lambda:"cmp %s,%s"%(a,b),lambda:"stlur %s,[%s,#%d]"%(a,c,I()),lambda:"ldadd %s,%s,[%s]"%(a,b,c),lambda:"stadd %s,[%s]"%(a,c),lambda:"%s %s,%s,%s"%(rg.choice(["mul","umulh","adcs","sbc"]),a,b,c),lambda:"extr %s,%s,%s,#%d"%(a,b,c,rg.randrange(0,66)),lambda:"lsr %s,%s,#%d"%(a,b,rg.randrange(0,66))][k]()
def main():
 seed,n=(int(sys.argv[1])if len(sys.argv)>1 else 1),(int(sys.argv[2])if len(sys.argv)>2 else 4000)
 rg,L=random.Random(seed),[]
 while len(L)<n:L.append(gen(rg))
 mine=[]
 for ins in L:
  mn,_,rest=ins.partition(" ")
  try:encodable(mn,_a64_ops(rest));mine.append(True)
  except AssertionError:mine.append(False)
 with tempfile.TemporaryDirectory()as tmp:
  err=errors(L,tmp)
  llvm=[i not in err for i in range(len(L))]
  unsound=[L[i]for i in range(len(L))if mine[i]and not llvm[i]]
  stricter=[L[i]for i in range(len(L))if llvm[i]and not mine[i]]
  def declared(x):
   mn,_,rest=x.partition(" ");ops=_a64_ops(rest)
   if mn=="ldp"and len(ops)>1 and ops[0]==ops[1]:return True
   if mn in("b","bl","cbz","cbnz","tbz","tbnz","dmb","wfe")or mn.startswith("b."):return True
   return any(o.strip()[:1]=="w"for o in ops)or mn not in _SHAPES and mn not in("dmb","dc")
  other=[x for x in stricter if not declared(x)]
  print("seed %d: %d instructions, LLVM accepts %d, a64 accepts %d; a64 accepts but LLVM rejects: %d; a64 stricter: %d ""(declared scope: %d, other: %d)"%(seed,n,sum(llvm),sum(mine),len(unsound),len(stricter),len(stricter)-len(other),len(other)))
  for d in(unsound+other)[:40]:print("  ",d)
  dis=unsound
  bad_k=[]
  for f in("leaf.s","route.s","secp_fe.s"):
   body=[l for l in expand(open(os.path.join(HERE,"kernels",f)).read()).split("\n")]
   src=os.path.join(tmp,f)
   open(src,"w").write("\t.text\n"+"\n".join(body)+"\n")
   p=subprocess.run(CLANG+[src,"-o",src+".o"],capture_output=True,text=True)
   nerr=sum(": error:"in l for l in p.stderr.split("\n"))
   dump=subprocess.run(([]if shutil.which("llvm-objdump")else["xcrun"])+["llvm-objdump","-d",src+".o"],capture_output=True,text=True).stdout
   ndis=sum(1 for l in dump.split("\n")if re.match(r"\s+[0-9a-f]+:\s+(([0-9a-f]{2} ){4}|[0-9a-f]{8}\s)",l))
   nprog=len(a64_load(open(os.path.join(HERE,"kernels",f)).read())[0])
   print("%s: LLVM errors %d; LLVM instructions %d, a64_load instructions %d"%(f,nerr,ndis,nprog))
   if nerr or ndis!=nprog:bad_k.append(f)
 return len(dis)+len(bad_k)
if __name__=="__main__":
 sys.exit(1 if main()else 0)
