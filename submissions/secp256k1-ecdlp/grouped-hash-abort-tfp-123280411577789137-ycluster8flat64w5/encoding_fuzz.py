import os,random,re,shutil,subprocess,sys,tempfile;HERE=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,os.path.join(HERE,"kernels"));from a64 import encodable,_a64_ops,expand,_SHAPES,_kind,a64_load;CLANG=["clang","-target","aarch64-linux-gnu","-march=armv8.7-a+ls64+lse","-c","-ferror-limit=0"]
if subprocess.run(["which","clang"],capture_output=True).returncode:CLANG=["xcrun"]+CLANG
def errors(aa,a6):
 a4=os.path.join(a6,"f.s");open(a4,"w").write("\t.text\n"+"".join("\t%s\n"%x for x in aa));p=subprocess.run(CLANG+[a4,"-o",a4+".o"],capture_output=True,text=True);a5={int(l.split(":")[1])-2 for l in p.stderr.split("\n")if": error:"in l}
 if a5:ok=[i for i in range(len(aa))if i not in a5];a5|={ok[i]for i in errors([aa[i]for i in ok],a6)}
 return a5
MNEMS=sorted(_SHAPES)+["dmb","dc","mulx","ldur","smulh"]
def gen_shape(rg):
 pool=[lambda:rg.choice(["x1","x9","x30","sp","xzr"]),lambda:rg.choice(["w1","w9","wzr"]),lambda:"#%d"%rg.choice([0,1,7,31,32,63,64,255,4095,4096,-8,65535]),lambda:"#0x%x"%rg.getrandbits(rg.choice([8,32,64])),lambda:rg.choice(["lsl #3","lsr #19","asr #1","ror #7","lsl #40","lsl #12"]),lambda:rg.choice(["eq","hs","lo","ne","zz"]),lambda:rg.choice(["[x2]","[x2,#16]","[sp,#-16]!","[x3,x4]","[x3,#8]"]),lambda:"lbl",lambda:rg.choice(["ish","zva","sy"])];mn=rg.choice(MNEMS);J=[rg.choice(pool)()for _ in range(rg.choice([0,1,2,2,3,3,3,4,4]))];return(mn+" "+",".join(J)).strip()
def gen(rg):
 if rg.random()<0.5:return gen_shape(rg)
 R=lambda:rg.choice(["x1","x2","x3","x9","x30","sp","xzr"]);I=lambda:rg.choice([rg.randrange(-1100,1100),rg.randrange(-1100,1100)//8*8,rg.randrange(0,33000)//8*8,rg.randrange(-300,300),32760,32768,504,512,-512,-520,255,256,-256,-257,0]);a,b,c,k=R(),R(),R(),rg.randrange(25);a9=lambda:rg.choice([rg.getrandbits(64),rg.getrandbits(16)<<rg.choice([0,16,32,48]),((1<<rg.randrange(1,64))-1)<<rg.randrange(0,8),rg.randrange(0,5000),rg.randrange(0,4096)<<12]);return[lambda:"ldp %s,%s,[%s,#%d]"%(a,b,c,I()),lambda:"stp %s,%s,[%s,#%d]!"%(a,b,c,I()),lambda:"ldp %s,%s,[%s],#%d"%(a,b,c,I()),lambda:"ldr %s,[%s,#%d]"%(a,c,I()),lambda:"str %s,[%s],#%d"%(a,c,I()),lambda:"ldr %s,[%s,#%d]!"%(a,c,I()),lambda:"ldr %s,[%s,%s,lsl #%d]"%(a,b,c,rg.choice([0,3,2])),lambda:"%s %s,%s,#%d"%(rg.choice(["add","sub","adds","subs"]),a,b,rg.choice([rg.randrange(-5000,5000),a9()])),lambda:"%s %s,%s,%s,%s #%d"%(rg.choice(["add","adds","sub","subs"]),a,b,c,rg.choice(["lsl","lsr"]),rg.randrange(0,66)),lambda:"%s %s,%s,#0x%x"%(rg.choice(["and","orr","eor"]),a,b,a9()),lambda:"tst %s,#0x%x"%(a,a9()),lambda:"mov %s,#0x%x"%(a,a9()),lambda:"mov %s,%s"%(a,b),lambda:"movk %s,#%d,lsl #%d"%(a,rg.randrange(0,70000),rg.choice([0,8,16,32,48])),lambda:"ubfx %s,%s,#%d,#%d"%(a,b,rg.randrange(0,66),rg.randrange(0,66)),lambda:"bfi %s,%s,#%d,#%d"%(a,b,rg.randrange(0,66),rg.randrange(0,66)),lambda:"bfxil %s,%s,#%d,#%d"%(a,b,rg.randrange(0,66),rg.randrange(0,66)),lambda:"cmp %s,#%d"%(a,rg.randrange(-5000,5000)),lambda:"cmp %s,%s"%(a,b),lambda:"stlur %s,[%s,#%d]"%(a,c,I()),lambda:"ldadd %s,%s,[%s]"%(a,b,c),lambda:"stadd %s,[%s]"%(a,c),lambda:"%s %s,%s,%s"%(rg.choice(["mul","umulh","adcs","sbc"]),a,b,c),lambda:"extr %s,%s,%s,#%d"%(a,b,c,rg.randrange(0,66)),lambda:"lsr %s,%s,#%d"%(a,b,rg.randrange(0,66))][k]()
def main():
 seed,n=(int(sys.argv[1])if len(sys.argv)>1 else 1),(int(sys.argv[2])if len(sys.argv)>2 else 4000);rg,L=random.Random(seed),[]
 while len(L)<n:L.append(gen(rg))
 a7=[]
 for ins in L:
  mn,_,a8=ins.partition(" ")
  try:encodable(mn,_a64_ops(a8));a7.append(True)
  except AssertionError:a7.append(False)
 with tempfile.TemporaryDirectory()as a6:
  a5=errors(L,a6);llvm=[i not in a5 for i in range(len(L))];unsound=[L[i]for i in range(len(L))if a7[i]and not llvm[i]];stricter=[L[i]for i in range(len(L))if llvm[i]and not a7[i]]
  def declared(x):
   mn,_,a8=x.partition(" ");J=_a64_ops(a8)
   if mn=="ldp"and len(J)>1 and J[0]==J[1]:return True
   if mn in("b","bl","cbz","cbnz","tbz","tbnz","dmb","wfe")or mn.startswith("b."):return True
   return any(o.strip()[:1]=="w"for o in J)or mn not in _SHAPES and mn not in("dmb","dc")
  other=[x for x in stricter if not declared(x)];print("seed %d: %d instructions, LLVM accepts %d, a64 accepts %d; a64 accepts but LLVM rejects: %d; a64 stricter: %d ""(declared scope: %d, other: %d)"%(seed,n,sum(llvm),sum(a7),len(unsound),len(stricter),len(stricter)-len(other),len(other)))
  for d in(unsound+other)[:40]:print("  ",d)
  dis=unsound;ab=[]
  for f in("leaf.s","route.s","secp_fe.s"):
   af=[l for l in expand(open(os.path.join(HERE,"kernels",f)).read()).split("\n")];a4=os.path.join(a6,f);open(a4,"w").write("\t.text\n"+"\n".join(af)+"\n");p=subprocess.run(CLANG+[a4,"-o",a4+".o"],capture_output=True,text=True);ac=sum(": error:"in l for l in p.stderr.split("\n"));dump=subprocess.run(([]if shutil.which("llvm-objdump")else["xcrun"])+["llvm-objdump","-d",a4+".o"],capture_output=True,text=True).stdout;ad=sum(1 for l in dump.split("\n")if re.match(r"\s+[0-9a-f]+:\s+(([0-9a-f]{2} ){4}|[0-9a-f]{8}\s)",l));ae=len(a64_load(open(os.path.join(HERE,"kernels",f)).read())[0]);print("%s: LLVM errors %d; LLVM instructions %d, a64_load instructions %d"%(f,ac,ad,ae))
   if ac or ad!=ae:ab.append(f)
 return len(dis)+len(ab)
if __name__=="__main__":sys.exit(1 if main()else 0)
