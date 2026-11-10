import os,random,re,shutil,subprocess,sys,tempfile;HERE=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,os.path.join(HERE,"kernels"));from a64 import encodable,_a64_ops,expand,_SHAPES,_kind,a64_load;CLANG=["clang","-target","aarch64-linux-gnu","-march=armv8.7-a+ls64+lse","-c","-ferror-limit=0"]
if subprocess.run(["which","clang"],capture_output=True).returncode:CLANG=["xcrun"]+CLANG
def errors(o,a4):
 J=os.path.join(a4,"f.s");open(J,"w").write("\t.text\n"+"".join("\t%s\n"%x for x in o));p=subprocess.run(CLANG+[J,"-o",J+".o"],capture_output=True,text=True);I={int(l.split(":")[1])-2 for l in p.stderr.split("\n")if": error:"in l}
 if I:ok=[i for i in range(len(o))if i not in I];I|={ok[i]for i in errors([o[i]for i in ok],a4)}
 return I
MNEMS=sorted(_SHAPES)+["dmb","dc","mulx","ldur","smulh"]
def gen_shape(o):pool=[lambda:o.choice(["x1","x9","x30","sp","xzr"]),lambda:o.choice(["w1","w9","wzr"]),lambda:"#%d"%o.choice([0,1,7,31,32,63,64,255,4095,4096,-8,65535]),lambda:"#0x%x"%o.getrandbits(o.choice([8,32,64])),lambda:o.choice(["lsl #3","lsr #19","asr #1","ror #7","lsl #40","lsl #12"]),lambda:o.choice(["eq","hs","lo","ne","zz"]),lambda:o.choice(["[x2]","[x2,#16]","[sp,#-16]!","[x3,x4]","[x3,#8]"]),lambda:"lbl",lambda:o.choice(["ish","zva","sy"])];J=o.choice(MNEMS);I=[o.choice(pool)()for _ in range(o.choice([0,1,2,2,3,3,3,4,4]))];return(J+" "+",".join(I)).strip()
def gen(o):
 if o.random()<0.5:return gen_shape(o)
 R=lambda:o.choice(["x1","x2","x3","x9","x30","sp","xzr"]);I=lambda:o.choice([o.randrange(-1100,1100),o.randrange(-1100,1100)//8*8,o.randrange(0,33000)//8*8,o.randrange(-300,300),32760,32768,504,512,-512,-520,255,256,-256,-257,0]);a,b,c,k=R(),R(),R(),o.randrange(25);J=lambda:o.choice([o.getrandbits(64),o.getrandbits(16)<<o.choice([0,16,32,48]),((1<<o.randrange(1,64))-1)<<o.randrange(0,8),o.randrange(0,5000),o.randrange(0,4096)<<12]);return[lambda:"ldp %s,%s,[%s,#%d]"%(a,b,c,I()),lambda:"stp %s,%s,[%s,#%d]!"%(a,b,c,I()),lambda:"ldp %s,%s,[%s],#%d"%(a,b,c,I()),lambda:"ldr %s,[%s,#%d]"%(a,c,I()),lambda:"str %s,[%s],#%d"%(a,c,I()),lambda:"ldr %s,[%s,#%d]!"%(a,c,I()),lambda:"ldr %s,[%s,%s,lsl #%d]"%(a,b,c,o.choice([0,3,2])),lambda:"%s %s,%s,#%d"%(o.choice(["add","sub","adds","subs"]),a,b,o.choice([o.randrange(-5000,5000),J()])),lambda:"%s %s,%s,%s,%s #%d"%(o.choice(["add","adds","sub","subs"]),a,b,c,o.choice(["lsl","lsr"]),o.randrange(0,66)),lambda:"%s %s,%s,#0x%x"%(o.choice(["and","orr","eor"]),a,b,J()),lambda:"tst %s,#0x%x"%(a,J()),lambda:"mov %s,#0x%x"%(a,J()),lambda:"mov %s,%s"%(a,b),lambda:"movk %s,#%d,lsl #%d"%(a,o.randrange(0,70000),o.choice([0,8,16,32,48])),lambda:"ubfx %s,%s,#%d,#%d"%(a,b,o.randrange(0,66),o.randrange(0,66)),lambda:"bfi %s,%s,#%d,#%d"%(a,b,o.randrange(0,66),o.randrange(0,66)),lambda:"bfxil %s,%s,#%d,#%d"%(a,b,o.randrange(0,66),o.randrange(0,66)),lambda:"cmp %s,#%d"%(a,o.randrange(-5000,5000)),lambda:"cmp %s,%s"%(a,b),lambda:"stlur %s,[%s,#%d]"%(a,c,I()),lambda:"ldadd %s,%s,[%s]"%(a,b,c),lambda:"stadd %s,[%s]"%(a,c),lambda:"%s %s,%s,%s"%(o.choice(["mul","umulh","adcs","sbc"]),a,b,c),lambda:"extr %s,%s,%s,#%d"%(a,b,c,o.randrange(0,66)),lambda:"lsr %s,%s,#%d"%(a,b,o.randrange(0,66))][k]()
def main():
 seed,n=(int(sys.argv[1])if len(sys.argv)>1 else 1),(int(sys.argv[2])if len(sys.argv)>2 else 4000);rg,L=random.Random(seed),[]
 while len(L)<n:L.append(gen(rg))
 a4=[]
 for ag in L:
  mn,_,a6=ag.partition(" ")
  try:encodable(mn,_a64_ops(a6));a4 .append(True)
  except AssertionError:a4 .append(False)
 with tempfile.TemporaryDirectory()as ae:
  ah=errors(L,ae);a7=[i not in ah for i in range(len(L))];J=[L[i]for i in range(len(L))if a4[i]and not a7[i]];I=[L[i]for i in range(len(L))if a7[i]and not a4[i]]
  def declared(x):
   mn,_,a6=x.partition(" ");a9=_a64_ops(a6)
   if mn=="ldp"and len(a9)>1 and a9[0]==a9[1]:return True
   if mn in("b","bl","cbz","cbnz","tbz","tbnz","dmb","wfe")or mn.startswith("b."):return True
   return any(o.strip()[:1]=="w"for o in a9)or mn not in _SHAPES and mn not in("dmb","dc")
  a5=[x for x in I if not declared(x)];print("seed %d: %d instructions, LLVM accepts %d, a64 accepts %d; a64 accepts but LLVM rejects: %d; a64 stricter: %d ""(declared scope: %d, other: %d)"%(seed,n,sum(a7),sum(a4),len(J),len(I),len(I)-len(a5),len(a5)))
  for d in(J+a5)[:40]:print("  ",d)
  ai=J;ab=[]
  for f in("leaf.s","route.s","secp_fe.s"):
   af=[l for l in expand(open(os.path.join(HERE,"kernels",f)).read()).split("\n")];a8=os.path.join(ae,f);open(a8,"w").write("\t.text\n"+"\n".join(af)+"\n");p=subprocess.run(CLANG+[a8,"-o",a8+".o"],capture_output=True,text=True);ac=sum(": error:"in l for l in p.stderr.split("\n"));dump=subprocess.run(([]if shutil.which("llvm-objdump")else["xcrun"])+["llvm-objdump","-d",a8+".o"],capture_output=True,text=True).stdout;ad=sum(1 for l in dump.split("\n")if re.match(r"\s+[0-9a-f]+:\s+(([0-9a-f]{2} ){4}|[0-9a-f]{8}\s)",l));aa=len(a64_load(open(os.path.join(HERE,"kernels",f)).read())[0]);print("%s: LLVM errors %d; LLVM instructions %d, a64_load instructions %d"%(f,ac,ad,aa))
   if ac or ad!=aa:ab.append(f)
 return len(ai)+len(ab)
if __name__=="__main__":sys.exit(1 if main()else 0)
