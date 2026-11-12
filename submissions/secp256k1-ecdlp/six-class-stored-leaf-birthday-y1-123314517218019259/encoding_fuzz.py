import os,subprocess,sys,tempfile
H=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,H+"/kernels")
from a64 import expand,a64_load
C=["clang","-target","aarch64-linux-gnu","-march=armv8.7-a+ls64+lse","-c"]
bad=0
with tempfile.TemporaryDirectory()as d:
 for f in("leaf.s","route.s","secp_fe.s"):
  t=open(H+"/kernels/"+f).read();a64_load(t);o=d+"/"+f;open(o,"w").write("\t.text\n"+expand(t)+"\n")
  p=subprocess.run(C+[o,"-o",o+".o"],capture_output=True,text=True);bad+=p.returncode!=0;print(f,"rc",p.returncode,p.stderr[:300])
sys.exit(1 if bad else 0)
