# Assembles each kernel (a64.py expand + a64_load checks) with llvm-mc. llvm-mc 18 lacks +cmpbr, so the 3 CMPBR
# lines (cbhs/cbhi in leaf.s) are assembled only as cmp + b.hs/b.hi stand-ins: no CMPBR encoding is checked here.
import os,re,subprocess as P,sys,tempfile
H=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,H+"/kernels")
from a64 import expand,a64_load
b=0
with tempfile.TemporaryDirectory()as d:
 for f in("leaf.s","route.s","secp_fe.s"):
  t=open(H+"/kernels/"+f).read();a64_load(t);o=d+"/"+f;open(o,"w").write("\t.text\n"+re.sub(r"(?m)^(\s*)cb(h[si]) (\w+),(\w+),",r"\1cmp \3,\4\n\1b.\2 ",expand(t))+"\n");p=P.run(["llvm-mc","-triple=aarch64","-mattr=+lse,+rand,+ls64,+rcpc-immo","-filetype=obj",o,"-o",o+".o"],capture_output=True,text=True);b|=p.returncode;print(f,p.returncode,p.stderr[:300])
sys.exit(b and 1)
