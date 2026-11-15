# Loader: proof.md 7 (Packaging).
import base64,hashlib,json,lzma,os,runpy,shutil,tempfile
SHA="64a7bb22cb15cec18f21722de4188e5970c57c64080a266ccec7852892145ec1"
def pk():
 h=os.path.dirname(os.path.abspath(__file__))
 for b in(h,os.path.join(h,"files"),"/snapshot/files"):
  if os.path.isfile(os.path.join(b,"code.b85")):return b
 raise FileNotFoundError("code.b85 not found")
B=pk();D=lzma.decompress(base64.b85decode(open(os.path.join(B,"code.b85"),"rb").read()))
R=tempfile.mkdtemp();os.makedirs(os.path.join(R,"kernels"))
for n in('certificate.py', 'kernels/a64.py', 'kernels/leaf.s', 'kernels/route.s', 'kernels/secp_fe.s', 'machine.py', 'proof.md', 'sites.py', 'submission.json', 'success.md'):shutil.copy(os.path.join(B,n),os.path.join(R,n))
H=hashlib.sha256(D).hexdigest();open(os.path.join(R,"experiment.py"),"wb").write(D)
json.dump(dict(files={"experiment.py":H},ok=H==SHA),open(os.path.join(R,"restored.json"),"w"))
runpy.run_path(os.path.join(R,"experiment.py"),run_name="__main__")
