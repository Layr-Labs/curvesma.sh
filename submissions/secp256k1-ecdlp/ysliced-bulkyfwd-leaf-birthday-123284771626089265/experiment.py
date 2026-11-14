# Loader (proof.md 7): Python evidence code ships as one lzma stream in base85 text (code.b85). This stub restores
# experiment.py's full text byte-identically into a temporary directory with
# copies of the plain files, checks each file's SHA-256 against SHA (also listed in proof.md 7), and
# runs the restored experiment.py, which reports the check as "restored". Decoding is evidence tooling, not part of
# the algorithm, so it is not charged.
import base64,hashlib,json,lzma,os,runpy,shutil,sys,tempfile
BLOBS=['code.b85']
SHA={"experiment.py": "a6655df929c080b209b167b3b9708f7301ae603c548a154df4892db2427b207f"}
def pk():
 h=os.path.dirname(os.path.abspath(__file__))
 for b in(h,os.path.join(h,"files"),"/snapshot/files"):
  if os.path.isfile(os.path.join(b,BLOBS[0])):return b
 raise FileNotFoundError(BLOBS[0]+" not found")
B=pk();F=json.loads(lzma.decompress(base64.b85decode(b"".join(open(os.path.join(B,n),"rb").read()for n in BLOBS))))
R=tempfile.mkdtemp();os.makedirs(os.path.join(R,"kernels"))
for n in('certificate.py', 'kernels/a64.py', 'kernels/leaf.s', 'kernels/route.s', 'kernels/secp_fe.s', 'machine.py', 'proof.md', 'sites.py', 'submission.json', 'success.md'):shutil.copy(os.path.join(B,n),os.path.join(R,n))
H={}
for n,t in F.items():d=t.encode("utf-8");H[n]=hashlib.sha256(d).hexdigest();open(os.path.join(R,n),"wb").write(d)
json.dump(dict(files=H,ok=H==SHA),open(os.path.join(R,"restored.json"),"w"))
runpy.run_path(os.path.join(R,"experiment.py"),run_name="__main__")
