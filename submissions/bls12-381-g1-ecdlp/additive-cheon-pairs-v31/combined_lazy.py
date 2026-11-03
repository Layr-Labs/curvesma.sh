import contextlib,io,json,runpy,sys
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root))
import kernel_check
buf=io.StringIO()
with contextlib.redirect_stdout(buf):runpy.run_path(str(root/'combined.py'),run_name='__main__')
out=json.loads(buf.getvalue())
out['kernel_check']=kernel_check.run({'seed':113,'samples':200})
from fractions import Fraction
import log_bound
out['exact_log_bound']=log_bound.certify(out['ledger']['W_instructions'],2000,Fraction(116277611990909175,10**15))
import h6_ladder
out['h6_hosted']=[h6_ladder.run_joint(8,8000,8,0,0),h6_ladder.run_joint(8,8000,8,0,8)]
out['scope_new']='Full amended ledger and finite implementation checks. Production success remains conditional on explicitly disclosed H1-H6; not production execution or independent mathematical certification.' 
print(json.dumps(out))
