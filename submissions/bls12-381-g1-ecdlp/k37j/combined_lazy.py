import contextlib,io,json,runpy,sys
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root))
import kernel_check, nocarry_witness, batch_audit
buf=io.StringIO()
with contextlib.redirect_stdout(buf):runpy.run_path(str(root/'evidence.py'),run_name='__main__')
out=json.loads(buf.getvalue())
out['kernel_check']=kernel_check.run({'seed':113,'samples':200})
out['kernel_check_plus']=kernel_check.run({'seed':113,'samples':200,'plus':True})
out['leaf_kernel']=kernel_check.run_leaf({'seed':113,'samples':200})
out['leaf_kernel_plus']=kernel_check.run_leaf({'seed':113,'samples':200,'plus':True})
out['beta_keys']=batch_audit.beta_run({'seed':11})
out['table_level']=batch_audit.run({'table':True,'seed':7})
out['nocarry_witness']=nocarry_witness.run({'seed':5,'samples':2000})
from fractions import Fraction
import log_bound
out['exact_log_bound']=log_bound.certify(out['ledger']['W_instructions'],2000,Fraction(115893955025065712,10**15))
out['scope_new']='Amended ledger + finite checks.'
print(json.dumps(out))
