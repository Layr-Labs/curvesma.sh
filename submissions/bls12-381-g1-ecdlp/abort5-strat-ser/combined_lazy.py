import contextlib,io,json,runpy,sys
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root))
import native_check
buf=io.StringIO()
with contextlib.redirect_stdout(buf):runpy.run_path(str(root/'combined.py'),run_name='__main__')
out=json.loads(buf.getvalue())
out['lazy_native']=native_check.run({'seed':113,'samples':300})
from fractions import Fraction
import log_bound
out['exact_log_bound']=log_bound.certify(out['ledger']['W_instructions'],2000,Fraction(116157101456605786,10**15))
import abort_check
_abort_cfg={"bits":23,"logd":5,"K":"0.50","trials":60,"seed":2026,"keybits":18,"nodebits":5,"headbits":6,"head_cap":64,"logX2":13,"logG":12,"draws":16,"seg":8,"max_tree":8,"strata":4}
out['abort_exec']=abort_check.run(_abort_cfg)
out['scope_new']='Full amended ledger and finite implementation checks. Production success remains conditional on explicitly disclosed source H1-H5; not production execution or independent mathematical certification.'
print(json.dumps(out))
