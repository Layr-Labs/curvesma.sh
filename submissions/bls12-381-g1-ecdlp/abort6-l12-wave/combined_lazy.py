import contextlib,io,json,runpy,sys
from pathlib import Path
from fractions import Fraction
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root))
import native_check,log_bound,abort_check
buf=io.StringIO()
with contextlib.redirect_stdout(buf):runpy.run_path(str(root/'combined.py'),run_name='__main__')
out=json.loads(buf.getvalue())
out['lazy_native']=native_check.run({'seed':113,'samples':300})
out['exact_log_bound']=log_bound.certify(out['ledger']['W_instructions'],2000,Fraction(116157807533101552,10**15))
out['abort_exec']=abort_check.run({"bits":23,"logd":5,"K":"0.50","trials":60,"seed":2026,"keybits":18,"nodebits":5,"headbits":6,"head_cap":64,"logX2":13,"logG":12,"draws":16,"seg":8,"max_tree":8,"strata":4})
out['scope_new']='Amended ledger and finite checks. Production success conditional on disclosed H1-H5.'
print(json.dumps(out))
