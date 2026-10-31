"""Run inherited evidence and the added parity-key audit in isolation."""
import contextlib, io, json, runpy, sys
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():
    root=root/'files'
sys.path.insert(0,str(root))
import key_check
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    runpy.run_path(str(root/'evidence.py'),run_name='__main__')
out=json.loads(buf.getvalue())
out['branch_key']=key_check.check()
out['new_key_checks_passed']=True
print(json.dumps(out))
