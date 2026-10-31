"""One isolated verifier run with two separately captured checks.
Reads only immutable files supplied by the verifier snapshot.
"""
import io,json,sys,contextlib
from pathlib import Path
cfg=json.load(sys.stdin)
def capture(filename, inputs):
    source=Path('/snapshot/files',filename).read_text()
    output=io.StringIO()
    prior=sys.stdin
    try:
        sys.stdin=io.StringIO(json.dumps(inputs))
        with contextlib.redirect_stdout(output):
            exec(compile(source,filename,"exec"),{"__name__":"__main__"})
    finally:
        sys.stdin=prior
    return json.loads(output.getvalue())
ledger=capture('accounting.py',{})
procedure=capture('experiment.py',cfg)
schedule=capture('forward_check.py',{})
print(json.dumps({"accounting":ledger,"procedure":procedure,"forward_schedule":schedule}))
