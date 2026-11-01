"""Run inherited and transfer checks without hiding a failed inherited result."""
import contextlib
import io
import json
from pathlib import Path
import sys

here = Path(__file__).resolve().parent
files = here / "files"
sys.path.insert(0, str(files if (files / "evidence16.py").is_file() else here))
import evidence16
import transfer_check

def captured_run(function, input_text):
    buffer = io.StringIO()
    previous_stdin = sys.stdin
    try:
        sys.stdin = io.StringIO(input_text)
        with contextlib.redirect_stdout(buffer):
            function()
    finally:
        sys.stdin = previous_stdin
    return json.loads(buffer.getvalue())

if __name__ == "__main__":
    source = captured_run(evidence16.main, sys.stdin.read())
    transfer = captured_run(transfer_check.main, "{}")
    ok = source.get("checks_ok") is True and all(transfer["checks"].values())
    print(json.dumps({"checks_ok": ok, "source": source, "transfer": transfer}))
    if not ok:
        raise SystemExit(1)
