"""Experiment entry point: the evidence code is evidence_core.py inside libextra.b85 (served to import by common12.py)."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE / "files" if (HERE / "files").is_dir() else HERE; sys.path.insert(0, str(ROOT))
import common12
from evidence_core import main
if __name__ == "__main__": main()
