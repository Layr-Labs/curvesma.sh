# Complete source + SHA256: stats16.b85/source_archive.
from pathlib import Path
import sys
_h=Path(__file__).resolve().parent
sys.path.insert(0,str(_h/'files'if(_h/'files').is_dir()else _h))
from owner10 import _load_source
_load_source(globals(),'counts16.py','stats16.b85')
