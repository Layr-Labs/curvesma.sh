"""Independent structural checks; no arithmetic-ledger parser cache."""
import json,sys,time
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root))
import evidence
cfg=json.loads(sys.stdin.read()); start=time.monotonic()
out={"mode":"structure"}
out["planner"]=evidence.planner_mode(cfg["planner"])
out["ec"]=evidence.ec_attack(cfg["ec"])
out["proto"]=evidence.proto_mode(cfg.get("proto",{}))
out["census"]=[evidence.census(dict(cfg["census"]))]
out["h2"]=evidence.h2_sizes()
out["seconds"]=round(time.monotonic()-start,2)
print(json.dumps(out))
