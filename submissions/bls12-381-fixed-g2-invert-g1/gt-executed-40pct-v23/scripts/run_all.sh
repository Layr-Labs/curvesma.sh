#!/bin/bash
# run_all.sh -- v23 measurement sweep (background)
set -x
cd "$(dirname "$0")"
python3 measure23.py '{"mode":"census","fixtures":512,"procs":16,"out":"census23.json"}'
python3 measure23.py '{"mode":"twoscale","fixtures":512,"procs":16,"out":"twoscale23.json"}'
python3 measure23.py '{"mode":"scan","Rg":65536,"depth":5,"segs":64,"steps":160000,"procs":16,"out":"scan_prod23.json"}'
python3 measure23.py '{"mode":"scan","Rg":4096,"depth":1,"segs":40,"steps":160000,"procs":16,"out":"scan_v20_23.json"}'
python3 measure23.py '{"mode":"scan","Rg":12,"depth":1,"class_only":true,"segs":40,"steps":160000,"procs":16,"out":"scan_classonly23.json"}'
python3 measure23.py '{"mode":"survival","Rg":[128,1024,8192],"pairs":400,"steps":40,"gaps":[1,5,20,40],"procs":9,"out":"survival23.json"}'
echo ALLDONE
