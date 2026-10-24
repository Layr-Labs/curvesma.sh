#!/bin/bash
# preflight_v10p.sh -- intake laundry, measured not guessed.
# The ARCHIVE is the WHOLE track submissions dir, so this measures the tar expansion of that dir,
# not the package size.  Cap 262144 B; target <= 240 KB for margin.
set -u
TRACK=/Users/andrewgordon/Documents/ChatGPT/curvesma.sh/curvesmash/submissions/p256-ecdlp
PKG=$TRACK/rhosigma-w2measured-40pct-v10p
echo "== files in package: $(find $PKG -type f | wc -l | tr -d ' ') (cap 30)"
echo "== pycache:"; find $TRACK -name __pycache__ -o -name '*.pyc' | head
find $TRACK -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null
find $TRACK -name '*.pyc' -delete 2>/dev/null
echo "== track dir listing:"; ls $TRACK
echo "== tar expansion of the whole track dir:"
EXP=$(tar -cf - -C "$TRACK" . 2>/dev/null | (cd /tmp && mkdir -p tarprobe && cd tarprobe && tar -xf - && du -sk . | awk '{print $1*1024}'))
echo "   expands to $EXP bytes (cap 262144, target 245760)"
python3 - "$EXP" <<'EOF'
import sys
exp=int(sys.argv[1]); cap=262144
print("   OK" if exp<=cap else "   OVER CAP by %d"%(exp-cap), "margin %d B"%(cap-exp))
EOF
echo "== title length:"; python3 - <<'EOF'
import json
d=json.load(open("/Users/andrewgordon/Documents/ChatGPT/curvesma.sh/curvesmash/submissions/p256-ecdlp/rhosigma-w2measured-40pct-v10p/submission.json"))
print("   title", len(d["title"]), "claims[0]", len(d["claims"][0]["statement"]), "keys", sorted(d.keys()))
EOF
echo "== evidence files referenced vs present:"
python3 - <<'EOF'
import json, os
pkg="/Users/andrewgordon/Documents/ChatGPT/curvesma.sh/curvesmash/submissions/p256-ecdlp/rhosigma-w2measured-40pct-v10p"
d=json.load(open(os.path.join(pkg,"submission.json")))
miss=[f for f in d["claims"][0]["evidence"] if not os.path.exists(os.path.join(pkg,f))]
have=set()
for r,_,fs in os.walk(pkg):
    for f in fs: have.add(os.path.relpath(os.path.join(r,f),pkg))
extra=sorted(have-set(d["claims"][0]["evidence"]))
print("   missing:", miss)
print("   present but unlisted:", extra)
EOF
