#!/usr/bin/env python3
"""Render the prose artifacts from render.json (built by build_certificate.py).

Single source of truth: every number in public-note.md, argument.md and
submission.json is substituted from render.json.  Any token left unfilled is a
hard error, so no hand-typed number can survive.
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
r = json.load(open(HERE + "/render.json"))
val = {k.upper(): str(v) for k, v in r.items()}

for src, dst in (("submission.json.tpl", "submission.json"),
                 ("public-note.md.tpl", "public-note.md"),
                 ("argument.md.tpl", "argument.md")):
    txt = open(HERE + "/" + src).read()
    for k, v in val.items():
        txt = txt.replace("__%s__" % k, v)
    left = sorted(set(re.findall(r"__[A-Z0-9_]+__", txt)))
    if left:
        print("UNFILLED in %s: %s" % (src, left))
        sys.exit(1)
    open(HERE + "/" + dst, "w").write(txt)
    print("%-18s -> %-18s %d bytes" % (src, dst, len(txt)))

note = open(HERE + "/public-note.md").read()
if len(note) < 6000:
    print("NOTE TOO SMALL"); sys.exit(1)
if note.count("Model:") != 1 or note.count("Harness:") != 1:
    print("MODEL/HARNESS HEADER NOT EXACTLY ONCE"); sys.exit(1)
for ln, expect in ((1, "Model: qwen3.8-flash-next"), (2, "Harness: Hermes")):
    if note.split("\n")[ln - 1].strip() != expect:
        print("HEADER LINE %d WRONG: %r" % (ln, note.split("\n")[ln - 1])); sys.exit(1)
print("render ok")
