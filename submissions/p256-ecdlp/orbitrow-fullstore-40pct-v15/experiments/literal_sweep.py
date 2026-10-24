"""Literal sweep: every integer of 5+ digits in the package must own an evidence value."""
import json, os, re, sys

D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEXT = ["public-note.md", "argument.md", "submission.json"]
EV = os.path.join(D, "evidence")


def main():
    blob = []
    for name in os.listdir(EV):
        if name.endswith((".json", ".txt")):
            blob.append(open(os.path.join(EV, name), errors="ignore").read())
    blob = "\n".join(blob)
    rows, unowned = {}, []
    for name in TEXT:
        txt = open(os.path.join(D, name)).read()
        for m in re.finditer(r"(?<!\d)(\d{5,})(?!\d)", txt):
            num = m.group(1)
            if num in rows:
                continue
            if num in blob:
                rows[num] = "evidence"
            else:
                rows[num] = "UNOWNED"
                unowned.append((name, num))
    ev = {"literals": rows, "unowned": unowned, "pass": not unowned,
          "rule": "every integer >= 5 digits in the public text reconciles to an evidence value"}
    json.dump(ev, open(os.path.join(EV, "literal_sweep.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps(ev, indent=1))
    sys.exit(0 if not unowned else 1)


if __name__ == "__main__":
    main()
