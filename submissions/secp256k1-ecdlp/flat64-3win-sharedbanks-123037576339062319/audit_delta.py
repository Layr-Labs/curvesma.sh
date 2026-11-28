"""Chain of custody for the LLVM instruction audit (encoding-audit.json) of kernels/leaf.s and kernels/secp_fe.s.

The LLVM audit (llvm-mc, flags in the json) covers the seven source templates of 8cc64cb7, whose sha256 values are
input_kernel_sha256. Five of them are shipped unchanged (generic.s, pointer_init.s, route.s, ystage.s, ytransport.s). leaf.s and
secp_fe.s are shipped changed (three-window table lookups through the T aperture). Their audited versions are restored from
runtime3.json as kernels/base_leaf.s and kernels/base_secp_fe.s. This module checks their hashes, recomputes the multiset difference
of the expanded instruction texts (branch labels normalised: their names differ per macro expansion), checks it against the lists
recorded in encoding-audit.json, and checks that the FORM (mnemonic and operand kinds) of every new text occurs among the
instructions of the audited sources. The new texts were not assembled by LLVM here: the form check, the encoder checks of a64.py
(immediate ranges, applied when a source is loaded) and the executions in experiment.py are the encoding evidence for them
(encoding-audit.json: shared_t_delta.llvm_status)."""
import hashlib, re
from collections import Counter

BASE_FILES = ("leaf.s", "secp_fe.s")
_BR = {"b", "bl", "cbz", "cbnz", "tbz", "tbnz", "adr", "cbhs", "cbhi", "cbeq", "cbne", "cblt", "cbge"}
_CONDS = {"eq", "ne", "hs", "lo", "cs", "cc", "mi", "pl", "vs", "vc", "hi", "ls", "ge", "lt", "gt", "le", "al", "nv"}


def _split(ins):
    mn, _, rest = ins.partition(" ")
    return mn, [x.strip() for x in re.split(r",(?![^\[]*\])", rest) if x.strip()]


def norm(ins):
    mn, ops = _split(ins)
    if (mn in _BR or mn.startswith("b.")) and ops:
        ops[-1] = "LABEL"
    return mn + (" " + ",".join(ops) if ops else "")


def form(ins):
    mn, ops = _split(ins)
    out = []
    for o in ops:
        if re.fullmatch(r"x\d+|xzr|sp", o): out.append("R")
        elif re.fullmatch(r"#-?(0x[0-9a-fA-F]+|\d+)", o): out.append("#I")
        elif o == "LABEL" or re.fullmatch(r"\d+[fb]", o): out.append("L")
        elif o.startswith("["): out.append("[M]")
        elif re.fullmatch(r"(lsl|lsr|asr|ror) #\d+", o): out.append("S")
        elif re.fullmatch(r"[A-Za-z_.][A-Za-z0-9_.]*", o) and o not in _CONDS: out.append("L")
        else: out.append(o)
    return mn + " " + ",".join(out)


_NORM = {}


def instrs(text, a64):
    out = []
    for l in a64.expand(text).split("\n"):
        l = l.split("//")[0].strip()
        if l and not l.endswith(":") and not l.startswith("."):
            n = _NORM.get(l)
            if n is None:
                n = _NORM[l] = norm(l)
            out.append(n)
    return out


def delta(base_text, ship_text, a64):
    ib, ish = instrs(base_text, a64), instrs(ship_text, a64)
    cb, cs = Counter(ib), Counter(ish)
    new = sorted(t for t in cs if cs[t] > cb.get(t, 0))
    removed = sorted(t for t in cb if cb[t] > cs.get(t, 0))
    return ib, ish, new, removed


def check(audit, K, a64, assembled_bytes, materialize):
    d = audit["shared_t_delta"]
    texts = {n: (K / n).read_text() for n in ("generic.s", "route.s", "pointer_init.s", "ystage.s", "ytransport.s")}
    others = [materialize(texts["generic.s"] + "\n" + texts["route.s"]), texts["pointer_init.s"], texts["ystage.s"] + "\n" + texts["ytransport.s"]]
    forms = set()
    for t in others: forms |= {form(i) for i in set(instrs(t, a64))}
    result = {}
    for source_name, expected in audit["input_kernel_sha256"].items():
        if source_name not in BASE_FILES:
            assert hashlib.sha256((K / source_name).read_bytes()).hexdigest() == expected, source_name
            continue
        base = (K / ("base_" + source_name)).read_text(); ship = (K / source_name).read_text()
        assert hashlib.sha256(base.encode()).hexdigest() == expected == d["files"][source_name]["base_sha256"], source_name
        assert hashlib.sha256(ship.encode()).hexdigest() == d["files"][source_name]["shipped_sha256"], source_name
        ib, ish, new, removed = delta(base, ship, a64)
        rec = d["files"][source_name]
        assert len(ib) == rec["base_insns"] == next(t for t in audit["templates"] if t["name"] == source_name.replace(".s", ""))["insns"], source_name
        assert len(ish) == rec["shipped_insns"], source_name
        assert new == rec["new_instruction_texts"], "stale new-text list: " + source_name
        assert removed == rec["removed_instruction_texts"], "stale removed-text list: " + source_name
        forms |= {form(i) for i in set(ib)}
        result[source_name] = dict(base_insns=len(ib), shipped_insns=len(ish), new_texts=len(new), removed_texts=len(removed))
    for source_name, rec in d["files"].items():
        missing = [t for t in rec["new_instruction_texts"] if form(t) not in forms]
        assert not missing, (source_name, missing)
        result[source_name]["new_forms_all_in_audited_sources"] = True
    assert audit["assembled_bytes_per_endpoint"] == d["audited_bytes_per_endpoint"] != assembled_bytes == d["assembled_bytes_per_endpoint_shipped"]
    return result
