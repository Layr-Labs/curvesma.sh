"""The declared run. Input on stdin: {"seed": 1, "seconds": 12}. Output: one JSON object; exit status 1 if a check fails.

What runs, in this order:
  certificate  machine.py imports certificate.py: every check of the certificate is evaluated. The lines it prints
               (each check with its value, the names of the measured constants, the allowances with no kernel behind
               them) are copied into the output under "certificate".
  cases        cases.py: every transport kernel on the memory model of machine.py, path by path, with controls.
  leaf         leafx.py: the leaf kernels (gate, tables, psel, yround, fwd, fe_inv, bwd, keys) and the big-integer
               model of the label; frames reach the leaf kernels only through yprod, yrel and ym0.
  toy          toy.py: the end-to-end run in two processor orders, with controls.
  lemma_j      exponential sums over all points of a small curve (the shape of Lemma J's bound; not its count).
"constants": every measured constant of the certificate, with the certificate's value and the number of executed
counts of this run whose formula depends on it (a different value would fail that many comparisons).
"compared": the number of executed instruction counts compared with their formulas, by kernel.
"restored", "restored_sha256": for each file restored from code.b94 (unpack.py), whether its SHA-256 equals the value
listed in unpack.py, and the SHA-256 this run computed.
"""
import json, os, sys, time
T0 = time.time(); CFG = json.load(sys.stdin); SEED = int(CFG.get("seed", 1))
def pkg():
    here = os.path.dirname(os.path.abspath(__file__))
    for base in (here, os.path.join(here, "files"), "/snapshot/files"):
        if os.path.isfile(os.path.join(base, "machine.py")): return base
    raise FileNotFoundError("machine.py not found")
HERE = pkg(); sys.path.insert(0, HERE)
import unpack; RESTORED, SHA256 = dict(unpack.OK), dict(unpack.GOT)
import machine as MC
MC.SEED = SEED; K = MC.K; T = [time.time()]
import cases as C
T.append(time.time())
import leafx, toy
L = leafx.run(SEED); T.append(time.time())
E = toy.end_to_end(SEED, leafx); T.append(time.time())
J = leafx.lemma_j(SEED); T.append(time.time())
O = C.finish()
OUT = dict(input=CFG, certificate=dict(passed=not K.FAILED, failed=K.FAILED, claim="%d.%015d" % divmod(int(K.CLAIM_E * 10**15), 10**15), lines=MC.CERT_TEXT.strip().split("\n")),
           constants=O["constants"], cases={k: v["checks"] if v["ok"] else v for k, v in O["cases"].items()}, controls=O["controls"], counts=O["counts"], compared=O["compared"],
           stat=O["stat"], leaf=L, toy=E, lemma_j=J, restored=RESTORED, restored_sha256=SHA256)
OUT["all_checks_pass"] = bool(not K.FAILED and O["ok"] and len(O["controls"]) == 19 and all(v[1] for v in O["constants"].values()) and L["ok"] and J["ok"] and (RESTORED is None or all(RESTORED.values())))
OUT["seconds"] = round(time.time() - T0, 2); OUT["seconds_by_part"] = [round(b - a, 2) for a, b in zip([T0] + T, T)]
json.dump(OUT, sys.stdout, separators=(",", ": "), default=lambda x: int(x) if x == int(x) else str(x)); print()
if not OUT["all_checks_pass"]: sys.exit(1)
