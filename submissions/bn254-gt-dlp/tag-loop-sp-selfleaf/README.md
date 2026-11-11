# Gated tag-tracing rho on BN254 GT with a 213-instruction step loop

time_log2 = **119.65** (11965/100). This is an **expected** total work bound; every run is capped at
2^126.68 units. Success is at least 0.39.

The package derives from e3771b98 (source b3aa659a). It changes the step loop and the distinguished-point store, and corrects the mean accounting. Parent scripts retain the same syntax trees and comments; indentation is compacted to fit the archive limit. The step loop goes from 217 to 213 instructions per certified trip, so a step costs
249 instead of 253:
- the limbs are read SP-relative, which leaves 22 resident;
- there are no scratch registers;
- window leaves point to themselves, and the window is counted, so there is no null-child test.

- [proof.md](proof.md): sections "Parent and changes", 5 and 7
- [machine.md](machine.md): "Step loop v2"
- [ttloop2.py](ttloop2.py): the declared initial run (loop, full driver at W=4, store, ledger); recorded output in ttloop2-evidence.json
- Distinguished-point store: a deterministic crit-bit trie keyed by the full canonical element (3,072 bits), with no hashing and no probe-failure term (machine.md)
- [lazy.py](lazy.py) / walks.py: the parent's experiment, kept as attributed historical evidence and not declared (its parts C/X assert the parent's figures); recorded outputs in development-evidence*.json
