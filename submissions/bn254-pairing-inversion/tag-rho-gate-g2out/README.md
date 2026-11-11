# bn254-pairing-inversion: gated two-level tag-walk rho with pairing output (P0, [x]Q0)

Claimed time_log2 = 118.863948. This is the expected total work over uniform challenges and fresh coins:
one fresh coin opens the run or stops it at once. On every run the work is at most 2^129.347893 units.
Success is at least 0.39 under the disclosed heuristics H1, H2 and H3. Peak memory is at most 2^62.877284
bytes in one address space, and there is no advice.

The walk (121 instructions per step, 178 more for a second-level tag), store, ledger and evidence come
from the GT-DLP result 4cea9e8c (source 95a0ca34). The change is the output: after g^x = z, the
algorithm outputs (P0, [x]Q0) with a counted G2 scalar multiplication (proof.md section 0, from pairing
package 8ebb2b12), charged as an extra 2^22 in the output row. The parent's full fixed allowance is retained; the upward-rounded score remains 118.863948.

- [Algorithm, pairing output, tag lemmas, success analysis, ledger](proof.md)
- [Instruction counts and storage](machine.md)
- [Checks](experiment.py): the declared run; mode all also runs the G2 output check
- [G2 output](g2out.py), [earlier recorded output run](g2out-evidence.json)
- [Interpreter and kernels](kern.py), [toy walk](toy.py), [offline toy runs](toy_offline.json)
