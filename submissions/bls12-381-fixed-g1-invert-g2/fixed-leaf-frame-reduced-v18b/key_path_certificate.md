# Current complete native key paths

Only the internal exact sign normal form and ordinary helper implementations changed from the source;the full192-byte key still has exact twelve-Frobenius fibers. key_contract.md gives the all-input normal-form/fiber theorem. S is now transient worker scratch,not a stored payload,so production correctness uses actualGT exponent verification as in proof.md.

## Exhaustive scalar sign paths

The original prefix loads Br and ORs its six limbs into x14;only if zero does it load Bi into x2..x7. It retains this first-nonzero component and zero flag. Load p once,compute the low-word p-a by64-bit SUB,compare unsigned with a_low,and copy or canonically negate BOTH components. p_low odd prevents a nonzero tie;bothzero takes the canonical-zero copy. All chosen source limbs are loaded before their stores;exact alias/disjoint contracts and all callee/SP/LR preservation hold.

Br nonzero prefix17,zero-flag branch1,real negate/store9,imag load/test9,nonzero imag negate/store/RET10 totals46. If Bi zero,the final block is4,so40. Br zero/high Bi is20+1+3 zero real stores+6 subtract+3 stores+RET=34. Low-copy paths<=33. These cases partition ALL canonical Fp2,not favorable random paths. Thus46 throughRET is a safe deterministic bound. The actual source makes no half-modulus decision,though the unused public constant remains harmless data.

## Odd transform and complete parent

The odd helper saves/restores only its used x19..x22,x29/LR and maps arguments directly into working registers. Its frame/scratch layout remains336 and raw input immutable. Direct canonical negation loads all six limbs,OR-detects zero,stores zero or subtracts from full p. Its nonzero21/zero13 paths replace zero-buffer/general subtraction;no zero becomes p. The two remaining add/sub helpers each31 and two constructed fixed-s lookups each202 are unchanged. Own longest61 plus62+404=527 throughRET (zero paths shorter). Every called helper preserves the ordinary ABI.

The unchanged complete parent canonical_twelve_bfirst.S has an acyclic own longest path110,counting BL as one own instruction;retain111. It includes all prologue/calls/comparisons,even/odd selection,full192 output bytes,identity abort and epilogue. Combining safe independent helper upper bounds gives111+2*46+527=730. At most8 caller/identity instructions738,retained741 with3 spare. Identity output is invalid and not read;all nonidentity outputs initialize all192 bytes. No pair probability or sampled-path independence is assumed for these bounds.

Current independent shared-only interpreter checks actual constants,all conditional carry/zero operations,initialized reads,permitted writes,immutable rawPhi,ABI,immediates,sign<=46,odd<=527,fixed202,key<=730 against its independent full-p formula,including identity and field boundary cases. Its scoped nonzero field fixtures are not asserted to be actualGT samples. They do not replace the twelve-fiber theorem or certify production success,cost,table construction,hardware or qualification. Full recurrence2623/map2387 and full183-bit NH346 checks remain separate. Global work is only the current proof/architecture ledger;older benchmark totals are not claims of this certificate.
