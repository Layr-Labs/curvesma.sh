# Appendices (compact)
Parameters/ledger figures match `flip_ledger.py` / `probe.py` output (manifest_matches_ledger).
- Claim: 124761778695793762/10^15/10^15 = 2^124.761778695793762 group units, E=36060994515067490625304289300712755165625 scalar instructions.
- BWDP=1925/2, BWDY=7, BWDC=-229, bwd_count(MP=512,NY)=1925*MP*NY//2+7*NY-229.
- Label error [-21,12], 33 bad residues, 175-bit key, K2 mod 2^58, CAP_PG=90112, CAP_H=ceil(m_h)+2^41, gate K=7246216233269738109.
- Backward-pass layout: fully unrolled MP=512 pairs/Y, x19=32 hoisted once, no end-pointer loop/reload; borrow charged every pair.
- Success/capacity/tail/cold/source-recovery rows as in flip_ledger.py and probe checks. See code and proof.md.
