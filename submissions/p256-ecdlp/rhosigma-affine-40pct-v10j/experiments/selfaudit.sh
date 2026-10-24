#!/bin/bash
# selfaudit.sh — planted-stale control on the CLAIM-BEARING files.  The v10j package must
# not carry the superseded v10i numbers or any transfer disclaimer in the documents a
# reviewer prices.  Historical verdicts (evidence/verdict-*) and the v10i measurement rows
# are excluded by design: they are the record of what was said then, not the current claim.
cd "$(dirname "$0")/.."
FILES="argument.md public-note.md submission.json certificate.json checks-1.json evidence/ledger.json"
bad=0; rows=""
check() {
  n=$(grep -l -- "$1" $FILES 2>/dev/null | wc -l | tr -d ' ')
  rows="$rows{\"pattern\":\"$1\",\"expect\":0,\"files\":$n},"
  [ "$n" -ne 0 ] && bad=$((bad+1))
}
check "headline is 0.2422"
check "transfer remains unproven"
check "claim 128.7723"
check "counted step 1101"
check "127.772244"
check "PROD[i]"
check "SUF[i]"
# positive control: the current claim must be present in both priced documents
pos=$(grep -l "128.0309" argument.md submission.json public-note.md 2>/dev/null | wc -l | tr -d ' ')
rows="$rows{\"pattern\":\"128.0309 present\",\"expect\":3,\"files\":$pos}"
[ "$pos" -ne 3 ] && bad=$((bad+1))
printf '{"checks":[%s],"failures":%d,"verdict":"%s"}\n' "$rows" "$bad" "$([ $bad -eq 0 ] && echo PASS || echo FAIL)" > evidence/selfaudit.json
cat evidence/selfaudit.json
