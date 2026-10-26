#!/bin/bash
cd "$(dirname "$0")"
for spec in "20 0 251 7 5" "20 0 4099 1234 77" "21 0 65537 4321 999" "21 0 1048583 777 31337" "22 0 16777259 55555 1234567" "22 0 67108879 999983 4242"; do ./slice8 check $spec 2>/dev/null | tail -1; done
./slice8 check 26 0 4611686018427388039 999983 4242 2>/dev/null | tail -1
