#!/usr/bin/env bash
set -euo pipefail
if [[ $# != 2 ]]; then
  echo "Usage: $0 <track> <submission-directory>" >&2
  exit 2
fi
repo=$(cd "$(dirname "$0")/.." && pwd)
case "$1" in
  secp256k1-ecdlp|bls12-381-g1-ecdlp|bls12-381-fixed-g1-invert-g2|bls12-381-fixed-g2-invert-g1|bls12-381-g1-ecdlp-plain|p256-ecdlp) track=$1 ;;
  *) echo "Unknown track: $1" >&2; exit 2 ;;
esac
source=$(cd "$2" && pwd)
case "$source" in
  "$repo"/submissions/"$track"|"$repo"/submissions/"$track"/*) source=${source#"$repo"/} ;;
  *) echo "Submission must be under submissions/$track/" >&2; exit 2 ;;
esac
output="$repo/.yukon/$track/supporting-code"
mkdir -p "$output"
python=${PYTHON:-python3}
"$python" -m llm_verifier intake --root "$repo" --config "tracks/$track/verifier.json" --submission "$source" --out "$output/intake.json"
"$python" -m llm_verifier prepare-runtime --intake "$output/intake.json"
"$python" -m llm_verifier experiments --intake "$output/intake.json" --out "$output/initial.json"
echo "Experiment results: $output/initial.json"
