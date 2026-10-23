#!/usr/bin/env bash
set -euo pipefail
if [[ $# != 2 && $# != 3 ]] || [[ $# == 3 && "$3" != --preflight-only ]]; then
  echo "Usage: $0 <track> <submission-directory> [--preflight-only]" >&2
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
python=${PYTHON:-python3}
if [[ ${3:-} == --preflight-only ]]; then
  exec "$python" -m llm_verifier preflight --root "$repo" --config "tracks/$track/verifier.json" --submission "$source"
fi
output="$repo/.yukon/$track/supporting-code"
mkdir -p "$output"
"$python" -m llm_verifier intake --root "$repo" --config "tracks/$track/verifier.json" --submission "$source" --out "$output/intake.json"
"$python" -m llm_verifier prepare-runtime --intake "$output/intake.json"
"$python" -m llm_verifier experiments --intake "$output/intake.json" --out "$output/initial.json"
echo "Experiment results: $output/initial.json"
