#!/usr/bin/env bash
set -euo pipefail
if [[ $# != 2 ]]; then
  echo "Usage: $0 <track> <submission-directory>" >&2
  exit 2
fi
repo=$(cd "$(dirname "$0")/.." && pwd)
case "$1" in
  secp256k1-ecdlp|bls12-381-invert-g1|bls12-381-g1-ecdlp|bls12-381-miller-inversion) track=$1 ;;
  *) echo "Unknown track: $1" >&2; exit 2 ;;
esac
source=$(cd "$2" && pwd)
case "$source" in
  "$repo"/submissions/"$track"|"$repo"/submissions/"$track"/*) source=${source#"$repo"/} ;;
  *) echo "Submission must be under submissions/$track/" >&2; exit 2 ;;
esac
output="$repo/.yukon/$track/supporting-code"
mkdir -p "$output"
python=${PYTHON:-"$repo/.yukon/verifier/venv/bin/python"}
if ! command -v "$python" >/dev/null 2>&1; then
  echo "Verifier missing. Run yukon setup --track $track first (see tracks/docs/LOCAL_SETUP.md)." >&2
  exit 2
fi
"$python" -m llm_verifier intake --root "$repo" --config "tracks/$track/verifier.json" --submission "$source" --out "$output/intake.json"
"$python" -m llm_verifier prepare-runtime --intake "$output/intake.json"
"$python" -m llm_verifier experiments --intake "$output/intake.json" --out "$output/initial.json"
echo "Experiment results: $output/initial.json"
