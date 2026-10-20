"""Local dependency bootstrap and finalization; review/execution remain upstream phases."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import venv

ROOT = Path(__file__).resolve().parents[1]
REVISION = "64b3a2b01f21f1281383c673503964d2fc4c5c47"
REPOSITORY = "https://github.com/Layr-Labs/llm-verifier.git"
ENVIRONMENT = ROOT / ".yukon" / "verifier" / "venv"
PYTHON = ENVIRONMENT / "bin" / "python"


def run(*arguments):
    subprocess.run([str(PYTHON), "-m", "llm_verifier", *map(str, arguments)], cwd=ROOT, check=True)


def installed():
    if not PYTHON.is_file():
        return False
    check = subprocess.run([str(PYTHON), "-c", "from importlib.metadata import distribution; print(distribution('llm-verifier').read_text('direct_url.json'))"], capture_output=True, text=True)
    if check.returncode:
        return False
    try:
        metadata = json.loads(check.stdout)
        return metadata.get("url") == REPOSITORY and metadata.get("vcs_info", {}).get("commit_id") == REVISION
    except (ValueError, AttributeError):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["setup", "run"])
    tracks = json.loads((ROOT / "benchmark.json").read_text())["tracks"]
    parser.add_argument("track", choices=[track["name"] for track in tracks])
    args = parser.parse_args()
    output = ROOT / ".yukon" / args.track
    if args.phase == "setup":
        if sys.version_info < (3, 11):
            parser.exit(2, "Setup requires Python 3.11 or newer.\n")
        if not installed():
            venv.EnvBuilder(with_pip=True).create(ENVIRONMENT)
            subprocess.run([str(PYTHON), "-m", "pip", "install", f"llm-verifier @ git+{REPOSITORY}@{REVISION}"], check=True)
        run("intake", "--root", ROOT, "--config", f"tracks/{args.track}/verifier.json", "--out", output / "intake.json")
        print("Intake complete. No model review or candidate execution was performed.", flush=True)
        print("Next: follow tracks/docs/LOCAL_SETUP.md for the separate execution and review phases.")
        return
    # Do not leave a previous result looking like the outcome of a failed new run.
    for name in ("score.json", "dossier.json"):
        (output / name).unlink(missing_ok=True)
    if not installed():
        parser.exit(2, f"Pinned verifier missing. Run yukon setup --track {args.track} first.\n")
    required = ["intake.json", "initial.json", "review.json"]
    missing = [name for name in required if not (output / name).is_file()]
    if missing:
        parser.exit(2, f"Cannot finalize {args.track}: missing {', '.join(missing)}.\nSetup only creates intake.json; run does not execute experiments or call a reviewer.\nFollow tracks/docs/LOCAL_SETUP.md, or submit for hosted evaluation.\n")
    # The verifier validates artifact contents and bindings; this only selects its CLI flags.
    review = json.loads((output / "review.json").read_text())
    extra = []
    if review.get("data", {}).get("decision", {}).get("kind") == "followup":
        missing = [name for name in ["followup.json", "final-review.json"] if not (output / name).is_file()]
        if missing:
            parser.exit(2, f"Review requested follow-up; missing {', '.join(missing)}. See tracks/docs/LOCAL_SETUP.md.\n")
        extra = ["--followup", output / "followup.json", "--final-review", output / "final-review.json"]
    run("finalize", "--intake", output / "intake.json", "--initial", output / "initial.json",
        "--initial-review", output / "review.json", "--score", output / "score.json",
        "--dossier", output / "dossier.json", *extra)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        sys.exit(error.returncode)
