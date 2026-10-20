"""Challenge integration checks; do not execute candidates or call a provider."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("local_verifier", ROOT / "scripts/local_verifier.py")
local = importlib.util.module_from_spec(spec)
spec.loader.exec_module(local)


class LocalSetupTests(unittest.TestCase):
    def test_manifest_and_hosted_verifier_use_same_pin(self):
        for track in json.loads((ROOT / "benchmark.json").read_text())["tracks"]:
            name = track["name"]
            self.assertEqual(track["setupCommand"], ["python3", "scripts/local_verifier.py", "setup", name])
            self.assertEqual(track["benchmarkCommand"], ["python3", "scripts/local_verifier.py", "run", name])
            caller = (ROOT / ".github/workflows" / f"{name}.yml").read_text()
            self.assertIn(f"verify.yml@{local.REVISION}", caller)
            self.assertIn(f"verifier_revision: {local.REVISION}", caller)

    def invoke(self, files):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "benchmark.json").write_text(json.dumps({"tracks": [{"name": "secp256k1-ecdlp"}]}))
            output = root / ".yukon/secp256k1-ecdlp"
            output.mkdir(parents=True)
            for name, data in files.items():
                (output / name).write_text(json.dumps(data))
            (output / "score.json").write_text("stale score")
            stderr = io.StringIO()
            with patch.object(local, "ROOT", root), patch.object(local, "installed", return_value=True), patch.object(local, "run") as run, patch("sys.argv", ["local_verifier.py", "run", "secp256k1-ecdlp"]), contextlib.redirect_stderr(stderr):
                try:
                    local.main()
                    status = 0
                except SystemExit as error:
                    status = error.code
                self.assertFalse((output / "score.json").exists())
                return status, stderr.getvalue(), run.call_args

    def test_missing_review_is_actionable_and_clears_stale_score(self):
        status, message, call = self.invoke({"intake.json": {}})
        self.assertEqual(status, 2)
        self.assertIn("initial.json, review.json", message)
        self.assertIn("LOCAL_SETUP.md", message)
        self.assertIsNone(call)

    def test_direct_judgment_delegates_to_upstream_finalization(self):
        status, _, call = self.invoke({"intake.json": {}, "initial.json": {}, "review.json": {"data": {"decision": {"kind": "judgment"}}}})
        self.assertEqual(status, 0)
        self.assertEqual(call.args[0], "finalize")
        self.assertNotIn("--followup", call.args)

    def test_requested_followup_requires_both_artifacts(self):
        files = {"intake.json": {}, "initial.json": {}, "review.json": {"data": {"decision": {"kind": "followup"}}}}
        status, message, call = self.invoke(files)
        self.assertEqual(status, 2)
        self.assertIn("followup.json, final-review.json", message)
        self.assertIsNone(call)
        status, _, call = self.invoke({**files, "followup.json": {}, "final-review.json": {}})
        self.assertEqual(status, 0)
        self.assertIn("--followup", call.args)
        self.assertIn("--final-review", call.args)
