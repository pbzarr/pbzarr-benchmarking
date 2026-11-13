"""Cache-controlled trial timer: pixi run python -m unittest tests.test_timed_trials -v."""

import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
TIMER = ROOT / "workflows/lib/timed_trials.py"


def run_timer(*args):
    return subprocess.run(
        [sys.executable, str(TIMER), *map(str, args)],
        capture_output=True, text=True,
    )


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TimedTrialsTest(unittest.TestCase):
    def setUp(self):
        workspace = tempfile.TemporaryDirectory(prefix="timed-trials-test-")
        self.addCleanup(workspace.cleanup)
        self.dir = Path(workspace.name)
        self.data = self.dir / "data"
        (self.data / "nested").mkdir(parents=True)
        (self.data / "nested" / "chunk").write_bytes(b"x" * 4096)
        self.calls = self.dir / "calls"
        self.out = self.dir / "bench.tsv"
        # Each invocation appends one line, so the file counts runs.
        self.command = f"echo run >> {self.calls}; sleep 0.2"

    def test_warm_runs_one_untimed_warmup_then_records_each_trial(self):
        result = run_timer("--data", self.data, "--cache", "warm", "--trials", 3,
                           "--output", self.out, "--", "sh", "-c", self.command)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(self.calls.read_text().splitlines()), 4)
        rows = read_rows(self.out)
        self.assertEqual([row["trial"] for row in rows], ["1", "2", "3"])
        self.assertEqual({row["cache"] for row in rows}, {"warm"})
        for row in rows:
            self.assertGreaterEqual(float(row["s"]), 0.2)
            self.assertGreater(float(row["max_rss"]), 0)
            for field in ("io_in", "cpu_time", "mean_load"):
                self.assertGreaterEqual(float(row[field]), 0)

    @unittest.skipUnless(hasattr(os, "posix_fadvise"), "cold eviction needs posix_fadvise (Linux)")
    def test_cold_runs_no_warmup(self):
        result = run_timer("--data", self.data, "--cache", "cold", "--trials", 2,
                           "--output", self.out, "--", "sh", "-c", self.command)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(self.calls.read_text().splitlines()), 2)
        self.assertEqual({row["cache"] for row in read_rows(self.out)}, {"cold"})

    def test_peak_rss_sums_concurrent_processes(self):
        allocate = f"{sys.executable} -c 'b = bytearray(200 << 20); import time; time.sleep(1)'"
        result = run_timer("--data", self.data, "--cache", "warm", "--trials", 1,
                           "--output", self.out, "--", "sh", "-c",
                           f"{allocate} & {allocate} & wait")
        self.assertEqual(result.returncode, 0, result.stderr)
        # Two 200 MiB children alive together: the tree sum exceeds either one.
        self.assertGreater(float(read_rows(self.out)[0]["max_rss"]), 380)

    @unittest.skipIf(hasattr(os, "posix_fadvise"), "only platforms without posix_fadvise")
    def test_cold_refuses_when_eviction_is_unavailable(self):
        result = run_timer("--data", self.data, "--cache", "cold", "--trials", 1,
                           "--output", self.out, "--", "sh", "-c", self.command)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("posix_fadvise", result.stderr)
        self.assertFalse(self.calls.exists())

    def test_before_runs_untimed_ahead_of_every_trial(self):
        marker = self.dir / "before"
        result = run_timer("--data", self.data, "--cache", "warm", "--trials", 2,
                           "--before", f"echo cleanup >> {marker}; sleep 0.5",
                           "--output", self.out, "--", "sh", "-c", self.command)
        self.assertEqual(result.returncode, 0, result.stderr)
        # once before the warm-up and once before each timed trial
        self.assertEqual(len(marker.read_text().splitlines()), 3)
        for row in read_rows(self.out):
            self.assertLess(float(row["s"]), 0.5)

    def test_failing_command_fails_without_output(self):
        result = run_timer("--data", self.data, "--cache", "warm", "--trials", 2,
                           "--output", self.out, "--", "sh", "-c", "exit 3")
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.out.exists())


if __name__ == "__main__":
    unittest.main()
