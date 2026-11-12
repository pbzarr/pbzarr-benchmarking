"""Cohort benchmark aggregation: pixi run python -m unittest tests.test_aggregate_benchmarks -v."""

import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
AGGREGATE = ROOT / "workflows/create/aggregate_benchmarks.py"


class AggregateBenchmarksTest(unittest.TestCase):
    def test_scan_trials_split_rule_and_cache_mode(self):
        with tempfile.TemporaryDirectory() as workspace:
            results = Path(workspace)
            benchmarks = results / "ds/pbz/t4/n2/r1/benchmarks"
            benchmarks.mkdir(parents=True)
            for name in ("region_scan.cold", "region_scan.warm", "over10.warm"):
                cache = name.split(".")[1]
                (benchmarks / f"{name}.tsv").write_text(
                    "trial\tcache\ts\tmax_rss\tio_in\tcpu_time\tmean_load\n"
                    f"2\t{cache}\t1.0\t10\t0\t2\t200\n"
                    f"1\t{cache}\t1.5\t10\t0\t3\t200\n"
                )
            result = subprocess.run([sys.executable, str(AGGREGATE), str(results)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (results / "benchmarks.tsv").open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields = reader.fieldnames
                rows = list(reader)
        # The timer's own trial and cache columns become metadata, not metrics.
        self.assertNotIn("benchmark_cache", fields)
        self.assertNotIn("trial", fields)
        self.assertEqual(
            sorted((row["rule"], row["cache"], row["measurement"], row["s"]) for row in rows),
            [("over10", "warm", "1", "1.5"), ("over10", "warm", "2", "1.0"),
             ("region_scan", "cold", "1", "1.5"), ("region_scan", "cold", "2", "1.0"),
             ("region_scan", "warm", "1", "1.5"), ("region_scan", "warm", "2", "1.0")],
        )


if __name__ == "__main__":
    unittest.main()
