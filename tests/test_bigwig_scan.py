"""Real BigWig checks: pixi run python -m unittest discover -s tests -v."""

from contextlib import nullcontext
import importlib.util
import io
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import pyBigWig


ROOT = Path(__file__).resolve().parents[1]
SCANNER = ROOT / "workflows/create/bigwig_scan.py"


class BigWigScanTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = tempfile.TemporaryDirectory(prefix="bigwig-scan-test-")
        cls.addClassCleanup(cls.workspace.cleanup)
        cls.directory = Path(cls.workspace.name)
        cls.bigwig = cls.directory / "fixture.bw"
        with pyBigWig.open(str(cls.bigwig), "w") as bw:
            bw.addHeader([("chrA", 4096), ("chrB", 1000), ("chrEmpty", 500)], maxZooms=0)
            bw.addEntries(
                ["chrA", "chrA", "chrA", "chrB", "chrB"],
                [0, 900, 3000, 100, 800],
                ends=[100, 2300, 4096, 500, 1000],
                values=[2.0, 3.0, 5.0, 4.0, 2.0],
            )
        cls.zoom_bigwig = cls.directory / "zoom.bw"
        with pyBigWig.open(str(cls.zoom_bigwig), "w") as bw:
            bw.addHeader([("chrZoom", 10000)])
            starts = list(range(0, 10000, 10))
            bw.addEntries(
                ["chrZoom"] * 1000, starts,
                ends=[start + 1 for start in starts],
                values=[float(i % 7 + 1) for i in range(1000)],
            )

    def scan(self, bed=None, bigwig=None):
        args = [sys.executable, str(SCANNER), str(bigwig or self.bigwig)]
        if bed is not None:
            bed_file = self.directory / "regions.bed"
            bed_file.write_text(bed)
            args.append(str(bed_file))
        return subprocess.run(args, capture_output=True, text=True)

    def rows(self, result):
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr, "")
        return [
            (chrom, int(start), int(end), float(mean))
            for chrom, start, end, mean in (
                line.split("\t") for line in result.stdout.splitlines()
            )
        ]

    def test_full_scan_reports_each_chromosome_with_zero_filled_gaps(self):
        self.assertEqual(self.rows(self.scan()), [
            ("chrA", 0, 4096, 2.412109375),
            ("chrB", 0, 1000, 2.0),
            ("chrEmpty", 0, 500, 0.0),
        ])

    def test_worker_caps_inherited_threads_before_loading_pybigwig_and_numpy(self):
        thread_variables = (
            "OMP_NUM_THREADS", "GOTO_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "OPENBLAS_DEFAULT_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
        )
        # Observe the real imports in a fresh process. This catches caps that
        # use setdefault or run after NumPy has initialized its thread pool.
        probe = """
import os
import runpy
import sys

variables = sys.argv[3:]
seen = []

class ImportObserver:
    def find_spec(self, fullname, path=None, target=None):
        if fullname in ('pyBigWig', 'numpy'):
            seen.append(fullname)
            for name in variables:
                assert os.environ.get(name) == '1', (fullname, name, os.environ.get(name))
        return None

sys.meta_path.insert(0, ImportObserver())
scanner, bigwig = sys.argv[1:3]
sys.argv = [scanner, bigwig]
try:
    runpy.run_path(scanner, run_name='__main__')
except SystemExit as error:
    assert error.code == 0, error.code
assert 'pyBigWig' in seen, seen
assert 'numpy' in seen or not sys.modules['pyBigWig'].numpy, seen
"""
        for inherited in ("64", None):
            with self.subTest(inherited=inherited):
                env = os.environ.copy()
                for name in thread_variables:
                    if inherited is None:
                        env.pop(name, None)
                    else:
                        env[name] = inherited
                result = subprocess.run(
                    [sys.executable, "-c", probe, str(SCANNER), str(self.bigwig), *thread_variables],
                    env=env, capture_output=True, text=True, timeout=30,
                )
                self.assertEqual(self.rows(result), [
                    ("chrA", 0, 4096, 2.412109375),
                    ("chrB", 0, 1000, 2.0),
                    ("chrEmpty", 0, 500, 0.0),
                ])

    def test_regions_clip_long_intervals_and_preserve_order_and_duplicates(self):
        self.assertEqual(self.rows(self.scan(
            "chrB\t0\t1000\nchrA\t950\t2250\n"
            "chrA\t50\t3050\nchrA\t50\t3050\n"
        )), [
            ("chrB", 0, 1000, 2.0),
            ("chrA", 950, 2250, 3.0),
            ("chrA", 50, 3050, 4550 / 3000),
            ("chrA", 50, 3050, 4550 / 3000),
        ])

    def test_zoom_summaries_do_not_change_full_or_region_results(self):
        with pyBigWig.open(str(self.zoom_bigwig)) as bw:
            self.assertGreater(bw.header()["nLevels"], 0)
            self.assertNotEqual(bw.stats("chrZoom", 3, 9997, type="sum")[0], 3996.0)
        # 142 cycles of 1..7 plus 1..6 sum to 3997. Starting at 3 omits the first 1.
        self.assertEqual(self.rows(self.scan(bigwig=self.zoom_bigwig)), [
            ("chrZoom", 0, 10000, 0.3997),
        ])
        self.assertEqual(self.rows(self.scan("chrZoom 3 9997\n", self.zoom_bigwig)), [
            ("chrZoom", 3, 9997, 3996 / 9994),
        ])

    def test_each_request_uses_one_exact_sum_over_the_whole_interval(self):
        spec = importlib.util.spec_from_file_location("bigwig_scan", SCANNER)
        scanner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(scanner)
        bed = self.directory / "spy-regions.bed"
        bed.write_text("chrA 50 3050\nchrB 0 1000\nchrA 50 3050\n")
        for bed_path, expected_calls in (
            (None, [
                mock.call("chrA", 0, 4096, type="sum", exact=True),
                mock.call("chrB", 0, 1000, type="sum", exact=True),
                mock.call("chrEmpty", 0, 500, type="sum", exact=True),
            ]),
            (bed, [
                mock.call("chrA", 50, 3050, type="sum", exact=True),
                mock.call("chrB", 0, 1000, type="sum", exact=True),
                mock.call("chrA", 50, 3050, type="sum", exact=True),
            ]),
        ):
            with self.subTest(bed=bed_path), pyBigWig.open(str(self.bigwig)) as real_bw:
                spy = mock.Mock(wraps=real_bw)
                with mock.patch.object(scanner.pyBigWig, "open", return_value=nullcontext(spy)):
                    output = io.StringIO()
                    scanner.scan(self.bigwig, bed_path, output=output)
                self.assertEqual(spy.stats.call_args_list, expected_calls)
                spy.values.assert_not_called()
                self.assertEqual(len(output.getvalue().splitlines()), 3)
                if bed_path is not None:
                    spy.intervals.assert_not_called()

    def test_missing_bases_are_zero_including_fully_uncovered_queries(self):
        self.assertEqual(self.rows(self.scan(
            "chrA 100 900\nchrB 0 200\nchrEmpty 25 100\n"
        )), [
            ("chrA", 100, 900, 0.0),
            ("chrB", 0, 200, 2.0),
            ("chrEmpty", 25, 100, 0.0),
        ])

    def test_empty_bigwig_reports_zero(self):
        empty = self.directory / "empty.bw"
        with pyBigWig.open(str(empty), "w") as bw:
            bw.addHeader([("chrEmpty", 500)])
        self.assertEqual(self.rows(self.scan(bigwig=empty)), [("chrEmpty", 0, 500, 0.0)])

    def test_bed_headers_blank_lines_comments_and_extra_columns(self):
        self.assertEqual(self.rows(self.scan(
            "\n  # comment\ntrack name=test\nbrowser position chrA:1-100\n"
            "chrA\t0\t100\tfeature\t0\t+\n\t\n"
        )), [("chrA", 0, 100, 2.0)])
        self.assertEqual(self.rows(self.scan("# no regions\n\n")), [])

    def test_invalid_bed_fails_with_line_number(self):
        for region in (
            "chrMissing 0 10", "chrA 0 4097", "chrA 10 10",
            "chrA 10 9", "chrA -1 10", "chrA zero 10",
            "chrA 0 1.5", "chrA 0 4294967296", "chrA 0", "chrA",
        ):
            with self.subTest(region=region):
                result = self.scan("# header\n" + region + "\n")
                self.assertGreater(result.returncode, 0, result.stderr)
                self.assertIn("line 2", result.stderr)

    def test_missing_bigwig_and_bed_fail(self):
        result = self.scan(bigwig=self.directory / "missing.bw")
        self.assertGreater(result.returncode, 0)
        self.assertTrue(result.stderr)
        result = subprocess.run(
            [sys.executable, SCANNER, self.bigwig, self.directory / "missing.bed"],
            capture_output=True, text=True,
        )
        self.assertGreater(result.returncode, 0)
        self.assertTrue(result.stderr)

    def test_corrupt_data_fails_instead_of_reporting_zero(self):
        data = bytearray(self.bigwig.read_bytes())
        data_offset = struct.unpack_from("<Q", data, 16)[0]
        data[data_offset + 8:data_offset + 12] = b"\0\0\0\0"
        corrupt = self.directory / "corrupt.bw"
        corrupt.write_bytes(data)
        for bed in (None, "chrA 0 100\n"):
            with self.subTest(bed=bed):
                result = self.scan(bed, bigwig=corrupt)
                self.assertGreater(result.returncode, 0, result.stdout)
                self.assertTrue(result.stderr)

    def test_corrupt_index_fails_instead_of_reporting_zero(self):
        data = bytearray(self.bigwig.read_bytes())
        index_offset = struct.unpack_from("<Q", data, 24)[0]
        data[index_offset:index_offset + 4] = b"\0\0\0\0"
        corrupt = self.directory / "corrupt-index.bw"
        corrupt.write_bytes(data)
        result = self.scan("chrA 0 100\n", bigwig=corrupt)
        self.assertGreater(result.returncode, 0, result.stdout)
        self.assertTrue(result.stderr)

    def test_nonfinite_signal_fails_instead_of_becoming_zero_or_nonfinite_output(self):
        for value in (float("nan"), float("inf"), float("-inf")):
            with self.subTest(value=value):
                nonfinite = self.directory / "nonfinite.bw"
                with pyBigWig.open(str(nonfinite), "w") as bw:
                    bw.addHeader([("chrA", 100)])
                    bw.addEntries(["chrA"], [0], ends=[100], values=[value])
                for bed in (None, "chrA 0 100\n"):
                    result = self.scan(bed, bigwig=nonfinite)
                    self.assertGreater(result.returncode, 0, result.stdout)
                    self.assertEqual(result.stdout, "")
                    self.assertTrue(result.stderr)

    def test_wrong_argument_count_fails_with_usage(self):
        for args in ([], [self.bigwig, "regions.bed", "extra"]):
            with self.subTest(args=args):
                result = subprocess.run(
                    [sys.executable, SCANNER, *args], capture_output=True, text=True,
                )
                self.assertGreater(result.returncode, 0)
                self.assertIn("usage:", result.stderr.lower())


if __name__ == "__main__":
    unittest.main()
