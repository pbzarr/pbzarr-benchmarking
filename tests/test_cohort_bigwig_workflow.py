"""Check cohort Megadepth creation and parallel exact pyBigWig scans."""

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import pyBigWig


ROOT = Path(__file__).resolve().parents[1]
# Cold eviction needs posix_fadvise; elsewhere only warm scans can run.
CACHES = ("cold", "warm") if hasattr(os, "posix_fadvise") else ("warm",)
TRIALS = 2


class CohortBigWigWorkflowTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = tempfile.TemporaryDirectory(prefix="cohort bigwig workflow ")
        cls.addClassCleanup(cls.workspace.cleanup)
        cls.root = Path(cls.workspace.name).resolve()
        workflow_dir = cls.root / "workflows/create"
        workflow_dir.mkdir(parents=True)
        for name in ("analysis.smk", "bigwig.smk", "samples.smk", "bigwig_scan.py", "timed_trials.py"):
            shutil.copy(ROOT / "workflows/create" / name, workflow_dir / name)
        regions = cls.root / "regions with spaces.bed"
        regions.write_text("ref1\t1\t9\nrefEmpty\t1\t5\n")
        cls.regions = regions
        sizes = cls.root / "chrom.sizes"
        sizes.write_text("ref1\t12\nrefEmpty\t6\n")
        cls.trial = cls.root / "results/bigwig/t2/n2/r1"
        cls.trials = [cls.trial.parent / f"r{replicate}" for replicate in (1, 2, 3)]
        data = cls.trial / "data"
        data.mkdir(parents=True)
        for sample, intervals in {
            "first": [(0, 3, 2.0), (7, 12, 4.0)],
            "second": [(0, 12, 3.0)],
        }.items():
            with pyBigWig.open(str(data / f"{sample}.bw"), "w") as bw:
                bw.addHeader([("ref1", 12), ("refEmpty", 6)])
                bw.addEntries(
                    ["ref1"] * len(intervals),
                    [start for start, _, _ in intervals],
                    ends=[end for _, end, _ in intervals],
                    values=[value for _, _, value in intervals],
                )
        cls.snakefile = workflow_dir / "Snakefile"
        cls.snakefile.write_text(
            "from pathlib import Path\n"
            f"ROOT = Path({str(cls.root)!r})\n"
            "WORKFLOW_DIR = ROOT / 'workflows/create'\n"
            "RESULTS = str(ROOT / 'results')\n"
            f"config = {{'regions_bed': {str(regions)!r}, 'analysis_trials': {TRIALS}}}\n"
            f"CHROM_SIZES = {str(sizes)!r}\n"
            "def selected_samples(n):\n"
            "    return [{'sample_id': name} for name in ('first', 'second')][:int(n)]\n"
            "rule all:\n"
            "    input:\n"
            f"        {[str(trial / 'benchmarks' / f'{phase}.{cache}.tsv') for trial in cls.trials for phase in ('full_scan', 'region_scan') for cache in CACHES]!r}\n"
            "include: 'analysis.smk'\n"
        )
        prefix = Path(os.environ.get("CONDA_PREFIX", ROOT / ".pixi/envs/default"))
        cls.env = dict(os.environ, PATH=str(prefix / "bin") + os.pathsep + os.environ["PATH"])

    def snakemake(self, *args, snakefile=None, directory=None):
        return subprocess.run(
            ["snakemake", "--snakefile", str(snakefile or self.snakefile),
             "--directory", str(directory or self.root),
             "--cores", "4", "--runtime-source-cache-path", str(self.root / "cache"), *args],
            env=self.env, text=True, capture_output=True, timeout=90,
        )

    def test_dag_uses_python_helper_and_two_workers_per_scan(self):
        result = self.snakemake("--dry-run", "--printshellcmds", "--forceall")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        text = result.stdout + result.stderr
        self.assertIn("bigwig_scan.py", text)
        jobs = 6 * len(CACHES)
        self.assertEqual(text.count("xargs -0 -n 1 -P 2"), jobs)
        self.assertEqual(text.count("threads: 2"), jobs)
        self.assertEqual(text.count(" --cache "), jobs)
        self.assertNotIn("rule full_scan_regions:", text)
        self.assertIn(str(self.regions), text)
        self.assertNotIn("megadepth", text)
        self.assertNotIn("build_bigwig_scan", text)
        self.assertNotIn("bigwig_scan.c", text)
        self.assertNotIn("-lBigWig", text)

    def test_creation_dag_uses_four_thread_megadepth_per_sample(self):
        fixtures = ROOT / "tests/fixtures/single_sample"
        manifest = self.root / "samples.tsv"
        manifest.write_text(
            "sample_id\talignment\tindex\n"
            + "".join(
                f"{sample}\t{fixtures / 'reads.bam'}\t{fixtures / 'reads.bam.bai'}\n"
                for sample in ("first", "second")
            )
        )
        config_file = self.root / "config.json"
        config_file.write_text(json.dumps({
            "dataset": "fixture-rna", "n_samples": 2,
            "manifest": str(manifest), "reference": str(fixtures / "ref.fa"),
            "chrom_sizes": str(self.root / "chrom.sizes"),
        }))
        result = self.snakemake(
            "--configfile", str(config_file), "--dry-run", "--printshellcmds",
            snakefile=self.root / "workflows/create/bigwig.smk",
            directory=self.root / "creation",
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        text = result.stdout + result.stderr
        self.assertEqual(text.count("--bigwig --filter-out 1796 --threads 4"), 2)
        self.assertEqual(text.count("--fasta"), 2)
        for sample in ("first", "second"):
            self.assertIn(f"mv data/{sample}.all.bw data/{sample}.bw", text)
            self.assertIn(f"benchmarks/megadepth/{sample}.tsv", text)
        self.assertNotIn("mosdepth", text)

    def test_scan_replicates_reuse_first_creation_and_write_separate_results_and_metrics(self):
        result = self.snakemake()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        for phase, bounds, expected in (
            ("full_scan", [("ref1", "0", "12"), ("refEmpty", "0", "6")],
             {"first": 26 / 12, "second": 3.0}),
            ("region_scan", [("ref1", "1", "9"), ("refEmpty", "1", "5")],
             {"first": 1.5, "second": 3.0}),
        ):
            for trial in self.trials:
                for cache in CACHES:
                    output = trial / f"{phase}.{cache}"
                    for sample, mean in expected.items():
                        with self.subTest(replicate=trial.name, phase=phase, cache=cache, sample=sample):
                            with (output / f"{sample}.bw.tsv").open() as handle:
                                rows = list(csv.reader(handle, delimiter="\t"))
                            self.assertEqual([tuple(row[:3]) for row in rows], bounds)
                            self.assertAlmostEqual(float(rows[0][3]), mean)
                            self.assertEqual(float(rows[1][3]), 0.0)
                    with (trial / "benchmarks" / f"{phase}.{cache}.tsv").open() as handle:
                        rows = list(csv.DictReader(handle, delimiter="\t"))
                    self.assertEqual([row["trial"] for row in rows], [str(t) for t in range(1, TRIALS + 1)])
                    self.assertEqual({row["cache"] for row in rows}, {cache})
                    self.assertGreaterEqual(float(rows[0]["s"]), 0)
                    self.assertEqual((trial / f"{phase}.{cache}.log").read_text(), "")
        for trial in self.trials[1:]:
            self.assertFalse((trial / "data").exists())


if __name__ == "__main__":
    unittest.main()
