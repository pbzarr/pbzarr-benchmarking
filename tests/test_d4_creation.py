"""Exercise dense and native sparse D4 cohort creation on BAM and CRAM.

Run with: pixi run python -m unittest discover -s tests -v
"""

import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests/fixtures/single_sample"


class D4CreationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = tempfile.TemporaryDirectory(prefix="d4-creation-test-")
        cls.addClassCleanup(cls.workspace.cleanup)
        cls.directory = Path(cls.workspace.name)
        prefix = Path(os.environ.get("CONDA_PREFIX", ROOT / ".pixi/envs/default"))
        cls.d4tools = prefix / "bin/d4tools"
        env = dict(os.environ, PATH=str(prefix / "bin") + os.pathsep + os.environ["PATH"])
        genome = cls.directory / "chrom.sizes"
        genome.write_text("ref1\t100\nref2\t60\n")
        cls.outputs = {}
        cls.dry_runs = {}
        cls.sparse_outputs = {}
        cls.sparse_dry_runs = {}
        cls.env = env
        cls.snakemake = prefix / "bin/snakemake"
        cls.base_configs = {}

        for alignment_format, index_suffix in (("bam", ".bai"), ("cram", ".crai")):
            alignment = FIXTURES / f"reads.{alignment_format}"
            manifest = cls.directory / f"samples-{alignment_format}.tsv"
            manifest.write_text(
                "sample_id\talignment\tindex\n"
                + "".join(
                    f"{sample}\t{alignment}\t{alignment}{index_suffix}\n"
                    for sample in ("first", "second")
                )
            )
            base_config = {
                "dataset": "fixture-wgs",
                "n_samples": 2,
                "manifest": str(manifest),
                "reference": str(FIXTURES / "ref.fa"),
                "chrom_sizes": str(genome),
            }
            cls.base_configs[alignment_format] = base_config
            variants = [
                ("dense", base_config, False, True),
                ("sparse", {**base_config, "dataset": "fixture-RnA"}, True, True),
            ]
            if alignment_format == "bam":
                variants.append((
                    "sparse-uncompressed",
                    {**base_config, "d4_sparse": True, "d4_deflate": False},
                    True, False,
                ))
            for name, config, sparse, compressed in variants:
                workdir = cls.directory / f"{alignment_format}-{name}"
                dry_result = cls.run_workflow(workdir, config, dry_run=True)
                result = cls.run_workflow(workdir, config)
                for invocation in (dry_result, result):
                    if invocation.returncode:
                        raise AssertionError(invocation.stdout + invocation.stderr)
                output = workdir / "data/cohort.d4"
                dry_text = dry_result.stdout + dry_result.stderr
                if sparse:
                    cls.sparse_outputs[alignment_format, compressed] = output
                    cls.sparse_dry_runs[alignment_format, compressed] = dry_text
                else:
                    cls.outputs[alignment_format] = output
                    cls.dry_runs[alignment_format] = dry_text

        cls.regions = cls.directory / "regions.bed"
        cls.regions.write_text(
            "ref1\t0\t100\nref1\t0\t10\nref1\t12\t22\n"
            "ref1\t65\t70\nref1\t70\t100\nref2\t0\t60\n"
        )

    @classmethod
    def run_workflow(cls, workdir, config, dry_run=False):
        workdir.mkdir(exist_ok=True)
        config_file = workdir / "config.json"
        config_file.write_text(json.dumps(config))
        command = [
            str(cls.snakemake),
            "--snakefile", str(ROOT / "workflows/create/d4.smk"),
            "--directory", str(workdir), "--cores", "8",
            "--configfile", str(config_file),
            "--runtime-source-cache-path", str(workdir / "source-cache"),
        ]
        if dry_run:
            command += ["--dry-run", "--printshellcmds", "--forceall"]
        return subprocess.run(
            command, env=cls.env, text=True, capture_output=True, timeout=60,
        )

    def run_d4(self, *args):
        result = subprocess.run(
            [str(self.d4tools), *map(str, args)],
            capture_output=True, timeout=20,
        )
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        return result.stdout

    def test_creation_dag_uses_four_thread_mosdepth_and_merges_sample_tracks(self):
        for alignment_format, text in self.dry_runs.items():
            with self.subTest(alignment_format=alignment_format):
                self.assertEqual(text.count("mosdepth -t 4 --fasta"), 2)
                self.assertEqual(text.count("--d4 tmp/"), 2)
                self.assertIn(
                    "d4tools merge tmp/first/coverage.per-base.d4:first "
                    "tmp/second/coverage.per-base.d4:second data/cohort.d4", text,
                )
                self.assertNotIn("encode_sparse_d4", text)
                self.assertNotIn("d4tools create", text)

    def test_rna_dag_uses_native_sparse_alignment_encoding_and_merge(self):
        for (alignment_format, compressed), text in self.sparse_dry_runs.items():
            with self.subTest(alignment_format=alignment_format, compressed=compressed):
                self.assertEqual(text.count("d4tools create -t 4"), 2)
                self.assertEqual(text.count("--dict-file"), 2)
                self.assertIn("d4-zero.dict", text)
                self.assertIn(f"-r {FIXTURES / 'ref.fa'}", text)
                self.assertIn(str(FIXTURES / f"reads.{alignment_format}"), text)
                self.assertIn(
                    "d4tools merge tmp/d4/first.d4:first tmp/d4/second.d4:second "
                    "data/cohort.d4", text,
                )
                self.assertEqual(text.count("d4tools create -t 4 -z"), 2 if compressed else 0)
                self.assertNotIn("mosdepth", text)

    def test_explicit_false_keeps_rna_on_the_dense_path(self):
        result = self.run_workflow(
            self.directory / "rna-dense-override",
            {**self.base_configs["bam"], "dataset": "fixture-RNA", "d4_sparse": False},
            dry_run=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        text = result.stdout + result.stderr
        self.assertEqual(text.count("mosdepth -t 4 --fasta"), 2)
        self.assertNotIn("d4tools create", text)

    def test_sparse_settings_require_booleans(self):
        for field in ("d4_sparse", "d4_deflate"):
            with self.subTest(field=field):
                result = self.run_workflow(
                    self.directory / f"invalid-{field}",
                    {**self.base_configs["bam"], field: "false"}, dry_run=True,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("boolean", result.stdout + result.stderr)

    def test_sparse_merged_tracks_have_zero_bit_primary_tables_and_requested_compression(self):
        for (alignment_format, compressed), output in self.sparse_outputs.items():
            for sample in ("first", "second"):
                with self.subTest(alignment_format=alignment_format, compressed=compressed, sample=sample):
                    metadata = json.loads(
                        self.run_d4("framedump", output, f"{sample}/.metadata").rstrip(b"\0")
                    )
                    self.assertEqual(metadata["dictionary"], {"SimpleRange": {"low": 0, "high": 1}})
                    self.assertEqual(self.run_d4("framedump", output, f"{sample}/.ptab"), b"")
                    secondary = json.loads(
                        self.run_d4("framedump", output, f"{sample}/.stab/.metadata").rstrip(b"\0")
                    )
                    self.assertEqual(
                        secondary["compression"], {"Deflate": 5} if compressed else "NoCompression",
                    )

    def test_sparse_native_coverage_remains_distinct_from_mosdepth(self):
        # Native d4tools counts this fixture's skipped/deleted spans differently;
        # pin that behavior rather than accidentally switching to mosdepth input.
        expected = (
            "ref1\t0\t10\t0\nref1\t10\t12\t2\nref1\t12\t22\t3\n"
            "ref1\t22\t25\t2\nref1\t25\t30\t3\nref1\t30\t40\t1\n"
            "ref1\t40\t45\t2\nref1\t45\t50\t1\nref1\t50\t60\t2\n"
            "ref1\t60\t65\t1\nref1\t65\t70\t0\nref1\t70\t98\t1\n"
            "ref1\t98\t100\t0\nref2\t0\t60\t0\n"
        )
        for key, output in self.sparse_outputs.items():
            with self.subTest(alignment_format=key[0], compressed=key[1]):
                actual = self.run_d4("view", "-H", output).decode().splitlines()
                self.assertEqual(actual[0].split("\t")[-2:], ["first", "second"])
                for column in (3, 4):
                    track = "".join(
                        "\t".join(row[:3] + [row[column]]) + "\n"
                        for row in (line.split("\t") for line in actual[1:])
                    )
                    self.assertEqual(track, expected)
                stats = self.run_d4("stat", "--no-index", "-s", "mean", "-t", "1", output).decode()
                self.assertEqual(stats.splitlines(), ["ref1\t0\t100\t1.33\t1.33", "ref2\t0\t60\t0\t0"])

    def test_merged_tracks_preserve_mosdepth_splice_and_deletion_semantics(self):
        expected = (
            "ref1\t0\t10\t0\nref1\t10\t12\t2\nref1\t12\t22\t3\n"
            "ref1\t22\t30\t2\nref1\t30\t40\t1\nref1\t40\t45\t2\n"
            "ref1\t45\t55\t1\nref1\t55\t60\t2\nref1\t60\t65\t1\n"
            "ref1\t65\t70\t0\nref1\t70\t75\t1\nref1\t75\t85\t0\n"
            "ref1\t85\t98\t1\nref1\t98\t100\t0\nref2\t0\t60\t0\n"
        )
        for alignment_format, output in self.outputs.items():
            with self.subTest(alignment_format=alignment_format):
                actual = self.run_d4("view", "-H", output).decode().splitlines()
                self.assertEqual(actual[0].split("\t")[-2:], ["first", "second"])
                for column in (3, 4):
                    track = "".join(
                        "\t".join(row[:3] + [row[column]]) + "\n"
                        for row in (line.split("\t") for line in actual[1:])
                    )
                    self.assertEqual(track, expected)

    def test_region_means_cover_both_merged_samples(self):
        expected = {
            ("ref1", 0, 100): 1.13,
            ("ref1", 0, 10): 0.0,
            ("ref1", 12, 22): 3.0,
            ("ref1", 65, 70): 0.0,
            ("ref1", 70, 100): 18 / 30,
            ("ref2", 0, 60): 0.0,
        }
        for alignment_format, output in self.outputs.items():
            with self.subTest(alignment_format=alignment_format):
                text = self.run_d4(
                    "stat", "--no-index", "-t", "1", "-s", "mean", "-H",
                    "-r", self.regions, output,
                ).decode().splitlines()
                self.assertEqual(text[0].split("\t")[-2:], ["first", "second"])
                actual = {
                    (row[0], int(row[1]), int(row[2])): list(map(float, row[3:]))
                    for row in (line.split("\t") for line in text[1:])
                }
                self.assertEqual(actual.keys(), expected.keys())
                for region, mean in expected.items():
                    self.assertEqual(len(actual[region]), 2)
                    for value in actual[region]:
                        self.assertAlmostEqual(value, mean, places=6)

    def test_creation_retains_metrics_and_removes_temporary_sample_data(self):
        variants = [(key, output, "mosdepth") for key, output in self.outputs.items()]
        variants += [(key, output, "encode_sparse_d4") for key, output in self.sparse_outputs.items()]
        for key, output, creation_rule in variants:
            workdir = output.parents[1]
            with self.subTest(variant=key):
                for sample in ("first", "second"):
                    self.assertFalse((workdir / "tmp" / sample).exists())
                    self.assertFalse((workdir / "tmp/d4" / f"{sample}.d4").exists())
                for metric in (f"{creation_rule}/first", f"{creation_rule}/second", "merge"):
                    with (workdir / "benchmarks" / f"{metric}.tsv").open() as handle:
                        rows = list(csv.DictReader(handle, delimiter="\t"))
                    self.assertEqual(len(rows), 1)
                    self.assertGreaterEqual(float(rows[0]["s"]), 0)


if __name__ == "__main__":
    unittest.main()
