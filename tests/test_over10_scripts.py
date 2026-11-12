"""Per-base over_10 outputs agree with numpy: pixi run python -m unittest tests.test_over10_scripts -v."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pbzarr
import pyBigWig
import pyd4


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflows/create"
# The PBZ script imports the proposed-API shim as `pbzarr`.
PROPOSED_API = SCRIPTS / "proposed_api"
GENOME = {"chrA": 2_500_000, "chrB": 1_200_000, "chrC": 300_000}
SAMPLES = [f"S{i}" for i in range(5)]
# The 1kg benchmark layout: 1 Mi x 8 shards of 16 Ki x 8 transposed pcodec chunks.
CODECS = {
    "chunk_grid": {"name": "regular", "configuration": {"chunk_shape": [1048576, 8]}},
    "codecs": [{"name": "sharding_indexed", "configuration": {
        "chunk_shape": [16384, 8],
        "codecs": [
            {"name": "transpose", "configuration": {"order": [1, 0]}},
            {"name": "numcodecs.pcodec", "configuration": {"level": 8}},
        ],
        "index_codecs": [{"name": "bytes", "configuration": {"endian": "little"}}, {"name": "crc32c"}],
        "index_location": "end",
    }}],
}


class Over10ScriptsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        workspace = tempfile.TemporaryDirectory(prefix="over10-test-")
        cls.addClassCleanup(workspace.cleanup)
        cls.dir = Path(workspace.name)
        rng = np.random.default_rng(1)
        (cls.dir / "genome.txt").write_text("".join(f"{c}\t{n}\n" for c, n in GENOME.items()))
        depth = {}
        for sample in SAMPLES:
            depth[sample] = {}
            for chrom, length in GENOME.items():
                lengths = rng.integers(1, 5000, size=length // 100)
                values = rng.choice([0, 5, 9, 10, 11, 30], size=lengths.size,
                                    p=[.05, .05, .1, .1, .2, .5])
                depth[sample][chrom] = np.repeat(values, lengths)[:length].astype(np.int32)

        # samples with depth >= 10x at each base
        cls.covered = {c: sum(depth[s][c] >= 10 for s in SAMPLES) for c in GENOME}

        for sample in SAMPLES:
            records = []
            for chrom in GENOME:
                values = depth[sample][chrom]
                cuts = np.flatnonzero(np.diff(values)) + 1
                starts = np.concatenate(([0], cuts))
                ends = np.concatenate((cuts, [values.size]))
                records += [(chrom, int(a), int(b), int(values[a])) for a, b in zip(starts, ends)]
            bedgraph = cls.dir / f"{sample}.bedgraph"
            bedgraph.write_text("".join(f"{c}\t{a}\t{b}\t{v}\n" for c, a, b, v in records))
            with pyBigWig.open(str(cls.dir / f"{sample}.bw"), "w") as bw:
                bw.addHeader(list(GENOME.items()))
                nonzero = [r for r in records if r[3]]  # Megadepth writes no zero records
                bw.addEntries([r[0] for r in nonzero], [r[1] for r in nonzero],
                              ends=[r[2] for r in nonzero], values=[float(r[3]) for r in nonzero])
            subprocess.run(["d4tools", "create", "-g", str(cls.dir / "genome.txt"),
                            str(bedgraph), str(cls.dir / f"{sample}.d4")], check=True)
        subprocess.run(["d4tools", "merge", *[str(cls.dir / f"{s}.d4") for s in SAMPLES],
                        str(cls.dir / "cohort.d4")], check=True)
        pbzarr.create_store(cls.dir / "cohort.pbz")
        pbzarr.import_d4(cls.dir / "cohort.pbz", "depth",
                         [(str(cls.dir / f"{s}.d4"), s) for s in SAMPLES],
                         column_dim="sample", codecs=CODECS)

    def run_script(self, name, *args, env=None):
        result = subprocess.run([sys.executable, str(SCRIPTS / name), *map(str, args)],
                                capture_output=True, text=True, timeout=300, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_bigwig_writes_fraction(self):
        out = self.dir / "over_10.bw"
        self.run_script("over10_bigwig.py", out, 2, *[self.dir / f"{s}.bw" for s in SAMPLES])
        with pyBigWig.open(str(out)) as bw:
            self.assertEqual(bw.chroms(), GENOME)
            for chrom, length in GENOME.items():
                got = bw.values(chrom, 0, length, numpy=True)
                np.testing.assert_allclose(got, self.covered[chrom] / len(SAMPLES), rtol=1e-6)

    def test_d4_writes_sample_count(self):
        out = self.dir / "over_10.d4"
        self.run_script("over10_d4.py", self.dir / "cohort.d4", out, 2)
        d4 = pyd4.D4File(str(out))
        self.assertEqual(dict(d4.chroms()), GENOME)
        for chrom in GENOME:
            # pyd4 misreads a small fraction of input bases (known, ignored here)
            self.assertGreater((d4[chrom] == self.covered[chrom]).mean(), 0.99)

    def test_pbz_script_uses_proposed_api_names(self):
        source = (SCRIPTS / "over10_pbz.py").read_text()
        self.assertIn("import pbzarr as pbz", source)
        self.assertIn("from pbzarr.cluster import DaskClient", source)
        self.assertIn('ds["depth"]', source)
        self.assertIn("pbz.create(", source)

    def test_bigwig_and_d4_scripts_use_a_process_pool(self):
        for name in ("over10_bigwig.py", "over10_d4.py"):
            with self.subTest(script=name):
                source = (SCRIPTS / name).read_text()
                self.assertIn("from multiprocessing import Pool", source)
                self.assertNotIn("dask", source)

    def test_pbz_writes_fraction_track(self):
        out = self.dir / "over_10.pbz"
        env = dict(os.environ, PYTHONPATH=str(PROPOSED_API))
        self.run_script("over10_pbz.py", self.dir / "cohort.pbz", out, 2, env=env)
        # the installed pbzarr (format 0.4) must accept what the shim wrote
        track = pbzarr.open(out / "over_10")
        for chrom in GENOME:
            got = track.pbz.region(chrom).values
            np.testing.assert_allclose(got, self.covered[chrom] / len(SAMPLES), rtol=1e-6)


if __name__ == "__main__":
    unittest.main()
