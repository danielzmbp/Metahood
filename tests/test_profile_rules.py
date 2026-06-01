import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from test_workflow_dryrun import write_fixture_config


@unittest.skipUnless(shutil.which("snakemake"), "snakemake is not installed")
class ProfileRuleTests(unittest.TestCase):
    def test_get_percent_uses_declared_inputs(self):
        repo = Path(__file__).resolve().parents[1]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = write_fixture_config(root, data_as_list=True)
            out = root / "out"
            map_dir = out / "Elv" / "map"
            map_dir.mkdir(parents=True)
            (map_dir / "Sample_067_mapped_read.txt").write_text("25\n")

            sample = root / "data" / "Sample_067"
            (sample / "Sample_067_.fastq_fastp.json").write_text('{"summary": {"after_filtering": {"total_reads": 100}}}')

            target = "Elv/profile/mapping_percent.tsv"
            result = subprocess.run(
                [
                    "snakemake",
                    "--directory", str(out),
                    "--snakefile", str(repo / "Master.snake"),
                    "-k",
                    "--config", f"LOCAL_DIR={repo}", f"CONFIG_PATH={config}", f"EXEC_DIR={out}",
                    "--configfile", str(config),
                    "--cores", "1",
                    "--resources", "memG=1", "nb_map=200",
                    "--allowed-rules", "get_percent",
                    "--latency-wait", "1",
                    target,
                ],
                cwd=repo,
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
                text=True,
                capture_output=True,
                timeout=120,
            )

            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertEqual((out / target).read_text(), "Sample_067\t0.25\n")


if __name__ == "__main__":
    unittest.main()
