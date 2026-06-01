import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from test_workflow_dryrun import write_fixture_config


@unittest.skipUnless(shutil.which("snakemake"), "snakemake is not installed")
class AnnotationRuleTests(unittest.TestCase):
    def test_ip_best_hit_writes_final_orf_group(self):
        repo = Path(__file__).resolve().parents[1]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = write_fixture_config(root, data_as_list=True)
            out = root / "out"
            splits = out / "Elv" / "annotation" / "temp_splits"
            splits.mkdir(parents=True)
            for index in range(100):
                (splits / f"Batch_ip_{index}.tsv").write_text("")
            (splits / "Batch_ip_0.tsv").write_text(
                "orf1\tsigA\tmd5\t10\tIPR001\tPFAM\t1\t10\t5.0\textra\n"
                "orf1\tsigB\tmd5\t10\tIPR002\tPFAM\t1\t10\t1.0\textra\n"
                "orf2\tsigC\tmd5\t10\tIPR003\tPFAM\t1\t10\t2.0\textra\n"
            )

            target = "Elv/annotation/contigs_IP_best_hits.tsv"
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
                    "--allowed-rules", "IP_best_hit",
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
            lines = (out / target).read_text().splitlines()
            self.assertEqual(len(lines), 2)
            self.assertTrue(lines[0].startswith("orf1\tIPR002\t"))
            self.assertTrue(lines[1].startswith("orf2\tIPR003\t"))


if __name__ == "__main__":
    unittest.main()
