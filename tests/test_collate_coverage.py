import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class CollateCoverageTests(unittest.TestCase):
    def test_collate_coverage_preserves_sorted_sample_output(self):
        repo = Path(__file__).resolve().parents[1]
        script = repo / "scripts" / "collate_coverage.py"

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample_b = root / "Sample_B.contigs.cov"
            sample_a = root / "Sample_A.contigs.cov"
            output = root / "coverage.tsv"

            sample_b.write_text("ctg1\t0\t10\tgene1\t2.0000000\nctg1\t10\t20\tgene2\t4.5\n")
            sample_a.write_text("ctg1\t0\t10\tgene1\t1.0000000\nctg1\t10\t20\tgene2\t3\n")

            result = subprocess.run(
                [sys.executable, str(script), "-o", str(output), "-l", str(sample_b), str(sample_a), "--batch-size", "1"],
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
                text=True,
                capture_output=True,
                timeout=120,
            )

            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertEqual(output.read_text(), "cov\tSample_A\tSample_B\ngene1\t1.0\t2.0\ngene2\t3.0\t4.5\n")

    def test_collate_coverage_rejects_mismatched_features(self):
        repo = Path(__file__).resolve().parents[1]
        script = repo / "scripts" / "collate_coverage.py"

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample_a = root / "Sample_A.contigs.cov"
            sample_b = root / "Sample_B.contigs.cov"
            output = root / "coverage.tsv"

            sample_a.write_text("ctg1\t0\t10\tgene1\t1\n")
            sample_b.write_text("ctg1\t0\t10\tgene2\t2\n")

            result = subprocess.run(
                [sys.executable, str(script), "-o", str(output), "-l", str(sample_a), str(sample_b), "--batch-size", "1"],
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
                text=True,
                capture_output=True,
                timeout=120,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("expected 'gene1'", result.stderr)


if __name__ == "__main__":
    unittest.main()
