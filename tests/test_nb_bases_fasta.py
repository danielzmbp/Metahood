import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class NbBasesFastaTests(unittest.TestCase):
    def test_counts_multiple_fasta_inputs(self):
        repo = Path(__file__).resolve().parents[1]
        script = repo / "scripts" / "nb_bases_fasta.py"

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample = root / "Sample_001"
            sample.mkdir()
            r1 = sample / "Sample_001_R1.fa"
            r2 = sample / "Sample_001_R2.fa"
            out = root / "nucleotides.tsv"
            r1.write_text(">r1\nACGT\n")
            r2.write_text(">r2\nAC\n")

            result = subprocess.run(
                [sys.executable, str(script), "-i", str(r1), str(r2), "-o", str(out)],
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                out.read_text(),
                "Normalisation\tSample_001\nNb_nucleotides\t6\n",
            )


if __name__ == "__main__":
    unittest.main()
