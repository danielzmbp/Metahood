import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


@unittest.skipUnless(shutil.which("bedtools") and shutil.which("samtools"), "bedtools and samtools are required")
class CoverageEquivalenceTests(unittest.TestCase):
    def test_combined_coverage_matches_separate_bedtools_outputs(self):
        repo = Path(__file__).resolve().parents[1]
        script = repo / "scripts" / "combined_bedtools_coverage.py"

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genome = root / "genome.tsv"
            genome.write_text("contig1\t100\ncontig2\t50\n")

            beds = {
                "contigs": root / "contigs.bed",
                "orf": root / "orf.bed",
                "contigs_C10K": root / "contigs_C10K.bed",
                "split_semibin2": root / "split_semibin2.bed",
            }
            beds["contigs"].write_text("contig1\t0\t100\tcontig1\ncontig2\t0\t50\tcontig2\n")
            beds["orf"].write_text("contig1\t0\t10\torf1\ncontig1\t10\t20\torf2\ncontig2\t0\t10\torf3\n")
            beds["contigs_C10K"].write_text("contig1\t0\t50\tcontig1.0\ncontig1\t50\t100\tcontig1.1\n")
            beds["split_semibin2"].write_text("contig1\t0\t50\tcontig1_1\ncontig1\t50\t100\tcontig1_2\n")

            sam = root / "reads.sam"
            unsorted_bam = root / "reads.unsorted.bam"
            bam = root / "reads.bam"
            sam.write_text(
                "@HD\tVN:1.6\tSO:coordinate\n"
                "@SQ\tSN:contig1\tLN:100\n"
                "@SQ\tSN:contig2\tLN:50\n"
                "read1\t0\tcontig1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n"
                "read2\t0\tcontig1\t11\t60\t10M\t*\t0\t0\tCCCCCCCCCC\tIIIIIIIIII\n"
                "read3\t0\tcontig1\t21\t60\t10M\t*\t0\t0\tGGGGGGGGGG\tIIIIIIIIII\n"
                "read4\t0\tcontig2\t1\t60\t10M\t*\t0\t0\tTTTTTTTTTT\tIIIIIIIIII\n"
            )
            subprocess.run(["samtools", "view", "-bS", str(sam), "-o", str(unsorted_bam)], check=True, capture_output=True, text=True)
            subprocess.run(["samtools", "sort", "-o", str(bam), str(unsorted_bam)], check=True, capture_output=True, text=True)

            old_dir = root / "old"
            new_dir = root / "new"
            old_dir.mkdir()
            new_dir.mkdir()
            old_outputs = {name: old_dir / f"{name}.cov" for name in beds}
            new_outputs = {name: new_dir / f"{name}.cov" for name in beds}

            for name, bed in beds.items():
                with old_outputs[name].open("w") as out_handle:
                    subprocess.run(
                        ["bedtools", "coverage", "-a", str(bed), "-b", str(bam), "-g", str(genome), "-mean", "-sorted"],
                        check=True,
                        stdout=out_handle,
                        text=True,
                    )

            command = [sys.executable, str(script), "--bam", str(bam), "--genome", str(genome)]
            for name, bed in beds.items():
                command.extend(["--target", name, str(bed), str(new_outputs[name])])
            subprocess.run(command, check=True, capture_output=True, text=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})

            for name in beds:
                self.assertEqual(new_outputs[name].read_text(), old_outputs[name].read_text(), name)


if __name__ == "__main__":
    unittest.main()
