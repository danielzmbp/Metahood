import os
import shutil
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path

from test_workflow_dryrun import write_fixture_config


def write_executable(path, content):
    path.write_text(content)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


@unittest.skipUnless(shutil.which("snakemake"), "snakemake is not installed")
class CoverageRuleTests(unittest.TestCase):
    def test_bedtools_scans_bam_once_for_all_feature_types(self):
        repo = Path(__file__).resolve().parents[1]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = write_fixture_config(root, data_as_list=True)
            out = root / "out"
            annotation = out / "Elv" / "annotation"
            mapping = out / "Elv" / "map"
            annotation.mkdir(parents=True)
            mapping.mkdir(parents=True)
            (mapping / "Sample_067_mapped_sorted.bam").write_text("fake bam\n")
            (annotation / "contigs_bedtools_target_definition.tsv").write_text("contig1\t100\n")
            (annotation / "contigs.bed").write_text("contig1\t0\t100\tcontig1\n")
            (annotation / "orf.bed").write_text("contig1\t10\t20\torf1\n")
            (annotation / "contigs_C10K.bed").write_text("contig1\t0\t50\tcontig1.0\n")
            (annotation / "split_semibin2.bed").write_text("contig1\t0\t50\tcontig1_1\n")

            fakebin = root / "bin"
            fakebin.mkdir()
            count_file = root / "bedtools.count"
            write_executable(
                fakebin / "bedtools",
                """#!/bin/sh
if [ "$1" != "coverage" ]; then
    echo "unexpected bedtools command: $*" >&2
    exit 1
fi
while [ "$#" -gt 0 ]; do
    if [ "$1" = "-a" ]; then
        shift
        bed="$1"
        break
    fi
    shift
done
count_file="$BEDTOOLS_COUNT_FILE"
count=$(cat "$count_file" 2>/dev/null || printf 0)
printf '%s\n' $((count + 1)) > "$count_file"
while IFS= read -r line; do
    printf '%s\t1.0\n' "$line"
done < "$bed"
""",
            )

            target = "Elv/map/Sample_067.contigs.cov"
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
                    "--allowed-rules", "bedtools",
                    "--latency-wait", "1",
                    target,
                ],
                cwd=repo,
                env={
                    **os.environ,
                    "PATH": f"{fakebin}{os.pathsep}{os.environ['PATH']}",
                    "BEDTOOLS_COUNT_FILE": str(count_file),
                    "PYTHONDONTWRITEBYTECODE": "1",
                },
                text=True,
                capture_output=True,
                timeout=120,
            )

            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertEqual(count_file.read_text(), "1\n")
            self.assertEqual((mapping / "Sample_067.contigs.cov").read_text(), "contig1\t0\t100\tcontig1\t1.0\n")
            self.assertEqual((mapping / "Sample_067.orf.cov").read_text(), "contig1\t10\t20\torf1\t1.0\n")
            self.assertEqual((mapping / "Sample_067.contigs_C10K.cov").read_text(), "contig1\t0\t50\tcontig1.0\t1.0\n")
            self.assertEqual((mapping / "Sample_067.split_semibin2.cov").read_text(), "contig1\t0\t50\tcontig1_1\t1.0\n")


if __name__ == "__main__":
    unittest.main()
