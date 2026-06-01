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
class MappingRuleTests(unittest.TestCase):
    def test_bwa_mem_produces_bam_without_intermediate_sam(self):
        repo = Path(__file__).resolve().parents[1]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = write_fixture_config(root, data_as_list=True)
            out = root / "out"
            contigs = out / "Elv" / "contigs"
            contigs.mkdir(parents=True)
            (contigs / "contigs.fa").write_text(">contig1\nACGT\n")
            (contigs / "contigs.fa.fai").write_text("contig1\t4\t9\t4\t5\n")
            (contigs / "index.done").write_text("")

            sample = root / "data" / "Sample_067"
            (sample / "Sample_067_R1_trimmed.fastq").write_text("@r1\nACGT\n+\n!!!!\n")
            (sample / "Sample_067_R2_trimmed.fastq").write_text("@r2\nTGCA\n+\n!!!!\n")

            fakebin = root / "bin"
            fakebin.mkdir()
            write_executable(
                fakebin / "bwa",
                """#!/bin/sh
if [ "$1" = "mem" ]; then
    printf '@HD\tVN:1.0\n'
    printf 'read1\t0\tcontig1\t1\t60\t4M\t*\t0\t0\tACGT\t!!!!\n'
else
    echo "unexpected bwa command: $*" >&2
    exit 1
fi
""",
            )
            write_executable(
                fakebin / "samtools",
                """#!/bin/sh
cmd="$1"
shift
case "$cmd" in
    view|sort)
        cat
        ;;
    *)
        echo "unexpected samtools command: $cmd $*" >&2
        exit 1
        ;;
esac
""",
            )

            target = "Elv/map/Sample_067_mapped_sorted.bam"
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
                    "--allowed-rules", "bwa_mem",
                    "--latency-wait", "1",
                    target,
                ],
                cwd=repo,
                env={**os.environ, "PATH": f"{fakebin}{os.pathsep}{os.environ['PATH']}", "PYTHONDONTWRITEBYTECODE": "1"},
                text=True,
                capture_output=True,
                timeout=120,
            )

            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            self.assertTrue((out / target).exists())
            self.assertFalse((out / "Elv" / "map" / "Sample_067.sam").exists())


if __name__ == "__main__":
    unittest.main()
