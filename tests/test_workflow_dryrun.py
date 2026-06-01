import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


def write_fixture_config(root, *, data_as_list=False):
    data = root / "data"
    sample = data / "Sample_067"
    sample.mkdir(parents=True)
    (sample / "Sample_067_R1.fastq").write_text("@r1\nACGT\n+\n!!!!\n")
    (sample / "Sample_067_R2.fastq").write_text("@r2\nTGCA\n+\n!!!!\n")

    output = root / "out"
    checkm = root / "checkm"
    (checkm / "hmms").mkdir(parents=True)
    (checkm / "hmms" / "checkm.hmm").write_text("fake hmm\n")

    data_yaml = f"[{data}]" if data_as_list else str(data)
    config = root / "config.yaml"
    config.write_text(
        textwrap.dedent(
            f"""
            threads: 1
            task_memory: 1
            Percent_memory: 0.5
            execution_directory: {output}
            data: {data_yaml}
            filtering: ""
            assembly:
              parameters: ""
              groups:
                Elv:
                  - "Sample_067"
            binning:
              ssa_unique_sample: true
              cobinning_samples:
                - "Sample_067"
              concoct:
                execution: 1
                contig_size: 1000
                max_bin_nb: 2000
              metabat2:
                execution: 0
                contig_size: 1500
            annotation:
              checkm: {checkm}
              diamond: {{}}
              cat_db: ""
              cat_path: ""
              ip_db: ""
              kraken_db: ""
              kofamscan:
                profiles: ""
                ko_list: ""
              virsorter: ""
              plasmidnet_install: ""
              genomad_db: ""
            profile_duplicate: false
            metaspades_ssa: false
            slurm_partitions:
              "":
                name: ""
                min_mem: ""
                max_mem: ""
                min_threads: ""
                max_threads: ""
            """
        ).lstrip()
    )
    return config


@unittest.skipUnless(shutil.which("snakemake"), "snakemake is not installed")
class WorkflowDryRunTests(unittest.TestCase):
    def run_metahood(self, step, config, cores="1"):
        repo = Path(__file__).resolve().parents[1]
        return subprocess.run(
            [sys.executable, str(repo / "Metahood.py"), step, str(config), "--cores", cores, "--dryrun"],
            cwd=repo,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
            text=True,
            capture_output=True,
            timeout=120,
        )

    def test_sample_qc_dryrun_accepts_scalar_data_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = write_fixture_config(Path(tmp), data_as_list=False)

            result = self.run_metahood("sample_qc", config)

            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)

    def test_master_dryrun_builds_minimal_dag_with_one_core(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = write_fixture_config(Path(tmp), data_as_list=True)

            result = self.run_metahood("all", config)

            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)


if __name__ == "__main__":
    unittest.main()
