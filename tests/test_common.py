import os
import sys
import tempfile
import unittest
from types import SimpleNamespace

sys.dont_write_bytecode = True

from scripts import common


class CommonHelperTests(unittest.TestCase):
    def test_scalar_paths_are_normalized_as_lists(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(common.expand_path_list(tmp), [os.path.realpath(tmp)])

    def test_default_slurm_partition_is_local(self):
        config = {}
        common.fill_default_values(config)

        self.assertFalse(common.has_slurm_partitions(config))
        self.assertEqual(common.get_slurm_partitions(config), [["", 0, 0, 0, 0]])
        self.assertEqual(config["nb_concurrent_map"], 4)

    def test_legacy_misspelled_map_key_is_supported(self):
        config = {"nb_concurent_map": 2}
        common.fill_default_values(config)

        self.assertEqual(config["nb_concurrent_map"], 2)

    def test_named_slurm_partition_is_parsed(self):
        config = {
            "slurm_partitions": {
                "compute": {
                    "name": "ei-compute",
                    "min_mem": "1",
                    "max_mem": "10",
                    "min_threads": "1",
                    "max_threads": "16",
                }
            }
        }

        self.assertTrue(common.has_slurm_partitions(config))
        self.assertEqual(
            common.get_slurm_partitions(config),
            [["ei-compute", 1000, 10000, 1, 16]],
        )

    def test_get_resource_real_handles_no_partitions(self):
        input_files = SimpleNamespace(size=0)

        self.assertEqual(
            common.get_resource_real(None, input_files, 4, 1, SLURM_PARTITIONS=[], mode="partition"),
            "",
        )
        self.assertEqual(
            common.get_resource_real(None, input_files, 4, 1, SLURM_PARTITIONS=[], mode="threads"),
            4,
        )

    def test_megahit_style_resource_scaling_uses_base2_floor(self):
        partitions = [["ei-compute", 1000, 600000, 1, 64]]

        self.assertEqual(
            common.get_resource_real(None, SimpleNamespace(size=0), 32, 1, SLURM_PARTITIONS=partitions, mode="mem", mult=8, min_size=8192),
            65536,
        )
        self.assertEqual(
            common.get_resource_real(None, SimpleNamespace(size=20_000_000_000), 32, 1, SLURM_PARTITIONS=partitions, mode="mem", mult=8, min_size=8192),
            160000,
        )
        self.assertEqual(
            common.get_resource_real(None, SimpleNamespace(size=100_000_000_000), 32, 1, SLURM_PARTITIONS=partitions, mode="mem", mult=8, min_size=8192),
            600000,
        )

    def test_replace_extensions_tracks_filter_and_fastq_trimming(self):
        read = "/data/Sample/Sample_R1.fastq.gz"

        self.assertEqual(
            common.replace_extensions(read, ""),
            "/data/Sample/Sample_R1_trimmed.fastq.gz",
        )
        self.assertEqual(
            common.replace_extensions(read, "/mask.fa"),
            "/data/Sample/Filtered_Sample_R1_trimmed.fastq.gz",
        )

    def test_get_fastp_json_tracks_filter_without_trimming(self):
        read = "/data/Sample/Sample_R1.fastq.gz"

        self.assertEqual(
            common.get_fastp_json(read, ""),
            "/data/Sample/Sample_.fastq.gz_fastp.json",
        )
        self.assertEqual(
            common.get_fastp_json(read, "/mask.fa"),
            "/data/Sample/Filtered_Sample_.fastq.gz_fastp.json",
        )


if __name__ == "__main__":
    unittest.main()
