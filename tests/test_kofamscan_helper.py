import unittest

from scripts import kofamscan_helper


class KofamscanHelperTests(unittest.TestCase):
    def test_intervals_overlap_matches_half_open_range_behavior(self):
        self.assertTrue(kofamscan_helper.intervals_overlap((1, 10), (9, 20)))
        self.assertFalse(kofamscan_helper.intervals_overlap((1, 10), (10, 20)))
        self.assertFalse(kofamscan_helper.intervals_overlap((1, 10), (11, 20)))


if __name__ == "__main__":
    unittest.main()
