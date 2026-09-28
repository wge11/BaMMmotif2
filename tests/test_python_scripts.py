"""Unit tests for the Python helper scripts in py/.

Run with:  python -m unittest discover -s tests
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(HERE, "..", "example")
sys.path.insert(0, os.path.join(HERE, "..", "py"))

import occur2bed  # noqa: E402

from utils import parse_meme, write_meme  # noqa: E402

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


class Occur2BedTest(unittest.TestCase):
    """Occurrences on both strands map back onto the right genomic bases."""

    def setUp(self):
        self.genome = "ACGTTGCAAGGCTTACCGATGACTCATTTGGCAACGTAGCTAGGCTAATCGGATCCATGACTCATAAC"
        self.region_start = 5
        self.region = self.genome[5:60]  # sequence scanned by BaMM!motif
        self.L = len(self.region)
        # the scanned string: forward strand, separator, reverse complement
        self.scanned = self.region + "N" + revcomp(self.region)

    def check(self, start0, width):
        pattern = self.scanned[start0 : start0 + width]
        strand = "+" if start0 < self.L else "-"
        chrom, start, end, strand, _ = occur2bed.occurrence_to_bed(
            f">chrT:{self.region_start}-{self.region_start + self.L}",
            self.L,
            strand,
            f"{start0 + 1}..{start0 + width}",
            "1e-5",
        )
        self.assertEqual(chrom, "chrT")
        genomic = self.genome[start:end]
        self.assertEqual(genomic if strand == "+" else revcomp(genomic), pattern)

    def test_forward_strand(self):
        for start0 in range(0, self.L - 8):
            self.check(start0, 8)

    def test_reverse_strand(self):
        for start0 in range(self.L + 1, 2 * self.L + 1 - 8):
            self.check(start0, 8)

    def test_header_without_coordinates(self):
        with self.assertRaises(ValueError):
            occur2bed.occurrence_to_bed(">chr1", 10, "+", "1..5", "0.1")


class MemeTest(unittest.TestCase):
    def test_parse_example(self):
        dataset = parse_meme(os.path.join(EXAMPLE, "PWM_peng10.meme"))
        self.assertEqual(dataset["alphabet"], "ACGT")
        self.assertAlmostEqual(dataset["bg_freq"][0], 0.274379)
        self.assertEqual(len(dataset["models"]), 6)
        model = dataset["models"][0]
        self.assertEqual(model["motif_length"], 12)
        self.assertEqual(len(model["pwm"]), 12)
        self.assertEqual(model["bg_freq"], dataset["bg_freq"])

    def test_round_trip(self):
        dataset = parse_meme(os.path.join(EXAMPLE, "PWM_peng10.meme"))
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.meme")
            write_meme(dataset, path)
            again = parse_meme(path)
        self.assertEqual(again["alphabet"], dataset["alphabet"])
        self.assertEqual(len(again["models"]), len(dataset["models"]))
        for a, b in zip(again["models"], dataset["models"]):
            self.assertEqual(a["model_id"], b["model_id"])
            for row_a, row_b in zip(a["pwm"], b["pwm"]):
                for x, y in zip(row_a, row_b):
                    self.assertAlmostEqual(x, y, places=4)


if __name__ == "__main__":
    unittest.main()
