import pathlib
import unittest

from primalbedtools.bedfiles import BedLine, BedLineParser, group_primer_pairs
from primalbedtools.primerpairs import PrimerPair


class TestPrimerPair(unittest.TestCase):
    # Read in basic bed file
    test_bedline = pathlib.Path("tests/test.bed")
    _test_headers, test_bedlines = BedLineParser.from_file(test_bedline)

    def test_group_primerpairs(self):
        # Test grouping of primer pairs
        primer_pairs = group_primer_pairs(self.test_bedlines)

        # Check correct number
        self.assertEqual(len(primer_pairs), 3)

    def test_primer_pair_creation(self):
        # Test creation of primer pairs
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "test_1_RIGHT_1", 1, "-", "ATGC")

        primer_pair = PrimerPair([fbedline], [rbedline])

        # Check correct attributes
        self.assertEqual(primer_pair.chrom, "chrom")
        self.assertEqual(primer_pair.pool, 1)
        self.assertEqual(primer_pair.amplicon_number, 1)
        self.assertEqual(primer_pair.prefix, "test")
        self.assertEqual(primer_pair.fbedlines, [fbedline])
        self.assertEqual(primer_pair.rbedlines, [rbedline])

    def test_primer_pair_creation_error_chromname(self):
        fbedline = BedLine("a", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "test_1_RIGHT_1", 1, "-", "ATGC")

        # Test error when chromname are different
        with self.assertRaises(ValueError):
            PrimerPair([fbedline], [rbedline])

    def test_primer_pair_creation_error_pool(self):
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "test_1_RIGHT_1", 2, "-", "ATGC")

        # Test error when pool are different
        with self.assertRaises(ValueError):
            PrimerPair([fbedline], [rbedline])

    def test_primer_pair_creation_error_amplicon_number(self):
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "test_2_RIGHT_1", 1, "-", "ATGC")

        # Test error when amplicon numbers are different
        with self.assertRaises(ValueError):
            PrimerPair([fbedline], [rbedline])

    def test_primer_pair_creation_error_no_forward_primers(self):
        rbedline = BedLine("chrom", 200, 220, "test_1_RIGHT_1", 1, "-", "ATGC")

        # Test error when no forward primers are present
        with self.assertRaises(ValueError):
            PrimerPair([], [rbedline])

    def test_primer_pair_creation_error_no_reverse_primers(self):
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")

        # Test error when no reverse primers are present
        with self.assertRaises(ValueError):
            PrimerPair([fbedline], [])


if __name__ == "__main__":
    unittest.main()
