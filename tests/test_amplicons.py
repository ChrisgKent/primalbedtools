import pathlib
import unittest

from primalbedtools.amplicons import Amplicon, create_amplicons
from primalbedtools.bedfiles import BedLine, BedLineParser, group_primer_pairs

TEST_BEDLINE = pathlib.Path(__file__).parent / "inputs/test.bed"
TEST_PROBE_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.probe.bed"


class TestAmplicon(unittest.TestCase):
    def setUp(self) -> None:
        self.test_bedline = TEST_BEDLINE
        self._test_headers, self.test_bedlines = BedLineParser.from_file(
            self.test_bedline
        )
        return super().setUp()

    def test_group_Amplicons(self):
        # Test grouping of primer pairs
        primer_pairs = group_primer_pairs(self.test_bedlines)

        # Check correct number
        self.assertEqual(len(primer_pairs), 3)

    def test_primer_pair_creation(self):
        # Test creation of primer pairs
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "test_1_RIGHT_1", 1, "-", "ATGC")

        primer_pair = Amplicon([fbedline], [rbedline])

        # Check correct attributes
        self.assertEqual(primer_pair.chrom, "chrom")
        self.assertEqual(primer_pair.pool, 1)
        self.assertEqual(primer_pair.amplicon_number, 1)
        self.assertEqual(primer_pair.prefix, "test")
        self.assertEqual(primer_pair.left, [fbedline])
        self.assertEqual(primer_pair.right, [rbedline])

    def test_primer_pair_creation_error_chromname(self):
        fbedline = BedLine("a", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "test_1_RIGHT_1", 1, "-", "ATGC")

        # Test error when chromname are different
        with self.assertRaises(ValueError):
            Amplicon([fbedline], [rbedline])

    def test_primer_pair_creation_error_pool(self):
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "test_1_RIGHT_1", 2, "-", "ATGC")

        # Test error when pool are different
        with self.assertRaises(ValueError):
            Amplicon([fbedline], [rbedline])

    def test_primer_pair_creation_error_amplicon_number(self):
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "test_2_RIGHT_1", 1, "-", "ATGC")

        # Test error when amplicon numbers are different
        with self.assertRaises(ValueError):
            Amplicon([fbedline], [rbedline])

    def test_primer_pair_creation_error_no_forward_primers(self):
        rbedline = BedLine("chrom", 200, 220, "test_1_RIGHT_1", 1, "-", "ATGC")

        # Test error when no forward primers are present
        with self.assertRaises(ValueError):
            Amplicon([], [rbedline])

    def test_primer_pair_creation_error_no_reverse_primers(self):
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")

        # Test error when no reverse primers are present
        with self.assertRaises(ValueError):
            Amplicon([fbedline], [])

    def test_create_Amplicons(self):
        # Create list of Amplicons
        amp = create_amplicons(self.test_bedlines)

        # check right amount of amp
        self.assertEqual(len(amp), 3)

    def test_ipool(self):
        amps = create_amplicons(self.test_bedlines)

        ipools = [amp.ipool for amp in amps]
        self.assertEqual(ipools, [0, 1, 0])

    def test_is_circular(self):
        amps = create_amplicons(self.test_bedlines)

        # Check Amplicon is not circular
        self.assertFalse(amps[0].is_circular)

        # Change primer coords
        amps[0].left[0].end = amps[0].right[0].end + 100
        amps[0].left[0].start = amps[0].right[0].start + 100

        # Check is now circular
        self.assertTrue(amps[0].is_circular)

    def test_coverage_start(self):
        amp = create_amplicons(self.test_bedlines)[0]

        self.assertEqual(amp.coverage_start, 78)

    def test_coverage_end(self):
        amp = create_amplicons(self.test_bedlines)[0]
        self.assertEqual(amp.coverage_end, 419)

    def test_to_amplicon_str(self):
        amp = create_amplicons(self.test_bedlines)[0]

        exp_str = "MN908947.3	47	447	SARS-CoV-2_1	1"
        self.assertEqual(amp.to_amplicon_str(), exp_str)

    def test_to_primertrim_str(self):
        amp = create_amplicons(self.test_bedlines)[0]
        exp_str = "MN908947.3	78	419	SARS-CoV-2_1	1"
        self.assertEqual(amp.to_primertrim_str(), exp_str)

    def test_get_regions(self):
        _headers, bedlines = BedLineParser.from_file(TEST_PROBE_BEDFILE)
        amp = create_amplicons(bedlines)[0]

        self.assertEqual(amp.left_region, (2010, 2030))
        self.assertEqual(amp.probe_region, (2035, 2060))
        self.assertEqual(amp.right_region, (2903, 2923))


if __name__ == "__main__":
    unittest.main()
