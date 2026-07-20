import logging
import pathlib
import unittest

from primalbedtools.amplicons import Amplicon, create_amplicons, do_pp_ol
from primalbedtools.bedfiles import BedLine, BedLineParser, group_primer_pairs
from tests.loghelpers import capture_logs

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

        # Test error when chromname are different, and that the offending
        # primers are named per-chromosome in the message.
        with self.assertRaisesRegex(ValueError, r"chrom a: test_1_LEFT_1"):
            Amplicon([fbedline], [rbedline])
        with self.assertRaisesRegex(ValueError, r"chrom chrom: test_1_RIGHT_1"):
            Amplicon([fbedline], [rbedline])

    def test_primer_pair_creation_error_pool(self):
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "test_1_RIGHT_1", 2, "-", "ATGC")

        # Test error when pool are different, and that each pool's members
        # are named on their own line.
        with self.assertRaisesRegex(ValueError, r"pool 1: test_1_LEFT_1"):
            Amplicon([fbedline], [rbedline])
        with self.assertRaisesRegex(ValueError, r"pool 2: test_1_RIGHT_1"):
            Amplicon([fbedline], [rbedline])

    def test_primer_pair_creation_error_amplicon_number(self):
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "test_2_RIGHT_1", 1, "-", "ATGC")

        # Test error when amplicon numbers are different, and offenders are named.
        with self.assertRaisesRegex(ValueError, r"amplicon 1: test_1_LEFT_1"):
            Amplicon([fbedline], [rbedline])
        with self.assertRaisesRegex(ValueError, r"amplicon 2: test_2_RIGHT_1"):
            Amplicon([fbedline], [rbedline])

    def test_primer_pair_creation_warning_different_prefix(self):
        # Different prefixes only warn (they don't raise); the debug record should
        # name which primers carry each prefix, and both are kept in the name.
        fbedline = BedLine("chrom", 100, 120, "aScheme_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "bScheme_1_RIGHT_1", 1, "-", "ATGC")

        with self.assertLogs("primalbedtools.amplicons", level="DEBUG") as cm:
            amplicon = Amplicon([fbedline], [rbedline])
        output = "\n".join(cm.output)

        self.assertIn("prefix aScheme: aScheme_1_LEFT_1", output)
        self.assertIn("prefix bScheme: bScheme_1_RIGHT_1", output)
        self.assertEqual(amplicon.prefix, "aScheme-bScheme")
        self.assertEqual(amplicon.amplicon_name, "aScheme-bScheme_1")

    def test_matching_prefix_logs_nothing(self):
        # The single-prefix path must stay silent and keep the bare prefix.
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")
        rbedline = BedLine("chrom", 200, 220, "test_1_RIGHT_1", 1, "-", "ATGC")

        with capture_logs() as records:
            amplicon = Amplicon([fbedline], [rbedline])

        self.assertEqual(records, [])
        self.assertEqual(amplicon.prefix, "test")
        self.assertEqual(amplicon.prefixes, ["test"])

    def test_primer_pair_creation_error_no_forward_primers(self):
        rbedline = BedLine("chrom", 200, 220, "test_1_RIGHT_1", 1, "-", "ATGC")

        # Test error when no forward primers are present; the present reverse
        # primer should be named.
        with self.assertRaisesRegex(
            ValueError, r"Reverse primers present: test_1_RIGHT_1"
        ):
            Amplicon([], [rbedline])

    def test_primer_pair_creation_error_no_reverse_primers(self):
        fbedline = BedLine("chrom", 100, 120, "test_1_LEFT_1", 1, "+", "ATGC")

        # Test error when no reverse primers are present; the present forward
        # primer should be named.
        with self.assertRaisesRegex(
            ValueError, r"Forward primers present: test_1_LEFT_1"
        ):
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

    def test_do_pp_ol_half_open_boundary(self):
        # Amplicon 1: [0, 30)
        a1_left = BedLine("chrom", 0, 10, "test_1_LEFT_1", 1, "+", "ATGC")
        a1_right = BedLine("chrom", 20, 30, "test_1_RIGHT_1", 1, "-", "ATGC")
        amp1 = Amplicon([a1_left], [a1_right])

        # Amplicon 2 starts exactly at 30, so no overlap for half-open intervals
        a2_left = BedLine("chrom", 30, 40, "test_2_LEFT_1", 1, "+", "ATGC")
        a2_right = BedLine("chrom", 50, 60, "test_2_RIGHT_1", 1, "-", "ATGC")
        amp2 = Amplicon([a2_left], [a2_right])

        self.assertFalse(do_pp_ol(amp1, amp2))


class TestPrefixDivergenceReporting(unittest.TestCase):
    """create_amplicons summarises prefix divergence once for the whole scheme."""

    def _divergent_bedlines(self, count, left_prefix="aScheme", right_prefix="bScheme"):
        bedlines = []
        for n in range(1, count + 1):
            start = n * 1000
            bedlines.append(
                BedLine(
                    "chrom",
                    start,
                    start + 20,
                    f"{left_prefix}_{n}_LEFT_1",
                    1,
                    "+",
                    "ATGC",
                )
            )
            bedlines.append(
                BedLine(
                    "chrom",
                    start + 100,
                    start + 120,
                    f"{right_prefix}_{n}_RIGHT_1",
                    1,
                    "-",
                    "ATGC",
                )
            )
        return bedlines

    def test_warns_once_for_whole_scheme(self):
        with capture_logs(level=logging.WARNING) as records:
            create_amplicons(self._divergent_bedlines(5))

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].levelno, logging.WARNING)
        self.assertEqual(records[0].name, "primalbedtools.amplicons")

        message = records[0].getMessage()
        self.assertIn("5 amplicon(s)", message)
        self.assertIn("aScheme-bScheme_1", message)

    def test_large_scheme_is_truncated(self):
        with capture_logs(level=logging.WARNING) as records:
            create_amplicons(self._divergent_bedlines(15))

        message = records[0].getMessage()
        # Guards against regressing to one message per amplicon
        self.assertEqual(len(message.splitlines()), 2)
        self.assertIn("(+5 more)", message)

    def test_grouped_by_prefix_set(self):
        bedlines = self._divergent_bedlines(2, "aScheme", "bScheme")
        bedlines += self._divergent_bedlines(2, "cScheme", "dScheme")
        # Renumber the second pattern so amplicon numbers stay unique
        for bedline in bedlines[4:]:
            bedline.amplicon_number += 10

        with capture_logs(level=logging.WARNING) as records:
            create_amplicons(bedlines)

        message = records[0].getMessage()
        self.assertEqual(len(message.splitlines()), 3)
        self.assertIn("prefixes aScheme, bScheme", message)
        self.assertIn("prefixes cScheme, dScheme", message)

    def test_clean_bedfile_is_silent(self):
        _headers, bedlines = BedLineParser.from_file(TEST_BEDLINE)

        with capture_logs() as records:
            create_amplicons(bedlines)

        self.assertEqual(records, [])


if __name__ == "__main__":
    unittest.main()
