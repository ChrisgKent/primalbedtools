import os
import random
import tempfile
import unittest

from primalbedtools.bedfiles import BedLine
from primalbedtools.scheme import DEFAULT_CSV_HEADERS, Scheme
from tests.infiles import (
    TEST_ATTRIBUTES_BEDFILE,
    TEST_BEDFILE,
    TEST_PROBE_BEDFILE,
)


class TestScheme(unittest.TestCase):
    def setUp(self) -> None:
        self.maxDiff = None
        return super().setUp()

    def test_read_from_file(self):
        """
        Test round trip io for files
        """
        scheme = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        file_str = TEST_ATTRIBUTES_BEDFILE.read_text()
        self.assertEqual(file_str, scheme.to_str())

    def test_read_from_str(self):
        """
        Test round trip io for str
        """
        file_str = TEST_ATTRIBUTES_BEDFILE.read_text()
        scheme = Scheme.from_str(file_str)
        self.assertEqual(file_str, scheme.to_str())

    def test_sort(self):
        """
        Tests sorting returns to known order
        """
        scheme = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        file_str = TEST_ATTRIBUTES_BEDFILE.read_text()
        # Shuffle the bedlines
        random.seed(42)  # Set seed for reproducible results
        random.shuffle(scheme.bedlines)

        # Check string is dif
        self.assertNotEqual(scheme.to_str(), file_str)

        # Check sorting returns to original
        scheme.sort_bedlines()
        self.assertEqual(file_str, scheme.to_str())

    def test_read_from_file_probe(self):
        """
        Test round trip io for files with probe
        """
        scheme = Scheme.from_file(str(TEST_PROBE_BEDFILE))
        file_str = TEST_PROBE_BEDFILE.read_text()

        self.assertEqual(file_str, scheme.to_str())

    def test_read_from_str_probe(self):
        """
        Test round trip io for str with probe
        """
        file_str = TEST_PROBE_BEDFILE.read_text()
        scheme = Scheme.from_str(file_str)
        self.assertEqual(file_str, scheme.to_str())

    def test_sort_probe(self):
        """
        Tests sorting returns to known order with probe
        """
        scheme = Scheme.from_file(str(TEST_PROBE_BEDFILE))
        file_str = TEST_PROBE_BEDFILE.read_text()
        # Shuffle the bedlines
        random.seed(42)  # Set seed for reproducible results
        random.shuffle(scheme.bedlines)

        # Check string is dif
        self.assertNotEqual(scheme.to_str(), file_str)

        # Check sorting returns to original
        scheme.sort_bedlines()
        self.assertEqual(file_str, scheme.to_str())

    def test_parse_headers(self):
        """
        Check the header can be parsed as expected
        """
        scheme = Scheme.from_file(str(TEST_PROBE_BEDFILE))

        attr_dict = scheme.header_dict

        self.assertDictEqual(
            attr_dict,
            {
                "/3BHQ_1/": "BlackHoleQuencher1",
                "/56-FAM/": "FAM",
                "/5HEX/": "HEX",
                "example multiplexed-qPCR assay": None,
            },
        )

    def test_contains_probes(self):
        scheme = Scheme.from_file(str(TEST_PROBE_BEDFILE))
        self.assertTrue(scheme.contains_probes)

        scheme = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        self.assertFalse(scheme.contains_probes)

    def test_to_csv(self):
        # Read the scheme in
        scheme = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))

        # Check include_headers=True, use_header_aliases=False
        csv_str = scheme.to_delim_str(include_headers=True, use_header_aliases=False)
        csv_line_list = csv_str.splitlines()
        # Check default headers are present
        test_headers = csv_line_list[0].split(",")
        for exp_header in DEFAULT_CSV_HEADERS:
            self.assertIn(exp_header, test_headers, f"{exp_header} not in first line")

        # Check attribute headers are there with no aliases
        self.assertIn("pw", test_headers, "pw not in first line")
        self.assertIn("gc", test_headers, "gc not in first line")

        # Check all bedlines are present
        self.assertEqual(len(scheme.bedlines) + 1, len(csv_line_list))

    def test_to_csv_aliases(self):
        # Read the scheme in
        scheme = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        # Check include_headers=True, use_header_aliases=True
        csv_str = scheme.to_delim_str(include_headers=True, use_header_aliases=True)
        csv_line_list = csv_str.splitlines()
        # Check default headers are present
        test_headers = csv_line_list[0].split(",")
        for exp_header in DEFAULT_CSV_HEADERS:
            self.assertIn(exp_header, test_headers, f"{exp_header} not in first line")

        # Check attribute headers are there with no aliases
        self.assertIn("pw", test_headers, "pw not in first line")
        self.assertIn("fractiongc", test_headers, "fractiongc not in first line")

        # Check all bedlines are present
        self.assertEqual(len(scheme.bedlines) + 1, len(csv_line_list))

    def test_to_csv_no_header(self):
        # Read the scheme in
        scheme = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        # Check include_headers=True, use_header_aliases=True
        csv_str = scheme.to_delim_str(include_headers=False, use_header_aliases=True)
        csv_line_list = csv_str.splitlines()
        # Check default headers are present
        test_headers = csv_line_list[0].split(",")
        for exp_header in DEFAULT_CSV_HEADERS:
            self.assertNotIn(
                exp_header, test_headers, f"{exp_header} found in first line"
            )
        # Check all bedlines are present
        self.assertEqual(len(scheme.bedlines), len(csv_line_list))

    def test_to_csv_does_not_mutate_default_headers(self):
        scheme = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        original_headers = list(DEFAULT_CSV_HEADERS)
        _csv_str = scheme.to_delim_str(include_headers=True, use_header_aliases=False)
        self.assertEqual(DEFAULT_CSV_HEADERS, original_headers)

    def test_to_csv_no_trailing_newline(self):
        scheme = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        self.assertFalse(scheme.to_delim_str().endswith("\n"))

    def test_to_csv_attribute_colliding_with_fixed_column_raises(self):
        bedline = BedLine("chr1", 100, 120, "test_1_LEFT_1", 1, "+", "ACGT")
        bedline.attributes = {"chrom": "elsewhere"}
        scheme = Scheme(headers=[], bedlines=[bedline])

        with self.assertRaisesRegex(ValueError, r"Attribute \(chrom\) collides"):
            scheme.to_delim_str()


class TestFromDelimStr(unittest.TestCase):
    """from_delim_str reverses to_delim_str."""

    def round_trip(self, bedfile):
        scheme = Scheme.from_file(str(bedfile))
        parsed = Scheme.from_delim_str(scheme.to_delim_str())
        return scheme, parsed

    def test_round_trips_bedlines(self):
        for bedfile in (TEST_ATTRIBUTES_BEDFILE, TEST_PROBE_BEDFILE, TEST_BEDFILE):
            with self.subTest(bedfile=bedfile.name):
                scheme, parsed = self.round_trip(bedfile)
                self.assertEqual(
                    [bl.to_bed() for bl in parsed.bedlines],
                    [bl.to_bed() for bl in scheme.bedlines],
                )

    def test_round_trips_attributes(self):
        scheme, parsed = self.round_trip(TEST_ATTRIBUTES_BEDFILE)

        # pw is coerced to float by the attributes setter, so the types survive
        self.assertEqual(parsed.bedlines[0].attributes, {"pw": 1.4, "gc": "0.35"})
        self.assertEqual(
            [bl.attributes for bl in parsed.bedlines],
            [bl.attributes for bl in scheme.bedlines],
        )

    def test_headers_are_not_carried(self):
        # The delimited format has nowhere to put them
        scheme, parsed = self.round_trip(TEST_ATTRIBUTES_BEDFILE)

        self.assertTrue(scheme.headers)
        self.assertEqual(parsed.headers, [])

    def test_aliased_columns_are_taken_literally(self):
        scheme = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        parsed = Scheme.from_delim_str(scheme.to_delim_str(use_header_aliases=True))

        # "gc" is aliased to "fractiongc" by a header, and without the headers
        # there is no map to reverse it
        self.assertIn("fractiongc", parsed.bedlines[0].attributes)
        self.assertNotIn("gc", parsed.bedlines[0].attributes)

    def test_empty_cells_mean_absent_attribute(self):
        csv_str = (
            "chrom,start,end,primername,pool,strand,sequence,pw\n"
            "chr1,100,120,test_1_LEFT_1,1,+,ACGT,1.4\n"
            "chr1,200,220,test_1_RIGHT_1,1,-,ACGT,\n"
        )
        parsed = Scheme.from_delim_str(csv_str)

        self.assertEqual(parsed.bedlines[0].attributes, {"pw": 1.4})
        self.assertEqual(parsed.bedlines[1].attributes, {})

    def test_blank_and_comment_rows_are_skipped(self):
        csv_str = (
            "# a comment\n"
            "chrom,start,end,primername,pool,strand,sequence\n"
            "\n"
            "chr1,100,120,test_1_LEFT_1,1,+,ACGT\n"
        )
        parsed = Scheme.from_delim_str(csv_str)

        self.assertEqual(len(parsed.bedlines), 1)

    def test_missing_required_column_raises(self):
        csv_str = "chrom,start,end,pool,strand,sequence\nchr1,100,120,1,+,ACGT\n"

        with self.assertRaisesRegex(
            ValueError, r"Missing required column\(s\): primername"
        ):
            Scheme.from_delim_str(csv_str)

    def test_ragged_row_raises(self):
        csv_str = (
            "chrom,start,end,primername,pool,strand,sequence\n"
            "chr1,100,120,test_1_LEFT_1,1,+\n"
        )

        with self.assertRaisesRegex(ValueError, r"Line 2 has 6 field\(s\), expected 7"):
            Scheme.from_delim_str(csv_str)

    def test_derived_column_mismatch_raises(self):
        csv_str = (
            "chrom,start,end,primername,pool,strand,sequence,amplicon_prefix\n"
            "chr1,100,120,test_1_LEFT_1,1,+,ACGT,EDITED\n"
        )

        with self.assertRaisesRegex(ValueError, r"amplicon_prefix \(EDITED\)"):
            Scheme.from_delim_str(csv_str)

    def test_matching_derived_columns_accepted(self):
        csv_str = (
            "chrom,start,end,primername,pool,strand,sequence,"
            "amplicon_prefix,amplicon_number,primer_class_str,primer_suffix\n"
            "chr1,100,120,test_1_LEFT_1,1,+,ACGT,test,1,LEFT,1\n"
        )
        parsed = Scheme.from_delim_str(csv_str)

        self.assertEqual(parsed.bedlines[0].primername, "test_1_LEFT_1")

    def test_empty_input_raises(self):
        with self.assertRaisesRegex(ValueError, r"No rows found"):
            Scheme.from_delim_str("")

    def write_csv(self, text, encoding="utf-8"):
        with tempfile.NamedTemporaryFile(
            "w", suffix=".csv", encoding=encoding, delete=False
        ) as f:
            f.write(text)
        self.addCleanup(os.unlink, f.name)
        return f.name

    def test_bom_is_tolerated(self):
        # Spreadsheets often save as "CSV UTF-8", which writes a leading BOM.
        # It binds to the first column name, hiding "chrom" from the parser.
        csv_str = (
            "chrom,start,end,primername,pool,strand,sequence\n"
            "chr1,100,120,test_1_LEFT_1,1,+,ACGT\n"
        )
        path = self.write_csv(csv_str, encoding="utf-8-sig")

        parsed = Scheme.from_delim_file(path)

        self.assertEqual(parsed.bedlines[0].chrom, "chr1")

    def test_attribute_keyed_like_a_bedline_property(self):
        # Resolving columns via getattr would read the property instead of the
        # stored attribute, silently dropping "weight" and rewriting "length"
        for key in ("weight", "length", "ipool", "amplicon_name"):
            with self.subTest(key=key):
                bedline = BedLine("chr1", 100, 120, "test_1_LEFT_1", 1, "+", "ACGT")
                bedline.attributes = {key: "0.75"}
                csv_str = Scheme(headers=[], bedlines=[bedline]).to_delim_str()

                parsed = Scheme.from_delim_str(csv_str)

                self.assertEqual(parsed.bedlines[0].attributes, {key: "0.75"})

    def test_values_needing_quoting_round_trip(self):
        for value in ("a,b", '"quoted', 'say"hi"'):
            with self.subTest(value=value):
                bedline = BedLine("chr1", 100, 120, "test_1_LEFT_1", 1, "+", "ACGT")
                bedline.attributes = {"note": value}
                csv_str = Scheme(headers=[], bedlines=[bedline]).to_delim_str()

                parsed = Scheme.from_delim_str(csv_str)

                self.assertEqual(parsed.bedlines[0].attributes, {"note": value})

    def test_duplicate_column_raises(self):
        csv_str = (
            "chrom,start,end,primername,pool,strand,sequence,pw,pw\n"
            "chr1,100,120,test_1_LEFT_1,1,+,ACGT,1.4,9.9\n"
        )

        with self.assertRaisesRegex(ValueError, r"Duplicate column name\(s\): pw"):
            Scheme.from_delim_str(csv_str)

    def test_crlf_is_tolerated(self):
        csv_str = (
            "chrom,start,end,primername,pool,strand,sequence,pw\r\n"
            "chr1,100,120,test_1_LEFT_1,1,+,ACGT,1.4\r\n"
        )

        parsed = Scheme.from_delim_str(csv_str)

        self.assertEqual(len(parsed.bedlines), 1)
        self.assertEqual(parsed.bedlines[0].attributes, {"pw": 1.4})
