import pathlib
import random
import unittest

from primalbedtools.scheme import DEFAULT_CSV_HEADERS, Scheme

TEST_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.bed"
TEST_V2_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.v2.bed"
TEST_WEIGHTS_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.weights.bed"
TEST_WEIGHTS_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.weights.bed"
TEST_ATTRIBUTES_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.attributes.bed"
TEST_PROBE_BEDFILE = pathlib.Path(__file__).parent / "inputs/test.probe.bed"


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
