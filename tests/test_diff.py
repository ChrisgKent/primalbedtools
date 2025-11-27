import unittest

from primalbedtools.diff import (
    diff_primernames,
    diff_sequence,
    ndiff_bedlines,
    unified_diff_bedlines,
)
from primalbedtools.scheme import Scheme
from tests.infiles import (
    TEST_ATTRIBUTES_BEDFILE,
)


class TestNDiff(unittest.TestCase):
    def setUp(self) -> None:
        self.scheme1 = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        self.scheme2 = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        return super().setUp()

    def test_no_change(self):
        # Ignore diff is empty
        self.assertEqual(
            list(
                ndiff_bedlines(
                    self.scheme1.bedlines, self.scheme1.bedlines, ignore_no_diff=True
                )
            ),
            [],
        )

    def test_ignore_header(self):
        """
        With ignore heading flag false detects a difference, with flag True no diff detected
        """
        self.scheme1.headers = ["# header1"]
        self.scheme2.headers = ["# header2"]

        # With ignore_header=False, should see differences
        diffs = list(
            ndiff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                header1=self.scheme1.headers,
                header2=self.scheme2.headers,
                ignore_header=False,
                ignore_no_diff=True,
            )
        )
        self.assertNotEqual(diffs, [])
        self.assertTrue(any(line.startswith("- # header1") for line in diffs))
        self.assertTrue(any(line.startswith("+ # header2") for line in diffs))

        # With ignore_header=True, should see no differences
        diffs = list(
            ndiff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                header1=self.scheme1.headers,
                header2=self.scheme2.headers,
                ignore_header=True,
                ignore_no_diff=True,
            )
        )
        self.assertEqual(diffs, [])

    def test_ignore_order(self):
        """
        With ignore_order flag false detects a difference, with flag True no diff detected
        """
        self.scheme2.bedlines.reverse()

        # With ignore_order=False, should see differences
        diffs = list(
            ndiff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                ignore_order=False,
                ignore_no_diff=True,
            )
        )
        self.assertNotEqual(diffs, [])

        # With ignore_order=True, should see no differences (since content is same)
        diffs = list(
            ndiff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                ignore_order=True,
                ignore_no_diff=True,
            )
        )
        self.assertEqual(diffs, [])

    def test_ignore_attr(self):
        """
        With ignore_attr flag false detects a difference, with flag True no diff detected
        """
        # add a new attribute of 'test=1' to scheme2
        self.scheme2.bedlines[0].attributes = {"test": 1}

        # With ignore_attr=False, should see differences
        diffs = list(
            ndiff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                ignore_attr=False,
                ignore_no_diff=True,
            )
        )
        self.assertNotEqual(diffs, [])

        # With ignore_attr=True, should see no differences
        diffs = list(
            ndiff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                ignore_attr=True,
                ignore_no_diff=True,
            )
        )
        self.assertEqual(diffs, [])

    def test_content_change(self):
        """
        Test that changes in content (e.g. start position) are detected.
        """
        # Change start position of first bedline
        original_start = self.scheme2.bedlines[0].start
        self.scheme2.bedlines[0].start = original_start + 10

        diffs = list(
            ndiff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                ignore_no_diff=True,
            )
        )
        self.assertNotEqual(diffs, [])
        # Should see removal of old line and addition of new line
        self.assertTrue(any(line.startswith("- ") for line in diffs))
        self.assertTrue(any(line.startswith("+ ") for line in diffs))

    def test_empty_inputs(self):
        """
        Test behaviour with empty inputs.
        """
        diffs = list(
            ndiff_bedlines(
                [],
                [],
                ignore_no_diff=True,
            )
        )
        self.assertEqual(diffs, [])

        # One empty
        diffs = list(
            ndiff_bedlines(
                self.scheme1.bedlines,
                [],
                ignore_no_diff=True,
            )
        )
        self.assertNotEqual(diffs, [])
        # Should only have deletions
        self.assertTrue(all(line.startswith("- ") for line in diffs))

    def test_different_lengths(self):
        """
        Test comparing schemes with different numbers of bedlines.
        """
        # Remove the last bedline from scheme2
        self.scheme2.bedlines.pop()

        diffs = list(
            ndiff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                ignore_no_diff=True,
            )
        )
        self.assertNotEqual(diffs, [])
        # Should see the removed line as a deletion
        self.assertTrue(any(line.startswith("- ") for line in diffs))


class TestUnifiedDiff(unittest.TestCase):
    def setUp(self) -> None:
        self.scheme1 = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        self.scheme2 = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        return super().setUp()

    def test_no_change(self):
        self.assertEqual(
            list(unified_diff_bedlines(self.scheme1.bedlines, self.scheme1.bedlines)),
            [],
        )

    def test_ignore_header(self):
        self.scheme1.headers = ["# header1"]
        self.scheme2.headers = ["# header2"]

        # With ignore_header=False
        diffs = list(
            unified_diff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                header1=self.scheme1.headers,
                header2=self.scheme2.headers,
                ignore_header=False,
            )
        )
        self.assertNotEqual(diffs, [])
        self.assertTrue(any(line.startswith("---") for line in diffs))
        self.assertTrue(any(line.startswith("+++") for line in diffs))

        # With ignore_header=True
        diffs = list(
            unified_diff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                header1=self.scheme1.headers,
                header2=self.scheme2.headers,
                ignore_header=True,
            )
        )
        self.assertEqual(diffs, [])

    def test_ignore_order(self):
        self.scheme2.bedlines.reverse()

        # With ignore_order=False
        diffs = list(
            unified_diff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                ignore_order=False,
            )
        )
        self.assertNotEqual(diffs, [])

        # With ignore_order=True
        diffs = list(
            unified_diff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                ignore_order=True,
            )
        )
        self.assertEqual(diffs, [])

    def test_ignore_attr(self):
        self.scheme2.bedlines[0].attributes = {"test": 1}

        # With ignore_attr=False
        diffs = list(
            unified_diff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                ignore_attr=False,
            )
        )
        self.assertNotEqual(diffs, [])

        # With ignore_attr=True
        diffs = list(
            unified_diff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
                ignore_attr=True,
            )
        )
        self.assertEqual(diffs, [])

    def test_content_change(self):
        original_start = self.scheme2.bedlines[0].start
        self.scheme2.bedlines[0].start = original_start + 10

        diffs = list(
            unified_diff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
            )
        )
        self.assertNotEqual(diffs, [])
        self.assertTrue(any(line.startswith("-") for line in diffs))
        self.assertTrue(any(line.startswith("+") for line in diffs))

    def test_empty_inputs(self):
        diffs = list(
            unified_diff_bedlines(
                [],
                [],
            )
        )
        self.assertEqual(diffs, [])

        diffs = list(
            unified_diff_bedlines(
                self.scheme1.bedlines,
                [],
            )
        )
        self.assertNotEqual(diffs, [])

    def test_different_lengths(self):
        self.scheme2.bedlines.pop()

        diffs = list(
            unified_diff_bedlines(
                self.scheme1.bedlines,
                self.scheme2.bedlines,
            )
        )
        self.assertNotEqual(diffs, [])


class TestDiffPrimerNames(unittest.TestCase):
    def setUp(self) -> None:
        self.scheme1 = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        self.scheme2 = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        return super().setUp()

    def test_no_change(self):
        diff1, diff2 = diff_primernames(self.scheme1.bedlines, self.scheme2.bedlines)
        self.assertEqual(diff1, set())
        self.assertEqual(diff2, set())

    def test_diff(self):
        removed_primer = self.scheme2.bedlines.pop(0)
        diff1, diff2 = diff_primernames(self.scheme1.bedlines, self.scheme2.bedlines)
        self.assertIn(removed_primer.primername, diff1)
        self.assertEqual(diff2, set())


class TestDiffSequence(unittest.TestCase):
    def setUp(self) -> None:
        self.scheme1 = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        self.scheme2 = Scheme.from_file(str(TEST_ATTRIBUTES_BEDFILE))
        return super().setUp()

    def test_no_change(self):
        diff1, diff2 = diff_sequence(self.scheme1.bedlines, self.scheme2.bedlines)
        self.assertEqual(diff1, set())
        self.assertEqual(diff2, set())

    def test_diff(self):
        original_seq = self.scheme2.bedlines[0].sequence
        self.scheme2.bedlines[0].sequence = "AAAA"
        diff1, diff2 = diff_sequence(self.scheme1.bedlines, self.scheme2.bedlines)
        self.assertIn(original_seq, diff1)
        self.assertIn("AAAA", diff2)


if __name__ == "__main__":
    unittest.main()
