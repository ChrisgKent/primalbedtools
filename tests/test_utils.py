import unittest

from primalbedtools.utils import (
    complement_seq,
    expand_ambiguous_bases,
    rc_seq,
    strip_all_white_space,
)


class TestValidate(unittest.TestCase):
    def test_rc_seq(self):
        self.assertEqual("CGAT", rc_seq("ATCG"))

        self.assertEqual(
            "TATTTGGAATCTGAAGGGGACCGGGAATTGG", rc_seq("CCAATTCCCGGTCCCCTTCAGATTCCAAATA")
        )

    def test_complement_seq(self):
        self.assertEqual("CGAT", complement_seq("GCTA"))

    def test_expand_ambiguous_bases(self):
        self.assertEqual(["A", "C"], expand_ambiguous_bases("M"))
        self.assertEqual(["A", "G"], expand_ambiguous_bases("R"))
        self.assertEqual(["A", "T"], expand_ambiguous_bases("W"))
        self.assertEqual(["C", "G"], expand_ambiguous_bases("S"))
        self.assertEqual(["C", "T"], expand_ambiguous_bases("Y"))
        self.assertEqual(["G", "T"], expand_ambiguous_bases("K"))
        self.assertEqual(["A", "C", "G"], expand_ambiguous_bases("V"))
        self.assertEqual(["A", "C", "T"], expand_ambiguous_bases("H"))
        self.assertEqual(["A", "G", "T"], expand_ambiguous_bases("D"))
        self.assertEqual(["C", "G", "T"], expand_ambiguous_bases("B"))
        self.assertEqual(["A", "C", "G", "T"], expand_ambiguous_bases("N"))

        # Test with a sequence containing multiple ambiguous bases
        self.assertEqual(["ACGT"], expand_ambiguous_bases("ACGT"))
        self.assertEqual(
            ["ACGA", "ACGC", "ACGG", "ACGT"], expand_ambiguous_bases("ACGN")
        )

        # Text non dna characters are raised
        with self.assertRaises(KeyError):
            expand_ambiguous_bases("ACGT/3SpC3/")

    def test_strip_all_white_space(self):
        self.assertEqual("", strip_all_white_space("\t"))
        self.assertEqual("", strip_all_white_space(" "))
        self.assertEqual("", strip_all_white_space("    "))
