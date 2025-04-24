import unittest

from primalbedtools.utils import complement_seq, rc_seq


class TestValidate(unittest.TestCase):
    def test_rc_seq(self):
        self.assertEqual("CGAT", rc_seq("ATCG"))

        self.assertEqual(
            "TATTTGGAATCTGAAGGGGACCGGGAATTGG", rc_seq("CCAATTCCCGGTCCCCTTCAGATTCCAAATA")
        )

    def test_complement_seq(self):
        self.assertEqual("CGAT", complement_seq("GCTA"))
